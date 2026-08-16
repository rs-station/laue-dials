"""
This file contains useful classes and functions for profiling and integration
"""

import logging

import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist, squareform

logger = logging.getLogger("laue-dials.algorithms.integration")


def stack_panels(pixels):
    """
    Normalize panel pixel data into a single (n_panels, rows, cols) array.

    A 2D array is treated as a single panel, so the single-panel call
    ``Integrator(pixels, centroids)`` is unchanged. A sequence of 2D arrays --
    one per panel, as returned by ``imageset.get_raw_data()`` -- is padded to
    the largest panel shape; the padding is marked invalid by
    :func:`panel_validity` and never enters the fit.

    Args:
        pixels (np.ndarray or sequence of np.ndarray): Panel pixel data.

    Returns:
        tuple: ``(stack, shapes)``, the padded (n_panels, rows, cols) array and
        the list of true (rows, cols) shapes.
    """
    if isinstance(pixels, np.ndarray) and pixels.ndim == 2:
        return pixels[None, ...], [pixels.shape]
    if isinstance(pixels, np.ndarray) and pixels.ndim == 3:
        return pixels, [pixels.shape[1:]] * len(pixels)

    arrays = [np.asarray(p) for p in pixels]
    if len(arrays) == 0:
        raise ValueError("no panel pixel data supplied")
    if any(a.ndim != 2 for a in arrays):
        raise ValueError("each panel must be a two-dimensional array of pixels")

    shapes = [a.shape for a in arrays]
    rows = max(s[0] for s in shapes)
    cols = max(s[1] for s in shapes)
    stack = np.zeros((len(arrays), rows, cols), dtype=arrays[0].dtype)
    for i, a in enumerate(arrays):
        stack[i, : a.shape[0], : a.shape[1]] = a
    return stack, shapes


def panel_validity(stack_shape, shapes, panel_masks=None):
    """
    Build the per-pixel validity array for a padded panel stack.

    A pixel is valid if it lies inside its panel's true extent and is not
    masked off by the detector mask.

    Args:
        stack_shape (tuple): Shape of the padded stack, (n_panels, rows, cols).
        shapes (list): True (rows, cols) shape of each panel.
        panel_masks (sequence, optional): One boolean mask per panel, True for
            good pixels. Entries may be None. If omitted, every real pixel is
            valid.

    Returns:
        np.ndarray: Boolean array with the shape of the padded stack.
    """
    n_panels, rows, cols = stack_shape
    valid = np.zeros((n_panels, rows, cols), dtype=bool)
    for i, (pr, pc) in enumerate(shapes):
        valid[i, :pr, :pc] = True

    if panel_masks is not None:
        for i, mask in enumerate(panel_masks):
            if mask is None:
                continue
            pr, pc = shapes[i]
            valid[i, :pr, :pc] &= np.asarray(mask, dtype=bool).reshape(pr, pc)
    return valid


class IntegratorBase:
    def __init__(
        self,
        pixels,
        centroids,
        panel_ids=None,
        panel_masks=None,
        neighbor_coords=None,
        radius=None,
        k=5,
        isigi_cutoff=3.0,
        epsilon=1e-6,
    ):
        """
        Args:
            pixels: Either a 2D pixel array (single panel) or a sequence of 2D
                arrays, one per detector panel.
            centroids (np.ndarray): (n, 2) predicted centroids as (x, y) in the
                pixel coordinates of the panel each reflection lies on.
            panel_ids (np.ndarray, optional): (n,) panel index per reflection.
                Defaults to panel 0 for every reflection.
            panel_masks (sequence, optional): One boolean detector mask per
                panel, True for good pixels. Masked pixels are excluded from
                the fit rather than dropping the reflection.
            neighbor_coords (np.ndarray, optional): (n, 2) or (n, 3) centroid
                coordinates in a frame shared by every panel, in pixel units.
                Used only to estimate the integration radius and to find each
                reflection's nearest strong neighbours. Panel-local pixel
                coordinates cannot serve for either on a multi-panel detector:
                every panel starts again at (0, 0), so the panels pile up on
                top of one another, the nearest-neighbour distances collapse
                and with them the estimated radius. Defaults to ``centroids``,
                which is correct for a single panel.
        """
        self.pixels, self.panel_shapes = stack_panels(pixels)
        self.n_panels = len(self.panel_shapes)
        self.panel_valid = panel_validity(
            self.pixels.shape, self.panel_shapes, panel_masks
        )

        centroids = np.asarray(centroids)
        if panel_ids is None:
            panel_ids = np.zeros(len(centroids), dtype=int)
        panel_ids = np.asarray(panel_ids).astype(int)
        if panel_ids.shape != (len(centroids),):
            raise ValueError("panel_ids must have one entry per centroid")
        if panel_ids.size and (
            panel_ids.min() < 0 or panel_ids.max() >= self.n_panels
        ):
            raise ValueError("panel_ids refer to a panel with no pixel data")
        self.panel_ids = panel_ids

        if neighbor_coords is None:
            neighbor_coords = centroids
        self.neighbor_coords = np.asarray(neighbor_coords, dtype=float)

        self.epsilon = epsilon
        if radius is None:
            radius = estimate_integration_radius(self.neighbor_coords)
        window = np.mgrid[-radius : radius + 1, -radius : radius + 1].reshape((2, -1)).T
        r = np.sqrt(np.square(window[:, 0]) + np.square(window[:, 1]))
        self.radius = radius
        self.window_mask = window[r <= radius]
        self.centroids = centroids[..., ::-1]
        self.n = len(self.centroids)
        self.intensity = None
        self.profile_scale = (
            np.ones(self.n)[:, None, None] * np.eye(2) * self.radius / 2.0
        )
        self.profile_loc = self.centroids.copy()
        self.background = np.ones((self.n, 1))

        window_idx = (
            np.round(self.centroids).astype("int")[:, None, :]
            + self.window_mask[None, :, :]
        )

        # Clamp each axis independently to panel bounds, and remember which
        # window pixels were actually inside. The clamped duplicates keep the
        # indexing legal but must not enter the fit: on a segmented detector a
        # large fraction of windows overhang a panel edge (on 50 px wedges with
        # a 9 px radius, about 38% of them), and folding the edge pixel in
        # several times biases both the background and the profile.
        h, w = self.pixels.shape[1:]
        in_bounds = (
            (window_idx[..., 0] >= 0)
            & (window_idx[..., 0] < h)
            & (window_idx[..., 1] >= 0)
            & (window_idx[..., 1] < w)
        )
        window_idx[..., 0] = np.clip(window_idx[..., 0], 0, h - 1)
        window_idx[..., 1] = np.clip(window_idx[..., 1], 0, w - 1)

        # order is [xy, refl, pixel]
        # you can index like self.pixels[tuple(self.window_idx)] -> array[refl, pixel]
        self.window_idx = window_idx.transpose(2, 0, 1)
        self.window_panel = np.broadcast_to(
            self.panel_ids[:, None], self.window_idx.shape[1:]
        )
        self.window_valid = (
            in_bounds
            & self.panel_valid[
                self.window_panel, self.window_idx[0], self.window_idx[1]
            ]
        )
        self.m = self.window_idx.shape[-1]

        n_valid = self.window_valid.sum(-1)
        if (n_valid == 0).any():
            raise ValueError(
                f"{int((n_valid == 0).sum())} reflection(s) have no usable pixels "
                "in their integration window; drop them before integrating."
            )
        self.n_valid = n_valid

        self.intensity = (self.windows * self.window_valid).sum(-1) / n_valid
        self.uncertainty = np.sqrt(np.maximum(self.intensity, 0.0))

        self.k = k
        self.isigi_cutoff = isigi_cutoff
        self.strong = np.ones(len(self.centroids), dtype=bool)  # Start with all strong

    @property
    def windows(self):
        return self.pixels[self.window_panel, self.window_idx[0], self.window_idx[1]]

    def fit(self, maxiter=2):
        obj = []
        for i in range(maxiter):
            obj.append(self.score)
            self.assign_knn()
            self.estimate_background()
            self.estimate_profiles()
            self.integrate()
            self.set_strong()
            if not self.strong.any():
                raise RuntimeError("No strong spots remaining after integration.")
            if len(obj) == 1:
                continue
            if obj[-1] > obj[-2]:
                break

    def predict(self):
        p = self.profile_values
        v = np.maximum(0.0, self.intensity[:, None]) * p + self.background
        return v

    def get_log_p_mdist(self):
        return mvn_log_pdf(
            self.xy, self.profile_loc, self.profile_scale, return_zscore=True
        )

    @property
    def profile_dist(self):
        log_p, mdist = self.get_log_p_mdist()
        return mdist

    @property
    def log_profile_values(self):
        log_p, _ = self.get_log_p_mdist()
        return log_p

    @property
    def profile_values(self):
        from scipy.special import softmax

        # Normalizing over the valid pixels only keeps sum(p) == 1 over the
        # pixels that are actually used, which is what makes the profile-fitted
        # intensity and its variance consistent, and gives p == 0 on every
        # invalid pixel so they drop out of the weighted sums downstream.
        log_p = np.where(self.window_valid, self.log_profile_values, -np.inf)
        p = softmax(log_p, axis=-1)
        return p

    @property
    def xy(self):
        retval = self.window_idx.transpose(1, 2, 0) + 0.5
        return retval

    def plot_profiles(
        self,
        ax=None,
        n_std=2.0,
        weak_color="w",
        strong_color="y",
        facecolor="none",
        **kwargs,
    ):
        import matplotlib.transforms as transforms
        from matplotlib import pyplot as plt
        from matplotlib.patches import Ellipse

        if ax is None:
            ax = plt.gca()

        retval = []
        for loc, cov, strong in zip(self.profile_loc, self.profile_scale, self.strong):
            ecolor = strong_color if strong else weak_color

            loc = loc[..., ::-1]
            cov = cov.swapaxes(-1, -2)
            pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
            # Using a special case to obtain the eigenvalues of this
            # two-dimensional dataset.
            ell_radius_x = np.sqrt(1 + pearson)
            ell_radius_y = np.sqrt(1 - pearson)
            ellipse = Ellipse(
                (0, 0),
                width=ell_radius_x * 2,
                height=ell_radius_y * 2,
                facecolor=facecolor,
                edgecolor=ecolor,
                **kwargs,
            )

            # Calculating the standard deviation of x from
            # the squareroot of the variance and multiplying
            # with the given number of standard deviations.
            scale_x = np.sqrt(cov[0, 0]) * n_std
            scale_y = np.sqrt(cov[1, 1]) * n_std
            mean_x, mean_y = loc

            transf = (
                transforms.Affine2D()
                .rotate_deg(45)
                .scale(scale_x, scale_y)
                .translate(mean_x, mean_y)
            )

            ellipse.set_transform(transf + ax.transData)
            retval.append(ax.add_patch(ellipse))
        return retval

    def single_profile_image(self, fg_values, fill_value=0.0):
        im = (
            np.ones((2 * self.radius + 1, 2 * self.radius + 1), dtype=fg_values.dtype)
            * fill_value
        )
        np.add.at(im, tuple(self.window_mask.T), fg_values)
        im = np.fft.fftshift(im)
        return im

    def fg_to_image(self, window_values, fill_value=0.0):
        """
        convert an array of the same shape as self.pixels[*self.window_idx] -> array(refls x pixels)
        to something the same shape as refls.pixels
        """
        im = np.ones_like(self.pixels) * fill_value
        np.add.at(
            im,
            (self.window_panel, self.window_idx[0], self.window_idx[1]),
            window_values,
        )
        return im

    def plot_image(self, pixels=None, autoscale=True, **kwargs):
        """
        autoscale uses skimage.exposure.adjust_log
        """
        if pixels is None:
            if self.n_panels > 1:
                raise NotImplementedError(
                    "plot_image draws a single panel; pass pixels=... to choose one"
                )
            pixels = self.pixels[0]
        from matplotlib import pyplot as plt

        if autoscale:
            from skimage import exposure

            pixels = exposure.adjust_log(pixels)
        plt.matshow(pixels, **kwargs)

    def plot_image_with_profiles(self):
        self.plot_image()
        for n_std in (1.0, 2.0, 3.0):
            self.plot_profiles(n_std=n_std)


class Integrator(IntegratorBase):
    @property
    def pixel_weights(self):
        v = self.predict()
        from scipy.stats import poisson

        w = -poisson.logpmf(self.windows, v)
        return np.where(self.window_valid, w, 0.0)

    @property
    def score(self):
        """Loss function value"""
        return self.pixel_weights.sum()

    def set_strong(self):
        """set self.strong"""
        self.strong = self.intensity >= self.isigi_cutoff * self.uncertainty

    def assign_knn(self):
        k = self.k
        n_strong = self.strong.sum()
        if n_strong <= 1:
            raise RuntimeError(
                f"Only {n_strong} strong spot(s) found; cannot perform KNN profile estimation."
            )
        k_actual = min(k, n_strong - 1)
        if k_actual < k:
            logger.warning(
                "Only %d strong spots available; using %d neighbors instead of %d.",
                n_strong,
                k_actual,
                k,
            )
        # Neighbours are found in the shared frame, so a reflection near a panel
        # edge pools with the spots physically next to it on the adjacent panel
        # rather than with whatever happens to sit at the same panel-local
        # coordinates several panels away.
        knn_idx = KDTree(self.neighbor_coords[self.strong]).query(
            self.neighbor_coords, k=k_actual + 1
        )[1]
        # Strong spots have themselves as first result (distance 0); non-strong
        # spots are not in the tree so their first result is already a neighbor.
        # Select the right k_actual columns for each case.
        self.knn = np.where(self.strong)[0][
            np.where(self.strong[:, None], knn_idx[:, 1:], knn_idx[:, :k_actual])
        ]

    def estimate_background(self):
        c = self.windows
        w = np.where(self.window_valid, self.profile_dist + self.epsilon, 0.0)
        I = self.intensity

        p = self.profile_values

        bg = np.average(c - I[:, None] * p, axis=-1, weights=w, keepdims=True)
        self.background = np.maximum(self.epsilon, bg)

    def estimate_profiles(self):
        c = self.windows
        xy = self.xy - self.centroids[:, None, :]
        bg = self.background

        p = np.exp(
            self.log_profile_values
        )  # normalized over all space, not the profile
        p = np.where(self.window_valid, p, 0.0)
        num = np.maximum(0.0, (c - bg))
        den = np.maximum(self.epsilon, self.intensity)
        w = num / den[:, None] * p

        w = w[self.knn].reshape((self.n, -1))
        xy = xy[self.knn].reshape((self.n, -1, 2))

        # Only update profiles for reflections with nonzero weights. All-zero
        # weights (e.g. from dead pixels or negative intensity) would cause
        # cov() to divide by zero. Keeping the previous profile estimate is
        # better than overwriting it with nan.
        has_signal = w.sum(-1) > 0
        if has_signal.any():
            pscale, ploc = cov(
                xy[has_signal], w[has_signal][..., None], return_mean=True
            )
            ploc = ploc.squeeze(-2)
            self.profile_loc[has_signal] = ploc + self.centroids[has_signal]
            # Tikhonov regularization: add a small multiple of the identity to
            # guarantee positive definiteness and prevent LinAlgError in
            # mvn_log_pdf when neighboring pixels are collinear.
            self.profile_scale[has_signal] = pscale + self.epsilon * np.eye(2)

    def integrate(self):
        c = self.windows
        b = self.background
        v = self.predict()
        p = self.profile_values
        # p is zero on invalid pixels, so they contribute nothing to either sum
        # and the weights stay normalized over the pixels that are used.
        w = p / v / np.sum(np.square(p) / v, axis=-1, keepdims=True)
        I = (c - b) * w
        self.intensity = I.sum(-1)
        SigI = v * w
        SigI = np.sqrt(np.sum(SigI, axis=-1))
        self.uncertainty = SigI


def estimate_integration_radius(centroids):
    """
    Estimate the default integration radius from the spacing of spot centroids.

    The radius is half the 20th percentile of nearest-neighbor centroid
    distances, rounded to the nearest integer. The same radius is used both for
    the integration window and for dilating the detector mask when discarding
    predictions that fall in masked regions, so the two stay consistent.

    The centroids must be in a frame shared by every panel. Panel-local pixel
    coordinates superimpose the panels, which drives the nearest-neighbor
    distances -- and the radius with them -- towards zero.

    Args:
        centroids (np.ndarray): (n, 2) array of centroid pixel coordinates.

    Returns:
        int: Estimated radius in pixels.
    """
    dmat = squareform(pdist(centroids))
    closest_spot_dist = np.sort(dmat, axis=0)[1]
    radius = 0.5 * np.percentile(closest_spot_dist, 20)
    return int(np.round(radius))


def detector_global_pixels(detector, panel_ids, spots):
    """
    Map panel-local centroids onto a single detector-wide pixel grid.

    ``xyzcal.px`` is panel-local: every panel starts again at (0, 0), so on a
    multi-panel detector the panels are superimposed and the coordinates are
    useless as scaling metadata -- on the LADI drum they pile 48 wedges on top
    of one another. This lays the panels out on a common grid instead.

    The panels are assumed to share a slow direction, which is what makes a
    two-dimensional layout meaningful at all. Writing ``s`` for the mean slow
    axis, each lab-frame point is split into its component along ``s`` (the
    slow coordinate) and its azimuth about ``s`` (the fast coordinate,
    converted to a distance with the mean panel radius). For a curved detector
    that is the unrolled arc length; for a flat one it is a smooth monotonic
    function of the fast coordinate. The branch cut is placed in the largest
    angular gap between panels, so a detector wrapping past 180 degrees -- the
    LADI drum covers about 304 -- does not wrap around on itself. The origin is
    the corner of the panel at the low end of both coordinates, so the result
    is non-negative and depends only on the detector model.

    Args:
        detector: dxtbx detector model, or any sequence of panels supporting
            get_slow_axis, get_origin, get_pixel_size, get_image_size and
            get_pixel_lab_coord.
        panel_ids (np.ndarray): (n,) panel index per centroid.
        spots (np.ndarray): (n, 2) panel-local centroids in pixels.

    Returns:
        np.ndarray or None: (n, 2) centroids on the detector-wide grid, in
        pixels. A single-panel detector is returned unchanged. None if the
        panels have no common slow direction, or if they close a full circle
        and so leave no gap to cut at.
    """
    spots = np.asarray(spots, dtype=float)[:, :2]
    panel_ids = np.asarray(panel_ids).astype(int)
    if len(detector) < 2:
        return spots.copy()

    slow = np.array([p.get_slow_axis() for p in detector], dtype=float)
    s = slow.mean(axis=0)
    norm = np.linalg.norm(s)
    if norm < 1e-9:
        return None
    s = s / norm

    px = np.array([p.get_pixel_size() for p in detector], dtype=float)
    qx, qy = float(px[:, 0].mean()), float(px[:, 1].mean())
    size = np.array([p.get_image_size() for p in detector], dtype=float)

    centres = np.array(
        [p.get_pixel_lab_coord((w / 2.0, h / 2.0)) for p, (w, h) in zip(detector, size)],
        dtype=float,
    )
    cq = centres - np.outer(centres @ s, s)
    radii = np.linalg.norm(cq, axis=1)
    radius = float(radii.mean())
    if radius < 1e-9:
        return None

    e1 = cq[0] / radii[0]
    e2 = np.cross(s, e1)
    e2 = e2 / np.linalg.norm(e2)
    # Orient the azimuth so that it increases along the panels' fast axis, i.e.
    # in the same direction as the panel-local x of xyzcal.px. Without this the
    # sign is whichever way round np.cross happens to come out, and the global
    # coordinate can run backwards against the local one.
    fast = np.array(detector[0].get_fast_axis(), dtype=float)
    if np.dot(fast - np.dot(fast, s) * s, e2) < 0:
        e2 = -e2

    theta_panel = np.arctan2(cq @ e2, cq @ e1)

    # Put the branch cut in the widest angular gap between panels, so the
    # occupied arc is contiguous however far round it goes.
    order = np.argsort(theta_panel)
    ordered = theta_panel[order]
    gaps = np.diff(np.append(ordered, ordered[0] + 2 * np.pi))
    widest = int(np.argmax(gaps))
    if gaps[widest] <= 0:
        return None
    cut = ordered[widest] + 0.5 * gaps[widest]

    lab = np.array(
        [
            detector[int(p)].get_pixel_lab_coord((float(x), float(y)))
            for p, (x, y) in zip(panel_ids, spots)
        ],
        dtype=float,
    )
    along = lab @ s
    q = lab - np.outer(along, s)
    theta = np.mod(np.arctan2(q @ e2, q @ e1) - cut, 2 * np.pi)

    # Anchor on the detector, not on the reflections, so repeated runs and
    # different images of the same detector share one coordinate system.
    theta_ref = np.mod(theta_panel - cut, 2 * np.pi) - 0.5 * size[:, 0] * qx / radius
    along_ref = np.array([np.dot(p.get_origin(), s) for p in detector], dtype=float)

    x = radius * (theta - theta_ref.min()) / qx
    y = (along - along_ref.min()) / qy
    return np.column_stack([x, y])


def unmasked_prediction_selection(np_mask, x, y, radius, img_row_size):
    """
    Determine which predicted centroids fall outside the detector mask,
    dilated by the given radius.

    If np_mask has no bad (False) pixels at all -- e.g. no external mask
    file was supplied -- a genuine 1px border is marked as invalid around
    the detector edge, and the dilation radius is reduced by 1 to
    compensate for that added border. This is needed because
    skimage.morphology.isotropic_dilation relies on
    scipy.ndimage.distance_transform_edt, which has no real background
    reference point when there are no bad pixels at all, and would
    otherwise spuriously mask a small region near pixel (0, 0).

    Args:
        np_mask (np.ndarray): Boolean detector mask with shape (n_rows,
            n_cols); True for valid pixels.
        x (np.ndarray): Integer pixel x-coordinates of predicted centroids.
        y (np.ndarray): Integer pixel y-coordinates of predicted centroids.
        radius (int): Radius in pixels to dilate the detector mask by.
        img_row_size (int): Number of pixels per detector row, used to
            flatten (x, y) coordinates into the flattened mask.

    Returns:
        np.ndarray: Boolean array, True for centroids to keep.
    """
    from skimage.morphology import isotropic_dilation

    bad_pixels = ~np_mask

    dilation_radius = radius
    if not bad_pixels.any():
        bad_pixels[0, :] = True
        bad_pixels[-1, :] = True
        bad_pixels[:, 0] = True
        bad_pixels[:, -1] = True
        dilation_radius = max(radius - 1, 0)

    expanded_mask = ~isotropic_dilation(bad_pixels, dilation_radius)
    expanded_mask_flat = expanded_mask.flatten()
    return expanded_mask_flat[x + img_row_size * y]


def cov(m, aweights=None, return_mean=False, ddof=0):
    """
    A batched version of np.cov to estimate the sample covariance matrix where m has leading batch dims.
    """
    if ddof not in (0, 1):
        raise ValueError(f"ddof can only be 0 or 1, but received {ddof}")

    if aweights is None:
        loc = np.mean(m, axis=-2, keepdims=True)
        if ddof == 0:
            denom = m.shape[-2]
        elif ddof == 1:
            denom = m.shape[-2] - 1
    else:
        loc, w_sum = np.average(
            m, axis=-2, weights=aweights * np.ones_like(m), keepdims=True, returned=True
        )
        if ddof == 0:
            denom = w_sum
        elif ddof == 1:
            denom = (
                w_sum - np.square(aweights).sum(-2, keepdims=True) / w_sum
            )  # This is for ddof=1 version

    X = m - loc
    if aweights is not None:
        X_T = (X * aweights).swapaxes(-1, -2)
    else:
        X_T = X.swapaxes(-1, -2)

    S = X_T @ X / denom
    if return_mean:
        return S, loc
    return S


def mvn_log_pdf(x, loc, scale, return_zscore=False):
    """a batched version of scipy.stats.multivariate_normal.log_pdf"""
    d = loc.shape[-1]
    diff = x - loc[..., None, :]
    Sinv = np.linalg.inv(scale)

    log_Z = -0.5 * d * np.log(2 * np.pi) - 0.5 * np.linalg.slogdet(scale)[1]
    zscore = (diff[..., None, :] @ Sinv[..., None, :, :] @ diff[..., :, None]).squeeze(
        (-1, -2)
    )
    log_p = log_Z[..., None] - 0.5 * zscore
    if return_zscore:
        return log_p, zscore
    return log_p
