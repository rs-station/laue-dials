"""
This file contains useful classes and functions for profiling and integration
"""

import logging

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.spatial import KDTree

logger = logging.getLogger("laue-dials.algorithms.integration")


class IntegratorBase:
    def __init__(
        self, pixels, centroids, radius=None, k=5, isigi_cutoff=3.0, epsilon=1e-6
    ):
        self.pixels = pixels
        self.epsilon = epsilon
        if radius is None:
            radius = estimate_integration_radius(centroids)
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

        self.window_idx = (
            np.round(self.centroids).astype("int")[:, None, :]
            + self.window_mask[None, :, :]
        )

        # Clamp each axis independently to image bounds
        h, w = self.pixels.shape
        self.window_idx[..., 0] = np.clip(self.window_idx[..., 0], 0, h - 1)
        self.window_idx[..., 1] = np.clip(self.window_idx[..., 1], 0, w - 1)

        # order is [xy, refl, pixel]
        # you can index like self.pixels[tuple(self.window_idx)] -> array[refl, pixel]
        self.window_idx = self.window_idx.transpose(2, 0, 1)
        self.m = self.window_idx.shape[-1]

        # Bookkeeping for summing overlapping windows in predict(). Map every
        # window entry onto a compact index over the distinct pixels touched
        # by any window, so per-pixel sums never need a detector-sized array.
        flat_idx = np.ravel_multi_index(tuple(self.window_idx), (h, w))
        self._pixel_ids, self._pixel_inverse = np.unique(flat_idx, return_inverse=True)
        self._pixel_inverse = self._pixel_inverse.reshape(flat_idx.shape)
        # Clamping repeats edge pixels within a single window. Those repeats
        # share xy and therefore profile value, so only the first occurrence
        # of each pixel per reflection contributes to the per-pixel sums.
        order = np.argsort(flat_idx, axis=-1, kind="stable")
        sorted_idx = np.take_along_axis(flat_idx, order, axis=-1)
        first_sorted = np.ones_like(sorted_idx, dtype=bool)
        first_sorted[:, 1:] = sorted_idx[:, 1:] != sorted_idx[:, :-1]
        self._first_occurrence = np.empty_like(first_sorted)
        np.put_along_axis(self._first_occurrence, order, first_sorted, axis=-1)
        self._n_covering = self.pixel_sum(np.ones(flat_idx.shape))

        self.intensity = self.windows.mean(-1)
        self.uncertainty = np.sqrt(self.windows.mean(-1))

        self.k = k
        self.isigi_cutoff = isigi_cutoff
        self.strong = np.ones(len(self.centroids), dtype=bool)  # Start with all strong

    @property
    def windows(self):
        return self.pixels[tuple(self.window_idx)]

    def fit(self, maxiter=10, tol=1e-3):
        """
        Iteratively refine the background, profiles, and intensities.

        The objective is evaluated after an iteration's updates have been
        applied, so the value compared always describes the current state.
        Fitting stops once an iteration fails to improve the objective by at
        least ``tol`` in relative terms, which covers both a plateau and an
        outright increase.

        Args:
            maxiter (int): Maximum number of iterations to run.
            tol (float): Minimum relative improvement in the objective needed
                to keep iterating. Defaults to 1e-3, i.e. 0.1 percent.
        """
        previous = None
        for _ in range(maxiter):
            self.assign_knn()
            self.estimate_background()
            self.estimate_profiles()
            self.integrate()
            self.set_strong()
            if not self.strong.any():
                raise RuntimeError("No strong spots remaining after integration.")
            score = self.score
            if previous is not None:
                improvement = (previous - score) / max(abs(previous), self.epsilon)
                if improvement < tol:
                    break
            previous = score

    def predict(self):
        """
        Predict the expected counts in every window pixel, accounting for
        overlapping reflections.

        Signal is additive, so a pixel covered by several windows receives the
        profile-weighted intensity of every reflection that covers it. The
        background is a single physical quantity per pixel, so it is the mean
        of the covering reflections' background estimates rather than their
        sum. For reflections i and j sharing pixel n this gives

            v[i, n] = v[j, n] = I[i] p[i, n] + I[j] p[j, n]
                                + 0.5 * bg[i] + 0.5 * bg[j]

        and it reduces to I[i] p[i, n] + bg[i] where windows do not overlap.

        Returns:
            np.ndarray: (n_refls, n_window_pixels) array of expected counts.
        """
        return self.predict_pixels()[self._pixel_inverse]

    def predict_pixels(self, p=None):
        """
        Expected counts at each distinct pixel touched by any window, ordered
        like ``self._pixel_ids``. See ``predict`` for the model.
        """
        if p is None:
            p = self.profile_values
        signal = np.maximum(0.0, self.intensity[:, None]) * p
        return self.pixel_sum(signal) + self.pixel_background()

    def pixel_sum(self, values):
        """
        Sum an (n_refls, n_window_pixels) array onto the distinct pixels
        touched by any window. Pixels repeated within one window by edge
        clamping are counted once per reflection.
        """
        first = self._first_occurrence
        return np.bincount(
            self._pixel_inverse[first],
            weights=values[first],
            minlength=len(self._pixel_ids),
        )

    def pixel_background(self):
        """Mean background of the reflections covering each distinct pixel."""
        bg = np.broadcast_to(self.background, self._pixel_inverse.shape)
        return self.pixel_sum(bg) / self._n_covering

    @property
    def pixel_counts(self):
        """Observed counts at each distinct pixel touched by any window."""
        return self.pixels.ravel()[self._pixel_ids]

    def _profile_cache(self):
        """
        Evaluate the profiles on every window pixel, reusing the previous
        evaluation while the profile parameters are unchanged.

        The parameters are compared by value against copies taken at the last
        evaluation, so the cache stays correct whether they are reassigned or
        edited in place. The comparison is cheap next to the evaluation, which
        is otherwise repeated several times per iteration. The cached arrays
        are read-only so callers cannot corrupt them.

        Returns:
            dict: ``log_p`` and ``mdist`` from ``mvn_log_pdf``, and ``p``,
                the profile normalized over each window.
        """
        from scipy.special import softmax

        cache = getattr(self, "_profiles", None)
        if (
            cache is None
            or not np.array_equal(cache["loc"], self.profile_loc)
            or not np.array_equal(cache["scale"], self.profile_scale)
        ):
            log_p, mdist = mvn_log_pdf(
                self.xy, self.profile_loc, self.profile_scale, return_zscore=True
            )
            cache = {
                "loc": self.profile_loc.copy(),
                "scale": self.profile_scale.copy(),
                "log_p": log_p,
                "mdist": mdist,
                "p": softmax(log_p, axis=-1),
            }
            for name in ("log_p", "mdist", "p"):
                cache[name].setflags(write=False)
            self._profiles = cache
        return cache

    def get_log_p_mdist(self):
        cache = self._profile_cache()
        return cache["log_p"], cache["mdist"]

    @property
    def profile_dist(self):
        return self._profile_cache()["mdist"]

    @property
    def log_profile_values(self):
        return self._profile_cache()["log_p"]

    @property
    def profile_values(self):
        return self._profile_cache()["p"]

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
        np.add.at(im, tuple(self.window_idx), window_values)
        return im

    def plot_image(self, pixels=None, autoscale=True, **kwargs):
        """
        autoscale uses skimage.exposure.adjust_log
        """
        if pixels is None:
            pixels = self.pixels
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
    """
    Profile-fitting integrator that deblends overlapping reflections.

    With the profiles and background held fixed, the expected counts are
    linear in the intensities,

        c[n] ~ sum_j I[j] p[j, n] + bg[n],

    and the intensities are found by weighted least squares with weights
    1 / v[n] from ``predict``. Where windows do not overlap this is the
    familiar profile-fitting estimator,

        I[i] = sum_n p[i, n] (c[n] - bg[n]) / v[n] / sum_n p[i, n]^2 / v[n].

    How the intensities are found is set by ``overlap_method``:

    ``"joint"`` (default) solves the normal equations exactly. They are
    block diagonal over clusters of mutually overlapping reflections, so each
    cluster is solved as a small dense system, and the uncertainty is the
    marginal standard deviation from the inverse normal matrix, which grows
    when reflections are blended.

    ``"legacy"`` applies the single-reflection estimator to every window
    regardless of overlaps. Overlapping reflections then still downweight
    shared pixels through ``predict``, but a neighbor's counts are included
    in the intensity. It is kept for comparison with earlier results.

    Args:
        overlap_method (str): ``"joint"`` or ``"legacy"``.
    """

    overlap_methods = ("joint", "legacy")

    def __init__(self, *args, overlap_method="joint", **kwargs):
        if overlap_method not in self.overlap_methods:
            raise ValueError(
                f"overlap_method must be one of {self.overlap_methods}, "
                f"not {overlap_method!r}"
            )
        self.overlap_method = overlap_method
        super().__init__(*args, **kwargs)

    @property
    def pixel_weights(self):
        """Poisson negative log-likelihood of each distinct window pixel."""
        from scipy.stats import poisson

        return -poisson.logpmf(self.pixel_counts, self.predict_pixels())

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
        knn_idx = KDTree(self.centroids[self.strong]).query(
            self.centroids, k=k_actual + 1
        )[1]
        # Strong spots have themselves as first result (distance 0); non-strong
        # spots are not in the tree so their first result is already a neighbor.
        # Select the right k_actual columns for each case.
        self.knn = np.where(self.strong)[0][
            np.where(self.strong[:, None], knn_idx[:, 1:], knn_idx[:, :k_actual])
        ]

    def estimate_background(self):
        w = self.profile_dist
        p = self.profile_values

        # Remove the signal of every reflection covering each pixel, not just
        # this one, so neighboring spots do not inflate the background.
        signal = self.pixel_sum(self.intensity[:, None] * p)[self._pixel_inverse]
        bg = np.average(self.windows - signal, axis=-1, weights=w, keepdims=True)
        self.background = np.maximum(self.epsilon, bg)

    def estimate_profiles(self):
        """
        Re-estimate each profile from the pixels of its nearest strong spots.

        This is the M-step of expectation maximization for a Poisson mixture.
        The counts in each pixel are shared among the reflections and the
        background covering it in proportion to their predicted contribution,

            w[i, n] = c[n] * max(0, I[i]) p[i, n] / v[n],

        so a neighbor's signal is attributed to the neighbor rather than
        dragging this profile toward it. The profile is then the weighted mean
        and covariance of pixel positions. Weighting by the data alone, rather
        than by the data times the current profile, keeps the covariance from
        shrinking at every iteration.
        """
        xy = self.xy - self.centroids[:, None, :]
        signal = np.maximum(0.0, self.intensity)[:, None] * self.profile_values
        # Clamped repeats within a window are the same pixel; count them once.
        w = self._first_occurrence * self.windows * signal / self.predict()

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
        if self.overlap_method == "joint":
            self._integrate_joint()
        else:
            self._integrate_legacy()

    def _integrate_legacy(self):
        """Fit each reflection on its own window, ignoring its neighbors."""
        c = self.windows
        b = self.background
        v = self.predict()
        p = self.profile_values
        w = p / v / np.sum(np.square(p) / v, axis=-1, keepdims=True)
        I = (c - b) * w
        self.intensity = I.sum(-1)
        SigI = v * w
        SigI = np.sqrt(np.sum(SigI, axis=-1))
        self.uncertainty = SigI

    def _integrate_joint(self):
        """Solve the weighted least-squares normal equations exactly."""
        p = self.profile_values
        first = self._first_occurrence
        rows = self._pixel_inverse[first]
        cols = np.broadcast_to(np.arange(self.n)[:, None], first.shape)[first]
        A = sp.csr_matrix(
            (p[first], (rows, cols)), shape=(len(self._pixel_ids), self.n)
        )

        weight = 1.0 / self.predict_pixels(p)
        N = (A.T @ A.multiply(weight[:, None])).tocoo()
        rhs = A.T @ (weight * (self.pixel_counts - self.pixel_background()))

        self.intensity, variance = solve_block_diagonal(N, rhs, self.epsilon)
        self.uncertainty = np.sqrt(variance)


def solve_block_diagonal(N, rhs, epsilon=1e-6):
    """
    Solve a sparse symmetric positive definite system that splits into many
    small independent blocks, and return the diagonal of its inverse.

    The blocks are the connected components of the sparsity graph of ``N``.
    Blocks of equal size are stacked and solved together with dense batched
    linear algebra. Each diagonal is inflated by a relative ``epsilon`` so
    that reflections with indistinguishable profiles, such as two predictions
    landing on the same pixel, give a finite but very uncertain solution
    rather than a singular matrix.

    Args:
        N (scipy.sparse.coo_matrix): (n, n) symmetric positive definite matrix.
        rhs (np.ndarray): (n,) right-hand side.
        epsilon (float): Relative ridge added to the diagonal.

    Returns:
        tuple: (x, diag_inv), the solution of N x = rhs and diag(N^-1).
    """
    n = N.shape[0]
    _, label = connected_components(N, directed=False)
    size = np.bincount(label)
    # Position of every unknown within its block, in index order.
    order = np.argsort(label, kind="stable")
    pos = np.empty(n, dtype=int)
    pos[order] = np.arange(n) - np.repeat(np.cumsum(size) - size, size)

    x = np.empty(n)
    diag_inv = np.empty(n)
    for s in np.unique(size):
        blocks_of_size = np.flatnonzero(size == s)
        slot = np.full(len(size), -1)
        slot[blocks_of_size] = np.arange(len(blocks_of_size))

        members = np.flatnonzero(slot[label] >= 0)
        stack = np.zeros((len(blocks_of_size), s, s))
        b = np.zeros((len(blocks_of_size), s))
        b[slot[label[members]], pos[members]] = rhs[members]

        entries = slot[label[N.row]] >= 0
        r, c = N.row[entries], N.col[entries]
        np.add.at(stack, (slot[label[r]], pos[r], pos[c]), N.data[entries])

        d = np.arange(s)
        stack[:, d, d] *= 1.0 + epsilon
        inv = np.linalg.inv(stack)
        x[members] = (inv @ b[..., None])[slot[label[members]], pos[members], 0]
        diag_inv[members] = inv[slot[label[members]], pos[members], pos[members]]
    return x, diag_inv


def estimate_integration_radius(centroids):
    """
    Estimate the default integration radius from the spacing of spot centroids.

    The radius is half the 20th percentile of nearest-neighbor centroid
    distances, rounded to the nearest integer. The same radius is used both for
    the integration window and for dilating the detector mask when discarding
    predictions that fall in masked regions, so the two stay consistent.

    Args:
        centroids (np.ndarray): (n, 2) array of centroid pixel coordinates.

    Returns:
        int: Estimated radius in pixels.
    """
    # Query k=2 because the nearest result is the point itself, at distance
    # zero. A KDTree prunes to the nearest neighbor directly, where a full
    # pairwise distance matrix would be O(n^2) in both time and memory --
    # hundreds of MB for a densely predicted image.
    nn_dist = KDTree(centroids).query(centroids, k=2)[0][:, 1]
    radius = 0.5 * np.percentile(nn_dist, 20)
    return int(np.round(radius))


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
