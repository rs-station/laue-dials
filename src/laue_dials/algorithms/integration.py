"""
This file contains useful classes and functions for profiling and integration
"""

import logging

import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist, squareform

logger = logging.getLogger("laue-dials.algorithms.integration")


class IntegratorBase:
    def __init__(self, pixels, centroids, radius=None, k=5, isigi_cutoff=3.0):
        self.pixels = pixels
        if radius is None:
            dmat = squareform(pdist(centroids))
            closest_spot_dist = np.sort(dmat, axis=0)[1]
            radius = 0.5 * np.percentile(closest_spot_dist, 20)
            radius = int(np.round(radius))
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
        self.profile_loc = self.centroids
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

        self.intensity = self.windows.mean(-1)
        self.uncertainty = np.sqrt(self.windows.mean(-1))

        self.k = k
        self.isigi_cutoff = isigi_cutoff
        self.strong = np.ones(len(self.centroids), dtype=bool)  # Start with all strong

    @property
    def windows(self):
        return self.pixels[tuple(self.window_idx)]

    def fit(self, maxiter=2):
        obj = []
        self.assign_knn()
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

    def rectify(self, loc, scale):
        """Cast normal random variable to positive support under Sivia's prior"""
        from scipy.stats import truncnorm

        a = -loc / scale
        b = np.inf
        loc = truncnorm.mean(a, b, loc, scale)
        scale = truncnorm.std(a, b, loc, scale)
        return loc, scale

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

        p = softmax(self.log_profile_values, axis=-1)
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
    @property
    def pixel_weights(self):
        v = self.predict()
        from scipy.stats import poisson

        w = -poisson.logpmf(self.windows, v)
        return w

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
        knn = KDTree(self.centroids[self.strong]).query(
            self.centroids, k=k_actual + 1
        )  # k includes self
        knn = knn[1][:, 1:]  # remove self
        self.knn = np.where(self.strong)[0][knn]

    def estimate_background(self):
        c = self.windows
        w = self.profile_dist
        I = self.intensity

        p = self.profile_values

        bg = np.average(c - I[:, None] * p, axis=-1, weights=w, keepdims=True)
        self.background = np.maximum(1e-6, bg)

    def estimate_profiles(self):
        c = self.windows
        xy = self.xy - self.centroids[:, None, :]
        bg = self.background

        p = np.exp(
            self.log_profile_values
        )  # normalized over all space, not the profile
        # Avoid division by zero; np.maximum below zeros out results for intensity <= 0
        safe_intensity = np.where(self.intensity == 0, 1.0, self.intensity)
        w = np.maximum(0.0, (c - bg) / safe_intensity[:, None]) * p

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
            self.profile_scale[has_signal] = pscale

    def integrate(self):
        c = self.windows
        b = self.background
        v = self.predict()
        p = self.profile_values
        w = p / v / np.sum(np.square(p) / v, axis=-1, keepdims=True)
        I = (c - b) * w
        self.intensity = I.sum(-1)
        SigI = v * w
        SigI = np.sqrt(np.sum(SigI, axis=-1))
        self.intensity, self.uncertainty = self.rectify(self.intensity, SigI)


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


def test_cov(ddof=0):
    d = 10
    n = 1_000
    b = 5

    from scipy.stats import multivariate_normal

    loc = np.random.rand(d)
    L = np.tril(np.random.rand(d * d).reshape((d, d)))
    S = L.T @ L

    m = multivariate_normal.rvs(loc, S, size=(b, n))
    if b == 1:
        m = m[None, ...]

    expected = np.stack([np.cov(i.T, ddof=ddof) for i in m])
    result = cov(m, ddof=ddof)
    assert np.allclose(expected, result)

    aweights = np.random.rand(b * n).reshape((b, n))

    expected = np.stack(
        [np.cov(i.T, aweights=a, ddof=ddof) for i, a in zip(m, aweights)]
    )
    result = cov(m, aweights=aweights[..., None], ddof=ddof)
    assert np.allclose(expected, result)


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


def test_mvn_log_pdf():
    d = 10
    b = 5
    m = 50
    X = np.random.random(b * m * d).reshape((b, m, d))  # test points

    from scipy.stats import multivariate_normal

    loc = np.random.rand(b * d).reshape(b, d)
    L = np.tril(np.random.rand(b * d * d).reshape((b, d, d))) + np.eye(d)
    S = L.swapaxes(-1, -2) @ L

    expected = np.stack(
        [multivariate_normal.logpdf(i, j, k) for i, j, k in zip(X, loc, S)]
    )
    result = mvn_log_pdf(X, loc, S)
    assert np.allclose(expected, result)


if __name__ == "__main__":
    test_mvn_log_pdf()
    test_cov(ddof=0)
    test_cov(ddof=1)
