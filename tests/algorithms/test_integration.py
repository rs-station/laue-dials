import numpy as np
import pytest

from laue_dials.algorithms.integration import (Integrator, cov,
                                               estimate_integration_radius,
                                               mvn_log_pdf,
                                               unmasked_prediction_selection)


# ---------------------------------------------------------------------------
# Synthetic-data helpers
# ---------------------------------------------------------------------------
def make_synthetic_image(
    shape=(200, 200),
    centroids=None,
    amplitude=300.0,
    sigma=2.0,
    background=5.0,
    seed=0,
):
    """
    Build a synthetic detector image containing Gaussian spots on a Poisson
    background.

    Args:
        shape (tuple): (height, width) of the image.
        centroids (np.ndarray): (n, 2) array of spot centers in (x, y) pixel
            coordinates. If None, a 5x5 grid is generated.
        amplitude (float): Peak counts of each Gaussian spot.
        sigma (float): Standard deviation of each Gaussian spot in pixels.
        background (float): Mean background counts.
        seed (int): Seed for the random number generator.

    Returns:
        tuple: (pixels, centroids) where pixels is a float32 image and
            centroids is the (n, 2) array of spot centers.
    """
    rng = np.random.default_rng(seed)
    h, w = shape
    if centroids is None:
        xs = np.linspace(0.15 * w, 0.85 * w, 5)
        ys = np.linspace(0.15 * h, 0.85 * h, 5)
        centroids = np.array([(x, y) for x in xs for y in ys], dtype="float32")
    centroids = np.asarray(centroids, dtype="float32")

    pixels = rng.poisson(background, size=(h, w)).astype("float32")
    yy, xx = np.mgrid[0:h, 0:w]
    for cx, cy in centroids:
        g = amplitude * np.exp(-(((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * sigma**2)))
        pixels += rng.poisson(g).astype("float32")
    return pixels, centroids


# ---------------------------------------------------------------------------
# cov / mvn_log_pdf unit tests
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("ddof", [0, 1])
def test_cov(ddof):
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


def test_cov_invalid_ddof():
    m = np.random.rand(2, 20, 3)
    with pytest.raises(ValueError):
        cov(m, ddof=2)


def test_cov_return_mean():
    b, n, d = 3, 500, 4
    m = np.random.rand(b, n, d)
    S, loc = cov(m, return_mean=True)
    assert S.shape == (b, d, d)
    assert loc.shape == (b, 1, d)
    assert np.allclose(loc.squeeze(-2), m.mean(axis=-2))


def test_mvn_log_pdf():
    d = 10
    b = 5
    m = 50
    X = np.random.random(b * m * d).reshape((b, m, d))

    from scipy.stats import multivariate_normal

    loc = np.random.rand(b * d).reshape(b, d)
    L = np.tril(np.random.rand(b * d * d).reshape((b, d, d))) + np.eye(d)
    S = L.swapaxes(-1, -2) @ L

    expected = np.stack(
        [multivariate_normal.logpdf(i, j, k) for i, j, k in zip(X, loc, S)]
    )
    result = mvn_log_pdf(X, loc, S)
    assert np.allclose(expected, result)


def test_mvn_log_pdf_return_zscore():
    d = 2
    b = 4
    m = 30
    X = np.random.random((b, m, d))
    loc = np.random.rand(b, d)
    L = np.tril(np.random.rand(b, d, d)) + np.eye(d)
    S = L.swapaxes(-1, -2) @ L

    log_p, zscore = mvn_log_pdf(X, loc, S, return_zscore=True)
    assert log_p.shape == (b, m)
    assert zscore.shape == (b, m)
    # The squared Mahalanobis distance is non-negative for a valid covariance.
    assert (zscore >= 0).all()
    # log_p must equal the value returned without the zscore flag.
    assert np.allclose(log_p, mvn_log_pdf(X, loc, S))


# ---------------------------------------------------------------------------
# Integrator geometry / construction unit tests
# ---------------------------------------------------------------------------
def test_integrator_explicit_radius_and_circular_window():
    pixels, centroids = make_synthetic_image()
    radius = 6
    integ = Integrator(pixels, centroids, radius=radius)

    assert integ.radius == radius
    # Window mask should contain only offsets within the circle of `radius`.
    r = np.sqrt(np.square(integ.window_mask).sum(-1))
    assert (r <= radius + 1e-9).all()
    # The mask should be non-trivial (more than a single pixel).
    assert len(integ.window_mask) > 1


def test_integrator_dynamic_radius_is_positive_int():
    pixels, centroids = make_synthetic_image()
    integ = Integrator(pixels, centroids)  # radius=None -> computed
    assert isinstance(integ.radius, (int, np.integer))
    assert integ.radius > 0


def test_integrator_window_indices_clamped_to_bounds():
    """Windows for spots at the image corners must stay within the detector."""
    h, w = 60, 80
    centroids = np.array(
        [[0, 0], [w - 1, h - 1], [0, h - 1], [w - 1, 0]], dtype="float32"
    )
    pixels, centroids = make_synthetic_image(
        shape=(h, w), centroids=centroids, sigma=1.5
    )
    integ = Integrator(pixels, centroids, radius=5)

    # window_idx is [xy, refl, pixel]; row 0 indexes image rows (0..h-1),
    # row 1 indexes image columns (0..w-1).
    assert integ.window_idx[0].min() >= 0
    assert integ.window_idx[0].max() <= h - 1
    assert integ.window_idx[1].min() >= 0
    assert integ.window_idx[1].max() <= w - 1
    # Indexing the pixel array with the clamped windows must not raise.
    assert np.isfinite(integ.windows).all()


# ---------------------------------------------------------------------------
# Integrator execution test
# ---------------------------------------------------------------------------
def test_integrator_fit_execution():
    """
    End-to-end execution of the fit loop on synthetic spots. This exercises
    assign_knn, estimate_background, estimate_profiles, integrate and
    set_strong together (kmdalton: the integrator test should cover execution,
    not just construction).
    """
    pixels, centroids = make_synthetic_image(seed=7)
    integ = Integrator(pixels, centroids)

    integ.fit()

    # All outputs are finite (no NaN/inf leaking from divisions).
    assert np.isfinite(integ.intensity).all()
    assert np.isfinite(integ.uncertainty).all()
    assert np.isfinite(integ.background).all()
    assert np.isfinite(integ.profile_scale).all()
    assert np.isfinite(integ.profile_loc).all()

    # Bright, well-separated spots should be recovered as positive intensities
    # and be flagged strong.
    assert (integ.intensity > 0).all()
    assert integ.uncertainty.shape == integ.intensity.shape
    assert integ.strong.any()


def test_integrator_predict_shapes_and_positivity():
    pixels, centroids = make_synthetic_image(seed=3)
    integ = Integrator(pixels, centroids)
    integ.fit()
    v = integ.predict()
    assert v.shape == integ.windows.shape
    # predict() sums max(0, I) * profile and averages background, both
    # non-negative.
    assert (v > 0).all()


def _reference_predict(integ):
    """Loop-based predict() that sums overlapping windows pixel by pixel."""
    p = integ.profile_values
    signal = np.maximum(0.0, integ.intensity)[:, None] * p
    bg = np.broadcast_to(integ.background, p.shape)
    rows, cols = integ.window_idx
    per_pixel = {}
    for i in range(integ.n):
        seen = set()
        for k in range(integ.m):
            key = (rows[i, k], cols[i, k])
            if key in seen:  # clamped repeat within one window
                continue
            seen.add(key)
            per_pixel.setdefault(key, []).append((signal[i, k], bg[i, k]))
    v = np.empty_like(p)
    for i in range(integ.n):
        for k in range(integ.m):
            contribs = per_pixel[(rows[i, k], cols[i, k])]
            v[i, k] = sum(s for s, _ in contribs) + np.mean([b for _, b in contribs])
    return v


def _randomize_state(integ, seed):
    rng = np.random.default_rng(seed)
    integ.intensity = rng.uniform(-50.0, 500.0, integ.n)
    integ.background = rng.uniform(1.0, 10.0, (integ.n, 1))


def test_predict_matches_reference_with_overlaps_and_edges():
    # Crowded spots, some near the image edge, so windows both overlap one
    # another and get clamped at the detector boundary.
    centroids = np.array(
        [[10, 10], [14, 12], [2, 30], [5, 33], [40, 40], [43, 38], [45, 44], [58, 2]],
        dtype="float32",
    )
    pixels, centroids = make_synthetic_image(shape=(60, 60), centroids=centroids)
    integ = Integrator(pixels, centroids, radius=5)
    _randomize_state(integ, seed=11)
    assert np.allclose(integ.predict(), _reference_predict(integ))


def test_predict_two_overlapping_reflections_matches_formula():
    centroids = np.array([[20, 20], [24, 20]], dtype="float32")
    pixels, centroids = make_synthetic_image(shape=(40, 40), centroids=centroids)
    integ = Integrator(pixels, centroids, radius=4)
    integ.intensity = np.array([100.0, 30.0])
    integ.background = np.array([[2.0], [6.0]])
    p = integ.profile_values
    v = integ.predict()

    # Match pixels shared by the two windows.
    pix0 = {tuple(ij): k for k, ij in enumerate(integ.window_idx[:, 0].T)}
    shared = [
        (pix0[tuple(ij)], k)
        for k, ij in enumerate(integ.window_idx[:, 1].T)
        if tuple(ij) in pix0
    ]
    assert shared
    for k0, k1 in shared:
        expected = 100.0 * p[0, k0] + 30.0 * p[1, k1] + 0.5 * 2.0 + 0.5 * 6.0
        assert np.isclose(v[0, k0], expected)
        assert np.isclose(v[1, k1], expected)

    # Unshared pixels see only their own reflection.
    only0 = np.setdiff1d(np.arange(integ.m), [k0 for k0, _ in shared])
    assert np.allclose(v[0, only0], 100.0 * p[0, only0] + 2.0)


def test_predict_without_overlaps_is_unchanged():
    pixels, centroids = make_synthetic_image(seed=3)
    integ = Integrator(pixels, centroids, radius=5)
    _randomize_state(integ, seed=4)
    expected = (
        np.maximum(0.0, integ.intensity[:, None]) * integ.profile_values
        + integ.background
    )
    assert np.allclose(integ.predict(), expected)


# ---------------------------------------------------------------------------
# Guards flagged in the review
# ---------------------------------------------------------------------------
def test_assign_knn_requires_at_least_two_strong_spots():
    pixels, centroids = make_synthetic_image(seed=5)
    integ = Integrator(pixels, centroids)

    integ.strong = np.zeros(integ.n, dtype=bool)
    with pytest.raises(RuntimeError):
        integ.assign_knn()

    integ.strong = np.zeros(integ.n, dtype=bool)
    integ.strong[0] = True
    with pytest.raises(RuntimeError):
        integ.assign_knn()


def test_assign_knn_reduces_neighbors_when_few_strong_spots():
    pixels, centroids = make_synthetic_image(seed=6)
    integ = Integrator(pixels, centroids, k=5)
    # Only three strong spots: k must fall back to n_strong - 1 = 2.
    integ.strong = np.zeros(integ.n, dtype=bool)
    integ.strong[[0, 1, 2]] = True
    integ.assign_knn()
    assert integ.knn.shape == (integ.n, 2)


def test_estimate_profiles_ignores_negative_weight_contributions():
    """
    Regression test for the profile-weight sign bug kmdalton flagged: when a
    neighbor's (counts - background) is negative *and* its intensity is
    negative, the old ``max(0, (c - bg) / I)`` kept a spurious positive weight.
    The fix clamps the numerator and denominator independently so such a
    neighbor contributes zero weight and does not corrupt the profile estimate.
    """
    pixels, centroids = make_synthetic_image(seed=8)
    integ = Integrator(pixels, centroids)

    # Make exactly two spots strong so the second one's *only* KNN neighbor is
    # the first (pathological) spot.
    integ.strong = np.zeros(integ.n, dtype=bool)
    integ.strong[[0, 1]] = True
    integ.assign_knn()
    assert integ.knn[1].tolist() == [0]

    integ.estimate_background()
    # Poison spot 0: negative intensity and pixels far below background.
    integ.intensity = integ.intensity.astype(float)
    integ.intensity[0] = -100.0
    integ.background = integ.background.copy()
    integ.background[0] = 1e9

    loc_before = integ.profile_loc[1].copy()
    scale_before = integ.profile_scale[1].copy()
    integ.estimate_profiles()

    # Spot 1's profile must be untouched, since its sole neighbor contributes
    # zero weight. No NaNs may appear anywhere.
    assert np.allclose(integ.profile_loc[1], loc_before)
    assert np.allclose(integ.profile_scale[1], scale_before)
    assert np.isfinite(integ.profile_scale).all()


def test_weight_clamping_is_nonnegative_for_negative_over_negative():
    """
    Direct demonstration of the fix: clamping the whole fraction (old code)
    admits a large spurious positive weight, while clamping numerator and
    denominator separately (new code) yields zero.
    """
    c = np.array([1.0])
    bg = np.array([1e9])
    intensity = np.array([-100.0])
    eps = 1e-6

    buggy = np.maximum(0.0, (c - bg) / intensity)
    fixed = np.maximum(0.0, c - bg) / np.maximum(eps, intensity)

    assert buggy[0] > 0  # the bug the fix guards against
    assert fixed[0] == 0


# ---------------------------------------------------------------------------
# Mask-dilation / radius-estimation helpers (used by laue.integrate)
# ---------------------------------------------------------------------------
def test_unmasked_prediction_selection_no_mask_buffers_all_edges():
    """
    With no bad pixels marked (i.e. no external mask file supplied),
    predictions near any detector edge should be excluded by a buffer of
    width `radius`, not just near pixel (0, 0).
    """
    w = h = 100
    np_mask = np.ones((h, w), dtype=bool)  # No bad pixels: simulates no mask file
    radius = 5

    points = {
        "left": (2, h // 2),
        "right": (w - 3, h // 2),
        "top": (w // 2, 2),
        "bottom": (w // 2, h - 3),
        "top_left_corner": (2, 2),
        "bottom_right_corner": (w - 3, h - 3),
        "interior": (w // 2, h // 2),
    }
    x = np.array([p[0] for p in points.values()])
    y = np.array([p[1] for p in points.values()])

    sel = unmasked_prediction_selection(np_mask, x, y, radius, img_row_size=w)
    kept = dict(zip(points.keys(), sel))

    for name, is_kept in kept.items():
        if name == "interior":
            assert is_kept, f"{name} should be kept"
        else:
            assert not is_kept, f"{name} should be removed"


def test_unmasked_prediction_selection_dilates_around_real_mask():
    """
    A real masked (bad) region should still be dilated by the requested
    radius, excluding nearby predictions while keeping distant ones.
    """
    w = h = 100
    np_mask = np.ones((h, w), dtype=bool)
    np_mask[40:60, 40:60] = False  # Bad region in the middle of the detector
    radius = 5

    points = {
        "near_bad_region": (35, 50),  # Exactly `radius` px from the bad region
        "far_from_bad_region": (10, 10),
        "inside_bad_region": (50, 50),
    }
    x = np.array([p[0] for p in points.values()])
    y = np.array([p[1] for p in points.values()])

    sel = unmasked_prediction_selection(np_mask, x, y, radius, img_row_size=w)
    kept = dict(zip(points.keys(), sel))

    assert kept["far_from_bad_region"]
    assert not kept["inside_bad_region"]
    assert not kept["near_bad_region"]


def test_unmasked_prediction_selection_removes_on_mask_centroid_at_zero_radius():
    """
    A centroid sitting on a masked pixel must be removed for any radius,
    including radius == 0 (no dilation): dilation only grows the bad region,
    so on-mask centroids are never spared.
    """
    w = h = 40
    np_mask = np.ones((h, w), dtype=bool)
    np_mask[20, 20] = False  # single bad pixel

    x = np.array([20, 10])
    y = np.array([20, 10])
    sel = unmasked_prediction_selection(np_mask, x, y, radius=0, img_row_size=w)

    assert not sel[0]  # on the masked pixel -> removed
    assert sel[1]  # far away -> kept


def test_estimate_integration_radius_matches_formula():
    """
    estimate_integration_radius should match the nearest-neighbor spacing
    formula used by IntegratorBase's default radius.
    """
    from scipy.spatial.distance import pdist, squareform

    rng = np.random.default_rng(0)
    x = rng.uniform(0, 1000, size=200)
    y = rng.uniform(0, 1000, size=200)

    centroids = np.column_stack([x, y])
    dmat = squareform(pdist(centroids))
    closest_spot_dist = np.sort(dmat, axis=0)[1]
    expected = int(np.round(0.5 * np.percentile(closest_spot_dist, 20)))

    assert estimate_integration_radius(centroids) == expected


# ---------------------------------------------------------------------------
# fit() stopping rule
# ---------------------------------------------------------------------------
def _counting_integrator(pixels, centroids, **kwargs):
    """An Integrator that records how many times the update block ran."""

    class Counting(Integrator):
        def __init__(self, *a, **kw):
            super().__init__(*a, **kw)
            self.n_updates = 0

        def integrate(self):
            super().integrate()
            self.n_updates += 1

    return Counting(pixels, centroids, **kwargs)


def test_fit_applies_exactly_maxiter_updates_when_never_converged():
    """
    With a tolerance that can never be satisfied from below, fit() should run
    the full budget. This pins the loop bound independently of the stopping
    rule, and would catch an off-by-one in how many updates get applied.
    """
    pixels, centroids = make_synthetic_image(seed=7)
    integ = _counting_integrator(pixels, centroids)

    integ.fit(maxiter=4, tol=-np.inf)

    assert integ.n_updates == 4


def test_fit_stops_once_improvement_falls_below_tol():
    """
    A tolerance of 100% relative improvement is unreachable, so the second
    iteration must trip the stopping rule regardless of the data. fit() should
    then stop well short of a large maxiter.
    """
    pixels, centroids = make_synthetic_image(seed=7)
    integ = _counting_integrator(pixels, centroids)

    integ.fit(maxiter=50, tol=1.0)

    # One update to establish a baseline score, one more to compare against it.
    assert integ.n_updates == 2


def test_fit_does_not_apply_an_update_past_the_stopping_point():
    """
    The objective is evaluated after an iteration's updates are applied, so a
    stalled iteration ends the loop immediately.

    On this data the score is minimised after a single update and worsens
    thereafter, so the second iteration must trip the rule and fit() must stop
    holding two updates. The previous implementation scored the state *before*
    updating, which detected the stall one iteration late and returned a third,
    unwanted update.
    """
    pixels, centroids = make_synthetic_image(seed=7)
    integ = _counting_integrator(pixels, centroids)

    integ.fit(maxiter=3)

    assert integ.n_updates == 2


def test_fit_maxiter_one_matches_a_single_manual_update():
    """A one-iteration fit is exactly one pass of the update block."""
    pixels, centroids = make_synthetic_image(seed=7)
    integ = Integrator(pixels, centroids)
    integ.fit(maxiter=1)

    again = Integrator(pixels, centroids)
    again.assign_knn()
    again.estimate_background()
    again.estimate_profiles()
    again.integrate()
    again.set_strong()

    assert np.isclose(integ.score, again.score)
    assert np.allclose(integ.intensity, again.intensity)
