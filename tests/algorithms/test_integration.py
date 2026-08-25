import numpy as np
import pytest

from laue_dials.algorithms.integration import (
    Integrator,
    cov,
    detector_global_pixels,
    estimate_integration_radius,
    mvn_log_pdf,
    unmasked_prediction_selection,
)


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
    # predict() = max(0, I) * profile + background, both non-negative.
    assert (v > 0).all()


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
# Multi-panel detectors
# ---------------------------------------------------------------------------
def make_panelled_image(
    n_panels=6,
    panel_width=50,
    panel_height=300,
    n_per_panel=10,
    sigma=1.8,
    background=100.0,
    seed=0,
    narrow_last=True,
):
    """
    Build a set of narrow panels, each carrying Gaussian spots of known total
    intensity, in the style of a segmented detector such as the LADI drum.

    Returns:
        tuple: (panels, centroids, panel_ids, shared, truth) where centroids
            are panel-local (x, y), shared is the same points in a frame common
            to every panel, and truth is the integrated intensity of each spot.
    """
    rng = np.random.default_rng(seed)
    widths = [panel_width] * n_panels
    if narrow_last:
        widths[-1] = panel_width - 1

    panels = [np.full((panel_height, w), background) for w in widths]
    xs, ys, pids, truth = [], [], [], []
    for p, w in enumerate(widths):
        gy, gx = np.mgrid[0:panel_height, 0:w]
        # deliberately include both panel edges
        x = np.linspace(1.5, w - 1.5, n_per_panel)
        y = np.sort(rng.uniform(20, panel_height - 20, n_per_panel))
        intensity = rng.uniform(200.0, 4000.0, n_per_panel)
        for xi, yi, Ii in zip(x, y, intensity):
            panels[p] += (
                Ii
                * np.exp(-((gx - xi) ** 2 + (gy - yi) ** 2) / (2 * sigma**2))
                / (2 * np.pi * sigma**2)
            )
            xs.append(xi)
            ys.append(yi)
            pids.append(p)
            truth.append(Ii)

    panels = [rng.poisson(p).astype("float32") for p in panels]
    centroids = np.stack([xs, ys], axis=-1).astype("float32")
    pids = np.array(pids)
    shared = centroids.astype(float).copy()
    shared[:, 0] += panel_width * pids
    return panels, centroids, pids, shared, np.array(truth)


def test_panel_local_coordinates_collapse_the_radius():
    """
    Panel-local centroids superimpose every panel at the same coordinates, so
    the nearest-neighbor spacing -- and the radius derived from it -- collapses.
    This is why the caller has to supply coordinates in a shared frame.
    """
    _, centroids, _, shared, _ = make_panelled_image()

    assert estimate_integration_radius(centroids) < estimate_integration_radius(shared)


def test_integrator_reads_each_reflection_from_its_own_panel():
    """
    With panel_ids supplied, intensities must track the truth. Integrating
    everything against panel 0 -- what the code did before -- must not.
    """
    panels, centroids, pids, shared, truth = make_panelled_image()

    integ = Integrator(
        panels,
        centroids,
        panel_ids=pids,
        neighbor_coords=shared,
        radius=6,
        k=5,
        isigi_cutoff=1.0,
    )
    integ.fit(maxiter=4)
    assert np.corrcoef(integ.intensity, truth)[0, 1] > 0.9

    panel_zero = Integrator(panels[0], centroids, radius=6, k=5, isigi_cutoff=1.0)
    panel_zero.fit(maxiter=4)
    assert abs(np.corrcoef(panel_zero.intensity, truth)[0, 1]) < 0.5


def test_pixels_beyond_the_panel_edge_are_excluded_not_clamped():
    """
    A window that overhangs a panel edge must not fold the edge pixel in
    several times: those positions are clamped so the indexing stays legal, but
    they carry no weight.
    """
    panels, centroids, pids, shared, truth = make_panelled_image()
    radius = 6

    kwargs = dict(
        panel_ids=pids, neighbor_coords=shared, radius=radius, k=5, isigi_cutoff=1.0
    )
    integ = Integrator(panels, centroids, **kwargs)
    clamped = Integrator(panels, centroids, **kwargs)
    clamped.window_valid = np.ones_like(clamped.window_valid)

    integ.fit(maxiter=4)
    clamped.fit(maxiter=4)

    width = panels[0].shape[1]
    edge = (centroids[:, 0] < radius) | (centroids[:, 0] > width - 1 - radius)
    assert edge.any()

    r_excluded = np.corrcoef(integ.intensity[edge], truth[edge])[0, 1]
    r_clamped = np.corrcoef(clamped.intensity[edge], truth[edge])[0, 1]
    assert r_excluded > r_clamped

    # every window pixel outside the panel must be marked invalid
    assert not integ.window_valid[edge].all()


def test_narrower_panel_padding_is_never_valid():
    """Panels of unequal width are padded to a common shape; the padding is
    not a pixel and must never be used."""
    panels, centroids, pids, shared, truth = make_panelled_image()
    integ = Integrator(
        panels, centroids, panel_ids=pids, neighbor_coords=shared, radius=6
    )

    last = len(panels) - 1
    assert panels[last].shape[1] < panels[0].shape[1]
    assert not integ.panel_valid[last, :, panels[last].shape[1] :].any()


def test_panel_masks_exclude_pixels_without_dropping_the_reflection():
    """A bad-pixel block inside a window costs that window those pixels; the
    reflection itself survives."""
    panels, centroids, pids, shared, truth = make_panelled_image()
    masks = [np.ones(p.shape, dtype=bool) for p in panels]

    i0 = int(np.where(pids == 2)[0][0])
    x, y = int(centroids[i0, 0]), int(centroids[i0, 1])
    masks[2][y - 1 : y + 2, x + 2 : x + 5] = False

    kwargs = dict(
        panel_ids=pids, neighbor_coords=shared, radius=6, k=5, isigi_cutoff=1.0
    )
    plain = Integrator(panels, centroids, **kwargs)
    masked = Integrator(panels, centroids, panel_masks=masks, **kwargs)

    assert masked.window_valid[i0].sum() < plain.window_valid[i0].sum()
    assert masked.n == plain.n

    masked.fit(maxiter=4)
    assert np.isfinite(masked.intensity).all()


def test_single_panel_path_is_unchanged_by_the_panel_machinery():
    """
    A 2D pixel array must behave exactly as before: one panel, every window
    pixel valid for interior spots, window_idx still (2, n, m) and in range.
    """
    pixels, centroids = make_synthetic_image(seed=4)
    integ = Integrator(pixels, centroids, radius=5)

    assert integ.n_panels == 1
    assert integ.window_idx.shape == (2, len(centroids), integ.m)
    assert integ.window_idx[0].min() >= 0
    assert integ.window_idx[0].max() <= pixels.shape[0] - 1
    assert integ.window_idx[1].min() >= 0
    assert integ.window_idx[1].max() <= pixels.shape[1] - 1
    assert (integ.panel_ids == 0).all()

    integ.fit()
    assert np.isfinite(integ.intensity).all()


def test_panel_ids_must_match_the_pixel_data():
    pixels, centroids = make_synthetic_image(seed=4)
    with pytest.raises(ValueError):
        Integrator(pixels, centroids, panel_ids=np.ones(len(centroids), dtype=int))
    with pytest.raises(ValueError):
        Integrator(pixels, centroids, panel_ids=np.zeros(len(centroids) + 1, dtype=int))


# ---------------------------------------------------------------------------
# Detector-wide centroid grid (xcal / ycal in the MTZ)
# ---------------------------------------------------------------------------
class FakePanel:
    """The handful of dxtbx Panel methods detector_global_pixels needs."""

    def __init__(
        self, fast, slow, origin, pixel_size=(0.25, 0.25), image_size=(50, 1800)
    ):
        self._fast = np.asarray(fast, dtype=float)
        self._slow = np.asarray(slow, dtype=float)
        self._origin = np.asarray(origin, dtype=float)
        self._pixel_size = tuple(pixel_size)
        self._image_size = tuple(image_size)

    def get_fast_axis(self):
        return tuple(self._fast)

    def get_slow_axis(self):
        return tuple(self._slow)

    def get_origin(self):
        return tuple(self._origin)

    def get_pixel_size(self):
        return self._pixel_size

    def get_image_size(self):
        return self._image_size

    def get_pixel_lab_coord(self, xy):
        return tuple(
            self._origin
            + xy[0] * self._pixel_size[0] * self._fast
            + xy[1] * self._pixel_size[1] * self._slow
        )


def make_drum(n_panels=60, width=50, height=1800, radius=200.0, pixel=0.25):
    """Flat wedges tangent to a cylinder, as FormatLADI models the LADI plate."""
    panels = []
    for i in range(n_panels):
        c0 = i * width
        gamma = (c0 + 0.5 * width) * pixel / radius
        dgamma = width * pixel / radius
        r_panel = 0.5 * radius * (1.0 + np.cos(0.5 * dgamma))
        fast = (np.cos(gamma), 0.0, -np.sin(gamma))
        slow = (0.0, -1.0, 0.0)
        centre = (r_panel * np.sin(gamma), 0.0, r_panel * np.cos(gamma))
        half = 0.5 * width * pixel
        origin = (
            centre[0] - half * fast[0],
            height * pixel,
            centre[2] - half * fast[2],
        )
        panels.append(FakePanel(fast, slow, origin, (pixel, pixel), (width, height)))
    return panels


def make_tiled(n_fast=4, n_slow=3, width=100, height=80, pixel=0.1, distance=150.0):
    """A flat detector tiled into modules, as most pixel-array detectors are."""
    panels = []
    for j in range(n_slow):
        for i in range(n_fast):
            origin = (i * width * pixel, -j * height * pixel, distance)
            panels.append(
                FakePanel(
                    (1, 0, 0), (0, -1, 0), origin, (pixel, pixel), (width, height)
                )
            )
    return panels


def test_global_pixels_single_panel_is_unchanged():
    panels = make_drum(n_panels=1)
    spots = np.array([[1.0, 2.0], [30.0, 900.0]])
    out = detector_global_pixels(panels, np.zeros(2, dtype=int), spots)
    assert np.allclose(out, spots)


def test_global_pixels_unrolls_a_drum_onto_plate_coordinates():
    """
    On a segmented cylinder the global x must be the unrolled plate column and
    the global y the panel row, so that neighbouring wedges are adjacent
    instead of superimposed.
    """
    width, n_panels = 50, 60
    panels = make_drum(n_panels=n_panels, width=width)
    rng = np.random.default_rng(0)
    panel_ids = rng.integers(0, n_panels, 500)
    spots = np.column_stack([rng.uniform(0, width, 500), rng.uniform(0, 1800, 500)])

    out = detector_global_pixels(panels, panel_ids, spots)
    unrolled = panel_ids * width + spots[:, 0]

    slope, intercept = np.polyfit(unrolled, out[:, 0], 1)
    assert abs(slope - 1.0) < 1e-3  # chord vs arc over one wedge
    assert abs(intercept) < 1e-2
    assert np.abs(out[:, 0] - np.polyval([slope, intercept], unrolled)).max() < 0.05
    assert np.allclose(out[:, 1], spots[:, 1], atol=1e-6)


def test_global_pixels_increase_with_panel_local_x():
    """The global coordinate must run the same way as xyzcal.px, not backwards."""
    width = 50
    panels = make_drum(n_panels=20, width=width)
    spots = np.array([[1.0, 100.0], [width - 1.0, 100.0]])
    out = detector_global_pixels(panels, np.array([7, 7]), spots)
    assert out[1, 0] > out[0, 0]


def test_global_pixels_do_not_wrap_past_half_a_turn():
    """
    The LADI drum covers about 304 degrees. A branch cut left at an arbitrary
    azimuth would fold the far end of the plate back onto the near end.
    """
    width, n_panels = 50, 85  # 85 * 50 * 0.25 mm / 200 mm = 5.3 rad
    panels = make_drum(n_panels=n_panels, width=width)
    panel_ids = np.arange(n_panels)
    spots = np.column_stack([np.full(n_panels, width / 2.0), np.full(n_panels, 900.0)])

    x = detector_global_pixels(panels, panel_ids, spots)[:, 0]
    assert (np.diff(x) > 0).all()  # strictly increasing, no fold
    assert x.max() - x.min() > 0.9 * n_panels * width


def test_global_pixels_lay_out_a_tiled_flat_detector():
    """Modules of a flat detector must separate in both directions."""
    n_fast, n_slow, width, height = 4, 3, 100, 80
    panels = make_tiled(n_fast, n_slow, width, height)
    panel_ids = np.arange(len(panels))
    spots = np.column_stack(
        [np.full(len(panels), width / 2.0), np.full(len(panels), height / 2.0)]
    )
    out = detector_global_pixels(panels, panel_ids, spots)

    x = out[:, 0].reshape(n_slow, n_fast)
    y = out[:, 1].reshape(n_slow, n_fast)
    # x separates the columns of modules and is the same down each column
    assert (np.diff(x, axis=1) > 0).all()
    assert np.allclose(x, x[0], atol=1e-6)
    # y separates the rows of modules and is the same across each row
    assert (np.diff(y, axis=0) > 0).all()
    assert np.allclose(y.T, y[:, 0], atol=1e-6)


def test_global_pixels_need_a_common_slow_direction():
    """Panels facing opposite ways have no two-dimensional layout."""
    panels = [
        FakePanel((1, 0, 0), (0, -1, 0), (0, 10, 100)),
        FakePanel((1, 0, 0), (0, 1, 0), (0, -10, 100)),
    ]
    assert (
        detector_global_pixels(panels, np.zeros(2, dtype=int), np.zeros((2, 2))) is None
    )
