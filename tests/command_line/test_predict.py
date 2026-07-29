import numpy as np

from laue_dials.command_line.predict import (
    estimate_integration_radius,
    unmasked_prediction_selection,
)


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


def test_estimate_integration_radius_matches_formula():
    """
    estimate_integration_radius should match the same nearest-neighbor
    spacing formula used by IntegratorBase's default radius.
    """
    from scipy.spatial.distance import pdist, squareform

    rng = np.random.default_rng(0)
    x = rng.uniform(0, 1000, size=200)
    y = rng.uniform(0, 1000, size=200)

    centroids = np.column_stack([x, y])
    dmat = squareform(pdist(centroids))
    closest_spot_dist = np.sort(dmat, axis=0)[1]
    expected = int(np.round(0.5 * np.percentile(closest_spot_dist, 20)))

    assert estimate_integration_radius(x, y) == expected
