import pytest

from laue_dials.algorithms.monochromatic import (find_spots, initial_index,
                                                 scan_varying_refine)


@pytest.fixture
def example_params_experiments_reflections():
    params = None
    expts = None
    refls = None
    return params, expts, refls


def test_find_spots():
    params, expts, _ = example_params_experiments_reflections
    refls = find_spots(params, expts)
    assert len(refls) > 0


def test_initial_index():
    params, expts, refls = example_params_experiments_reflections
    expts_indexed, refls_indexed = initial_index(params, expts, refls)
    assert len(expts_indexed) > 0
    assert len(refls_indexed) > 0


def test_scan_varying_refine():
    params, expts, refls = example_params_experiments_reflections
    expts_refined, refls_refined = scan_varying_refine(params, expts, refls)
    assert len(expts_refined) > 0
    assert len(refls_refined) > 0
