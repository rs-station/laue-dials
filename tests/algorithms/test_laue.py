import numpy as np
import pytest
from dials.array_family.flex import reflection_table
from dxtbx.model import ExperimentList

from laue_dials.algorithms.laue import (LaueAssigner, LaueBase, LauePredictor,
                                        gen_beam_models, remove_beam_models,
                                        store_wavelengths)


@pytest.fixture
def example_experiments_reflections():
    expts = ExperimentList.from_file("../data/optimized.expt", check_format=False)
    refls = reflection_table.from_file("../data/optimized.refl")
    return expts, refls


def test_gen_beam_models(example_experiments_reflections):
    expts, refls = example_experiments_reflections
    new_expts, new_refls = gen_beam_models(expts, refls)
    assert len(new_expts) == len(new_refls)


def test_remove_beam_models(example_experiments_reflections):
    expts, refls = example_experiments_reflections
    original_ids = refls["id"]
    new_expts, new_refls = gen_beam_models(expts, refls)
    final_expts = remove_beam_models(new_expts, original_ids[0])
    assert len(final_expts) == len(expts.imagesets())


def test_store_wavelengths(example_experiments_reflections):
    expts, refls = example_experiments_reflections
    new_expts, new_refls = gen_beam_models(expts, refls)
    final_refls = store_wavelengths(expts, refls)
    assert "wavelength" in final_refls.keys()


def test_LaueBase():
    s0 = np.array([0, 0, 1])
    cell = (1, 1, 1, 90, 90, 90)
    R = np.eye(3)
    lam_min, lam_max, dmin = 0.8, 1.2, 1.0
    spacegroup = "P1"
    laue_base = LaueBase(s0, cell, R, lam_min, lam_max, dmin, spacegroup)
    assert laue_base is not None


def test_LaueAssigner():
    s0 = np.array([0, 0, 1])
    s1 = np.array([[1, 1, 1], [2, 2, 2]])
    cell = (1, 1, 1, 90, 90, 90)
    R = np.eye(3)
    lam_min, lam_max, dmin = 0.8, 1.2, 1.0
    spacegroup = "P1"
    laue_assigner = LaueAssigner(s0, s1, cell, R, lam_min, lam_max, dmin, spacegroup)
    assert laue_assigner is not None


def test_LauePredictor():
    s0 = np.array([0, 0, 1])
    cell = (1, 1, 1, 90, 90, 90)
    R = np.eye(3)
    lam_min, lam_max, dmin = 0.8, 1.2, 1.0
    spacegroup = "P1"
    laue_predictor = LauePredictor(s0, cell, R, lam_min, lam_max, dmin, spacegroup)
    assert laue_predictor is not None
