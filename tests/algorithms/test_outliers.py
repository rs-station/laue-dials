from os.path import dirname

import pytest
from dials.array_family import flex
from dxtbx.model import ExperimentList

from laue_dials.algorithms.outliers import gen_kde


@pytest.fixture
def example_experiment_list_reflection_table(request):
    current_dir = str(dirname(dirname(request.path)))
    elist = ExperimentList.from_file(
        current_dir + "/data/poly_refined.expt", check_format=False
    )
    refls = flex.reflection_table.from_file(current_dir + "/data/poly_refined.refl")
    return elist, refls


def test_gen_kde(example_experiment_list_reflection_table):
    elist, refls = example_experiment_list_reflection_table
    normalized_resolution, lams, kde = gen_kde(elist, refls)

    assert len(normalized_resolution) == len(lams)
    assert kde is not None
