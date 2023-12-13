"""
This file contains functions for monochromatic processing
"""


def find_spots(params, expts):
    """
    Find strong reflections on images given a set of experiments.

    Args:
        params (libtbx.phil.scope_extract): A phil scope extract containing the parameters for the DIALS spotfinding code.
        expts (dxtbx.model.ExperimentList): A list of experiment objects.

    Returns:
        refls (dials.array_family.flex.reflection_table): A reflection table containing the found strong reflections.
    """
    from dials.command_line.find_spots import do_spotfinding

    refls = do_spotfinding(expts, params)
    return refls


def initial_index(params, expts, refls):
    """
    Indexes a dataset at a single wavelength using FFT3D.

    Args:
        params (libtbx.phil.scope_extract): A phil scope extract containing the parameters for the DIALS indexing code.
        expts (dxtbx.model.ExperimentList): A list of imported experiment objects.
        refls (dials.array_family.flex.reflection_table): A reflection table containing strong reflections.

    Returns:
        expts_indexed (dxtbx.model.ExperimentList): An ExperimentList containing the indexing solution geometry.
        refls_indexed (dials.array_family.flex.reflection_table): A reflection table containing reflections with indexed data.
    """
    from dials.command_line.index import index

    expts_indexed, refls_indexed = index(expts, refls, params)
    return expts_indexed, refls_indexed


def scan_varying_refine(params, expts, refls):
    """
    Performs scan-varying geometric refinement over a sequence of images.

    Args:
        params (phil): A phil scope containing the needed input for dials.refine.
        expts (dxtbx.model.ExperimentList): A list of experiment objects.
        refls (dials.array_family.flex.reflection_table): A reflection table containing strong reflections.

    Returns:
        expts_refined (dxtbx.model.ExperimentList): An ExperimentList containing the refined geometric solution.
        refls_refined (dials.array_family.flex.reflection_table): A reflection table containing reflections with refined data.
    """
    from dials.command_line.refine import run_dials_refine

    expts_refined, refls_refined, _, _ = run_dials_refine(expts, refls, params)
    return expts_refined, refls_refined
