#!/usr/bin/env python
"""
This script optimizes Miller indices and wavelengths jointly.
"""

import logging
import sys
import time

import libtbx.phil
from dials.array_family.flex import reflection_table
from dials.util import show_mail_handle_errors
from dials.util.options import (ArgumentParser,
                                reflections_and_experiments_from_files)
from dxtbx.model import ExperimentList

logger = logging.getLogger("laue-dials.command_line.optimize_indexing")

help_message = """

This program takes initial monochromatic estimates of the geometry in a
DIALS experiment list and reflection table, and optimizes the indexed
solution by allowing the wavelength to vary. Only Miller indices, s1
vectors and wavelengths for scattered reflections are overwritten by
this program - all other geometry variables remain constant here.

The outputs are a pair of files (optimized.expt, optimized.refl), which
mimic the input files but with updated Miller indices, s1 vectors, and
wavelengths in the reflection table. Note that the experiment list is
unchanged entirely.

Examples::

    laue.optimize_indexing [options] monochromatic.expt monochromatic.refl
"""

# Set the phil scope
phil_scope = libtbx.phil.parse(
    """

output {
  experiments = 'optimized.expt'
    .type = str
    .help = "The output experiment list filename."

  reflections = 'optimized.refl'
    .type = str
    .help = "The output reflection table filename."

  log = 'laue.optimize_indexing.log'
    .type = str
    .help = "The log filename."
  }

n_macrocycles = 3
  .type = int(value_min=1)
  .help = "Number of macrocycles of index optimization to perform"

wavelengths {
  lam_min = 0.95
    .type = float(value_min=0.1)
    .help = "Minimum wavelength for beam spectrum"
  lam_max = 1.15
    .type = float(value_min=0.2)
    .help = "Maximum wavelength for beam spectrum"
  }

reciprocal_grid {
  d_min = 1.4
    .type = float(value_min=0.1)
    .help = "Minimum d-spacing for reflecting planes"
  }

""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


def index_image(params, refls, expts):
    #
    import gemmi
    import numpy as np
    from cctbx.sgtbx import space_group
    from cctbx.uctbx import unit_cell
    from dials.array_family import flex
    from dxtbx.model import ExperimentList
    from tqdm import trange

    from laue_dials.algorithms.laue import LaueAssigner

    if type(expts) != ExperimentList:
        expts = ExperimentList([expts])

    s0 = np.array(expts[0].beam.get_s0())

    # Write to reflection file
    refls["wavelength"] = flex.double(len(refls))
    refls["miller_index"] = flex.miller_index(len(refls))
    refls["harmonics"] = flex.bool([False] * len(refls))

    for i in range(len(expts.imagesets())):
        # Get experiment data from experiment objects
        experiment = expts[i]
        cryst = experiment.crystal
        spacegroup = gemmi.SpaceGroup(
            cryst.get_space_group().type().universal_hermann_mauguin_symbol()
        )

        # Get reflections on this image
        idx = refls["id"] == int(experiment.identifier)
        subrefls = refls.select(idx)

        # Get unit cell params
        cell_params = cryst.get_unit_cell().parameters()
        cell = gemmi.UnitCell(*cell_params)

        # Generate s vectors
        s1 = subrefls["s1"].as_numpy_array()

        # Get U matrix
        U = np.asarray(cryst.get_U()).reshape(3, 3)

        # Generate assigner object
        la = LaueAssigner(
            s0,
            s1,
            cell,
            U,
            params.wavelengths.lam_min,
            params.wavelengths.lam_max,
            params.reciprocal_grid.d_min,
            spacegroup,
        )

        # Optimize Miller indices
        la.assign()
        for j in trange(params.n_macrocycles):
            la.reset_inliers()
            la.update_rotation()
            la.assign()
            la.reject_outliers()
            la.update_rotation()
            la.assign()

        # Update s1 based on new wavelengths
        s1[la._inliers] = la.s1 / la.wav[:, None]

        # Reset crystal parameters based on new geometry
        cryst.set_U(la.R.flatten())
        cryst.set_A(la.RB.flatten())
        cryst.set_B(la.B.flatten())
        cryst.set_space_group(space_group(la.spacegroup.hall))
        cryst.set_unit_cell(unit_cell(la.cell.parameters))

        # Write data to reflections
        refls["s1"].set_selected(idx, flex.vec3_double(s1))
        refls["miller_index"].set_selected(
            idx,
            flex.miller_index(la._H.astype("int").tolist()),
        )
        refls["wavelength"].set_selected(
            idx,
            flex.double(la._wav.tolist()),
        )
        refls["harmonics"].set_selected(
            idx,
            flex.bool(la._harmonics.tolist()),
        )

    # Return reindexed expts, refls
    return expts, refls


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    # Parse arguments
    usage = "laue.optimize_indexing [options] monochromatic.expt monochromatic.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=True)

    # Configure logging
    console = logging.StreamHandler(sys.stdout)
    fh = logging.FileHandler(params.output.log, mode="w", encoding="utf-8")
    loglevel = logging.INFO

    logger.addHandler(fh)
    logger.addHandler(console)
    logging.captureWarnings(True)
    warning_logger = logging.getLogger("py.warnings")
    warning_logger.addHandler(fh)
    warning_logger.addHandler(console)
    dials_logger = logging.getLogger("dials")
    dials_logger.addHandler(fh)
    dials_logger.addHandler(console)
    dxtbx_logger = logging.getLogger("dxtbx")
    dxtbx_logger.addHandler(fh)
    dxtbx_logger.addHandler(console)
    xfel_logger = logging.getLogger("xfel")
    xfel_logger.addHandler(fh)
    xfel_logger.addHandler(console)

    logger.setLevel(loglevel)
    dials_logger.setLevel(loglevel)
    dxtbx_logger.setLevel(loglevel)
    xfel_logger.setLevel(loglevel)
    fh.setLevel(loglevel)

    # Log diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Print help if no input
    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        return

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    reflections = reflections[0]  # Get reflection table out of list

    # Sanity checks
    if len(experiments) == 0:
        parser.print_help()
        return

    # Get initial time for process
    start_time = time.time()

    # This will populate ['s1'] & ['rlp'] columns
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)

    # Loop over input files
    total_experiments = ExperimentList()
    total_reflections = reflection_table()
    for i in range(len(experiments)):
        # Reindex data
        logger.info(f"Reindexing experiment {i}.")
        j = int(experiments[i].identifier)
        indexed_experiments, indexed_reflections = index_image(
            params, reflections.select(reflections["id"] == j), experiments[i]
        )
        total_experiments.extend(indexed_experiments)
        total_reflections.extend(indexed_reflections)

    # Save experiments
    logger.info("Saving optimized experiments to %s", params.output.experiments)
    total_experiments.as_file(params.output.experiments)

    # Save reflections
    logger.info("Saving optimized reflections to %s", params.output.reflections)
    total_reflections.as_file(filename=params.output.reflections)

    # Final logs
    logger.info("")
    logger.info(
        "Time Taken for Total Processing = %f seconds", time.time() - start_time
    )


if __name__ == "__main__":
    run()
