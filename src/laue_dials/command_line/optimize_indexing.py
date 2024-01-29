#!/usr/bin/env python
"""
This script optimizes Miller indices and wavelengths jointly.
"""

import logging
import sys
import time
from itertools import repeat
from multiprocessing import Pool

import libtbx.phil
import numpy as np
from dials.array_family import flex
from dials.array_family.flex import reflection_table
from dials.util import show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dxtbx.model import ExperimentList

from laue_dials.utils.version import laue_version

# Print laue-dials + DIALS versions
laue_version()

logger = logging.getLogger("laue-dials.command_line.optimize_indexing")

help_message = """
This script optimizes Miller indices and wavelengths jointly.

This program takes initial monochromatic estimates of the geometry in a
DIALS experiment list and reflection table, and optimizes the indexed
solution by allowing the wavelength to vary. Only Miller indices, s1
vectors, and wavelengths for scattered reflections are overwritten by
this program - all other geometry variables remain constant here.

The outputs are a pair of files (optimized.expt, optimized.refl), which
mimic the input files but with updated Miller indices, s1 vectors, and
wavelengths in the reflection table. Note that the experiment list is
unchanged entirely.

Examples:

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

geometry {
  unit_cell = None
    .type = floats(size=6)
    .help = "Target unit cell for indexing."
}

filter_spectrum = True
  .type = bool
  .help = Whether to remove reflections outside of the provided wavelength limits.

keep_unindexed = False
  .type = bool
  .help = Whether to keep unindexed reflections.

nproc = 1
  .type = int
  .help = Number of parallel processes to run.

n_macrocycles = 3
  .type = int(value_min=1)
  .help = "Number of macrocycles of index optimization to perform."

wavelengths {
  lam_min = None
    .type = float(value_min=0.1)
    .help = "Minimum wavelength for beam spectrum."
  lam_max = None
    .type = float(value_min=0.2)
    .help = "Maximum wavelength for beam spectrum."
  }

reciprocal_grid {
  d_min = None
    .type = float(value_min=0.1)
    .help = "Minimum d-spacing for reflecting planes."
  }

""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


def index_image(params, refls, expts):
    """
    Reindex the given image reflections and update the experiment geometry.

    Args:
        params (libtbx.phil.scope_extract): Program parameters.
        refls (dials.array_family.flex.reflection_table): Reflection table for a single image.
        expts (dxtbx.model.ExperimentList): List of experiment objects.

    Returns:
        expts (dxtbx.model.experiment_list.ExperimentList): The optimized experiment list with updated crystal rotations.
        refls (dials.array_family.flex.reflection_table): The optimized reflection table with updated wavelengths.
    """
    import gemmi
    from cctbx.sgtbx import space_group
    from cctbx.uctbx import unit_cell
    from dials.array_family import flex
    from dxtbx.model import ExperimentList

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
        if params.geometry.unit_cell is not None:
            cell_params = params.geometry.unit_cell
            cell = gemmi.UnitCell(*cell_params)
            # Check compatibility with spacegroup
            if not cell.is_compatible_with_spacegroup(spacegroup):
                logger.warning(
                    f"WARNING: User-provided unit cell is incompatible with crystal space group on image {i}. Using crystal unit cell instead."
                )
                cell_params = cryst.get_unit_cell().parameters()
                cell = gemmi.UnitCell(*cell_params)
        else:
            cell_params = cryst.get_unit_cell().parameters()
            cell = gemmi.UnitCell(*cell_params)

        # Generate s vectors
        s1 = subrefls["s1"].as_numpy_array()

        # Get U matrix
        U = np.asarray(cryst.get_U()).reshape(3, 3)

        # Generate assigner object
        logger.info(f"Reindexing image {experiment.identifier}.")
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
        for j in range(params.n_macrocycles):
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

        # Get wavelengths
        spot_wavelengths = np.asarray(la._wav.tolist())

        # Write data to reflections
        refls["s1"].set_selected(idx, flex.vec3_double(s1))
        refls["miller_index"].set_selected(
            idx,
            flex.miller_index(la._H.astype("int").tolist()),
        )
        refls["harmonics"].set_selected(
            idx,
            flex.bool(la._harmonics.tolist()),
        )
        refls["wavelength"].set_selected(
            idx,
            flex.double(spot_wavelengths),
        )

    # Remove unindexed reflections
    if params.filter_spectrum:
        all_wavelengths = refls["wavelength"].as_numpy_array()
        keep = np.logical_and(
            all_wavelengths >= params.wavelengths.lam_min,
            all_wavelengths <= params.wavelengths.lam_max,
        )
        refls = refls.select(flex.bool(keep))

    # Remove unindexed reflections
    if not params.keep_unindexed:
        all_wavelengths = refls["wavelength"].as_numpy_array()
        keep = all_wavelengths > 0  # Unindexed reflections assigned wavelength of 0
        refls = refls.select(flex.bool(keep))

    # Return reindexed expts, refls
    return expts, refls


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    """
    Run the script to optimize Miller indices and wavelengths jointly.

    Args:
        args (list): Command-line arguments.
        phil: The phil scope for the program.

    Returns:
        None
    """
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

    params, options = parser.parse_args(args=args, show_diff_phil=False)

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

    # Check for valid parameter values
    if params.reciprocal_grid.d_min == None:
        logger.info("Please provide a d_min.")
        return
    elif params.wavelengths.lam_min == None or params.wavelengths.lam_max == None:
        logger.info(
            "Please provide upper and lower boundaries for the wavelength spectrum."
        )
        return
    elif params.wavelengths.lam_min > params.wavelengths.lam_max:
        logger.info("Minimum wavelength cannot be greater than maximum wavelength.")
        return

    # Get initial time for process
    start_time = time.time()

    # This will populate ['s1'] & ['rlp'] columns
    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)

    # Prepare parallel input
    ids = list(np.unique(reflections["id"]).astype(np.int32))
    expts_arr = []
    refls_arr = []
    for i in ids:  # Split DIALS objects into lists
        expts_arr.append(ExperimentList([experiments[i]]))
        refls_arr.append(reflections.select(reflections["id"] == i))
    inputs = list(zip(repeat(params), refls_arr, expts_arr))

    # Reindex data
    num_processes = params.nproc
    logger.info("Reindexing images.")
    with Pool(processes=num_processes) as pool:
        output = pool.starmap(index_image, inputs, chunksize=1)
    logger.info(f"All images reindexed.")

    # Convert reindexed data to DIALS objects
    total_experiments = ExperimentList()
    total_reflections = reflection_table()
    for i in ids:
        total_experiments.extend(output[i][0])
        total_reflections.extend(output[i][1])

    # Give all unindexed experiments wavelength 0
    hkls = np.asarray(total_reflections["miller_index"])
    sel = np.all(hkls == [0, 0, 0], axis=1)
    lams = np.asarray(total_reflections["wavelength"])
    lams[sel] = 0.0
    total_reflections["wavelength"] = flex.double(lams)

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
