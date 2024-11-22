#!/usr/bin/env python
"""
This script optimizes Miller indices and wavelengths jointly.
"""

import logging
import sys
import time
import math

import libtbx.phil
import numpy as np
from dials.util import show_mail_handle_errors
from dials.util.options import (ArgumentParser,
                                reflections_and_experiments_from_files)
from dxtbx.model import ExperimentList
from cctbx import crystal_orientation
from scitbx.matrix import sqr

from laue_dials.utils.version import laue_version


# Print laue-dials + DIALS versions
laue_version()

logger = logging.getLogger("laue-dials.command_line.align_crystals")

help_message = """
This script aligns the indexing solutions of all crystals in a data set.

This program handles finding symmetry operation relations between all
crystals in a data set and aligning all indexing solutions to be 
consistent with each other. Useful in the case of space groups with
multiple consistent indexing solutions.

The outputs are a pair of files (aligned.expt, aligned.refl), which
mimic the input files but with updated crystal orientations.

Examples:

    laue.align_crystals [options] monochromatic.expt monochromatic.refl
"""

# Set the phil scope
phil_scope = libtbx.phil.parse(
    """

output {
  experiments = 'aligned.expt'
    .type = str
    .help = "The output experiment list filename."

  reflections = 'aligned.refl'
    .type = str
    .help = "The output reflection table filename."

  log = 'laue.align_crystals.log'
    .type = str
    .help = "The log filename."
  }

image = 1
  .type = int
  .help = Which image in the data set to use as a reference (1-indexed).

""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    """
    Run the script to align crystal orientations.

    Args:
        args (list): Command-line arguments.
        phil: The phil scope for the program.

    Returns:
        None
    """
    # Parse arguments
    usage = "laue.align_crystals [options] monochromatic.expt monochromatic.refl"

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

    # Get initial time for process
    start_time = time.time()

    # Get reference crystal
    ref_cryst = experiments[params.image - 1].crystal
    ref_A = ref_cryst.get_A()
    ref_orientation = crystal_orientation.crystal_orientation(ref_A, crystal_orientation.basis_type.reciprocal)

    # Align all other crystals
    for cryst in experiments.crystals():
        # Get symmetry operation rotation to new orientation
        cryst_orientation = crystal_orientation.crystal_orientation(cryst.get_A(), crystal_orientation.basis_type.reciprocal)
#        sym_op_rot = np.reshape(ref_orientation.best_similarity_transformation(old_orientation, fractional_length_tolerance=np.inf, unimodular_generator_range=1),(3,3))
        sym_op_rot = sqr(cryst_orientation.best_similarity_transformation(ref_orientation, math.inf, 1))
        aligned_orientation = cryst_orientation.change_basis(sym_op_rot)
        cryst.set_U(aligned_orientation.get_U_as_sqr())

#        # Rotate crystal matrices
#        cryst.set_A((np.reshape(cryst.get_A(),(3,3))@sym_op_rot).flatten())

    # Convert reindexed data to DIALS objects
    #total_experiments = ExperimentList()
    #total_reflections = reflection_table()
    total_experiments = experiments
    #total_reflections = reflections
#    for i in ids:
#        if len(output[i][0]) == 0: # Empty experiment list
#            continue
#        total_experiments.extend(output[i][0])
#        total_reflections.extend(output[i][1])

    # Save experiments
    logger.info("Saving aligned experiments to %s", params.output.experiments)
    total_experiments.as_file(params.output.experiments)

#    # Save reflections
#    logger.info("Saving aligned reflections to %s", params.output.reflections)
#    total_reflections.as_file(filename=params.output.reflections)

    # Final logs
    logger.info("")
    logger.info(
        "Time Taken for Total Processing = %f seconds", time.time() - start_time
    )


if __name__ == "__main__":
    run()
