#!/usr/bin/env python
"""
This script performs spotfinding on an imported experiment.
"""
import logging
import sys
import time

import libtbx.phil
from dials.util import show_mail_handle_errors
from dials.util.options import ArgumentParser

from laue_dials.algorithms.monochromatic import find_spots
from laue_dials.utils.version import laue_version

logger = logging.getLogger("laue-dials.command_line.find_spots")

help_message = """
Perform spotfinding on an imported experiment.

This program takes a DIALS imported experiment list (generated with
dials.import) and generates a reflection table of strong spots
for the remainder of the pipeline. The output is a strong.refl file
that contains all the found strong spots for the experiment.

Examples:

    laue.find_spots [options] imported.expt
"""

# Set the phil scope
main_phil = libtbx.phil.parse(
    """
include scope dials.command_line.find_spots.phil_scope
""",
    process_includes=True,
)

output_phil = libtbx.phil.parse(
    """
output {
  log = 'laue.find_spots.log'
    .type = str
    .help = "The log filename."
}
"""
)

spotfinder_phil = libtbx.phil.parse(
    """
spotfinder {
  force_2d = True
}

output {
  shoeboxes = False
}
"""
)

working_phil = main_phil.fetch(sources=[output_phil, spotfinder_phil])


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    """
    Run the spotfinding script with the specified command-line arguments.

    Args:
        args (list): Command-line arguments.
        phil: The phil scope for the program.

    Returns:
        None
    """
    # Parse arguments
    usage = "laue.find_spots [options] imported.expt"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=False,
        read_experiments=True,
        check_format=True,
        epilog=help_message,
    )
    params, options = parser.parse_args(
        args=args, show_diff_phil=True, ignore_unhandled=True
    )

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

    # Print version information
    logger.info(laue_version())

    # Log diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Print help if no input
    if not params.input.experiments:
        parser.print_help()
        return

    # Import images into expt file
    imported_expts = params.input.experiments[0][1]

    # Get initial time for process
    time.time()

    # Find strong spots
    spotfinding_time = time.time()

    logger.info("")
    logger.info("*" * 80)
    logger.info("Finding strong spots")
    logger.info("*" * 80)

    find_spots(params, imported_expts)

    logger.info("")
    logger.info("Time Taken Spotfinding = %f seconds", time.time() - spotfinding_time)


if __name__ == "__main__":
    run()
