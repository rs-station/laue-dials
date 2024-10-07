#!/usr/bin/env python
"""
This script histograms the wavelengths present in a reflection table.
"""

import logging
import sys

import libtbx.phil
from dials.util import show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from matplotlib import pyplot as plt

from laue_dials.utils.version import laue_version

# Print laue-dials + DIALS versions
laue_version()

logger = logging.getLogger("laue-dials.command_line.plot_wavelengths")

help_message = """
This script generates a histogram of the wavelengths present in a DIALS reflection table.

This program takes a DIALS reflection table that contains a wavelength
column and produces a histogram of the wavelengths present in the
reflection table.

Examples:

    laue.plot_wavelengths [options] FILENAME.refl
"""

# Set the phil scope
phil_scope = libtbx.phil.parse(
    """

  show = True
    .type = bool
    .help = "Whether to display the generated histogram."

  save = False
    .type = bool
    .help = "Whether to save the generated histogram."

  output = "wavelengths.png"
    .type = str
    .help = "The filename for the generated histogram."

  refined_only = False
    .type = bool
    .help = "Only histogram wavelengths of refined reflections."

  n_bins = 30
    .type = int
    .help = "The number of bins to histogram wavelengths with."

  lam_min = None
    .type = int
    .help = "The lower limit for wavelength plot in Angstroms."

  lam_max = None
    .type = int
    .help = "The upper limit for wavelength plot in Angstroms."

  log = 'laue.plot_wavelengths.log'
    .type = str
    .help = "The log filename."

""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    """
    Run the script to histogram the wavelengths present in a reflection table.

    Args:
        args (list): Command-line arguments.
        phil: The phil scope for the program.

    Returns:
        None
    """
    # Parse arguments
    usage = "laue.plot_wavelengths [options] FILENAME.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=working_phil,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure logging
    console = logging.StreamHandler(sys.stdout)
    fh = logging.FileHandler(params.log, mode="w", encoding="utf-8")
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

    # Load reflection data
    if not params.input.reflections:
        parser.print_help()
        return

    refls, _ = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    refls = refls[0]

    #    refls = flex.reflection_table.from_file(params.input.reflections)
    if params.refined_only:
        refls = refls.select(refls.get_flags(refls.flags.used_in_refinement))

    if len(refls) == 0:
        logger.info("No reflections in table after filtering.")
        return

    # Read wavelengths
    try:
        lams = refls["wavelength"].as_numpy_array()
    except:
        logger.info(
            "Wavelength data could not be read. Please check to ensure there is a populated wavelength column in the reflection table."
        )
        return

    # Apply wavelength restrictions
    initial_count = len(lams)

    if params.lam_min is not None:
        lams = lams[lams >= params.lam_min]

    if params.lam_max is not None:
        lams = lams[lams <= params.lam_max]

    if initial_count != len(lams):
        logger.info(
            "Wavelength restrictions have removed some reflections. Plot will not reflect entirety of the data."
        )

    fig = plt.figure()
    plt.hist(lams, bins=params.n_bins)
    plt.title("Wavelength Spectrum")
    plt.xlabel("Wavelength (Angstroms)")
    plt.ylabel("Number of reflections")
    if params.save:
        fig.savefig(params.output)
    if params.show:
        plt.show()


if __name__ == "__main__":
    run()
