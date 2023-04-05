#!/usr/bin/env python
"""
SHORT_DESCRIPTION

Examples
--------

Usage Details
-------------
"""

import logging

import libtbx.phil
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser

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

    laue.optimize_indexing monochromatic.expt monochromatic.refl
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
""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


def index_image(expts, refls):
    print("To be implemented.")


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil, return_results=False):
    # Parse arguments
    usage = (
        "usage: laue.optimize_indexing [options] monochromatic.expt monochromatic.refl"
    )

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
    log.config(verbosity=options.verbose, logfile=params.output.log)

    # Log diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    # Sanity checks
    if len(experiments) == 0:
        parser.print_help()
        return

    # TODO: Optimize indexing routine here

    # Save experiments
    logger.info("Saving optimized experiments to %s", params.output.experiments)
    indexed_experiments.as_file(params.output.experiments)

    # Save reflections
    logger.info("Saving optimized reflections to %s", params.output.reflections)
    indexed_reflections.as_file(filename=params.output.reflections)


if __name__ == "__main__":
    run()
