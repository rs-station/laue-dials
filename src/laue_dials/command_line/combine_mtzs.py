#!/usr/bin/env python
"""
This script combines integrated MTZ files by image into one single MTZ file.
"""

import logging

import libtbx.phil
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser

from tqdm import trange
import reciprocalspaceship as rs

logger = logging.getLogger("laue-dials.command_line.combine_mtzs")

help_message = """

This program takes a set of MTZ files and combines them into a single output
MTZ file.

The output is an MTZ file consisting of all the data in the input MTZ files
suitable for merging in software such as careless.

Examples::

    laue.combine_mtzs [options] *.mtz
"""

# Set the phil scope
phil_scope = libtbx.phil.parse(
    """

output {
  filename = 'total_integrated.mtz'
    .type = str
    .help = "The output MTZ filename."

  log = 'laue.combine_mtzs.log'
    .type = str
    .help = "The log filename."
  }
""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    # Parse arguments
    usage = "laue.combine_mtzs [options] *.mtz"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=False,
        read_experiments=False,
        check_format=False,
        epilog=help_message,
    )

    params, options, files = parser.parse_args(
        args=args, show_diff_phil=True, ignore_unhandled=True, return_unhandled=True
    )

    # Configure logging
    log.config(verbosity=options.verbose, logfile=params.output.log)

    # Log diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Loop over input files
    logger.info("Beginning combination.")
    total_integrated_mtz = rs.read_mtz(files[0])  # First still
    for i in trange(1, len(files)):
        try:
            # Get MTZ data
            img_mtz = rs.read_mtz(files[i])
            img_mtz["BATCH"] = i + 1

            # Add MTZ data to total
            total_integrated_mtz = rs.concat([total_integrated_mtz, img_mtz])
        except:
            # Print integration failure to user
            j = i + 1
            logger.info("Image {j:06d} could not be integrated.")
            continue

    # Fix data type of BATCH
    total_integrated_mtz["BATCH"] = total_integrated_mtz["BATCH"].infer_mtz_dtype()

    # Save reflections
    logger.info("Saving combined MTZ data to %s", params.output.filename)
    total_integrated_mtz.write_mtz(params.output.filename)


if __name__ == "__main__":
    run()
