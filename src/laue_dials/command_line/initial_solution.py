#!/usr/bin/env python
"""
This script performs a monochromatic pipeline for an initial solution to feed
into the Laue pipeline.
"""
import logging
import time

import libtbx.phil
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser

from laue_dials.algorithms.monochromatic import (find_spots, initial_index)

logger = logging.getLogger("laue-dials.command_line.initial_solution")

help_message = """

This program takes a DIALS imported experiment list (generated with
dials.import) and generates a DIALS experiment list and reflection
table to use as an initial monochromatic solution to feed
into the remainder of the pipeline. The outputs are a pair of
files (monochromatic.expt, monochromatic.refl) that constitute a
monochromatic estimate of a geometric solution for the crystal and
experiment. Spotfinding and indexing are done in this script.

Examples:

    laue.initial_solution [options] imported.expt
"""

# Set the phil scope
master_phil = libtbx.phil.parse(
    """
laue_output {
  strong_filename = 'strong.refl'
    .type = str
    .help = "The output spotfinding reflection table filename."

  indexed {
    experiments = 'monochromatic.expt'
      .type = str
      .help = "The output indexed experiment list filename."

    reflections = 'monochromatic.refl'
      .type = str
      .help = "The output indexed reflection table filename."
  }

  log = 'laue.initial_solution.log'
    .type = str
    .help = "The log filename."
}

spotfinder {
  include scope dials.command_line.find_spots.phil_scope
}

indexer {
  include scope dials.command_line.index.phil_scope
}
""",
    process_includes=True,
)

spotfinder_phil = libtbx.phil.parse(
    """
spotfinder {
  spotfinder {
    threshold {
      dispersion {
        gain = 0.10
      }
    }

    filter {
      max_separation = 10
    }

    force_2d = True
  }

  output {
    shoeboxes = False
  }
}
"""
)

indexer_phil = libtbx.phil.parse(
    """
indexer {
  indexing {
    refinement_procotol {
      n_macrocycles = 10
    }
  }

  refinement {
    parameterisation {
      beam {
        fix = *all in_spindle_plane out_spindle_plane wavelength
      }

      detector {
        fix = *all position orientation distance
      }

      goniometer {
        fix = *all in_beam_plane out_beam_plane
      }

      scan_varying = True
    }

    reflections {
      outlier {
        algorithm = null auto mcd *tukey sauter_poon

        tukey {
          iqr_multiplier = 0.
        }

        minimum_number_of_reflections = 1
      }
    }
  }
}
"""
)

working_phil = master_phil.fetch(sources=[spotfinder_phil, indexer_phil])


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    # Parse arguments
    usage = "laue.initial_solution [options] image_*.mccd"

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
    log.config(verbosity=options.verbose, logfile=params.laue_output.log)

    # Log diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Import images into expt file
    imported_expts = params.input.experiments[0][1]

    # Get initial time for process
    start_time = time.time()

    # Find strong spots
    spotfinding_time = time.time()

    logger.info("")
    logger.info("*" * 80)
    logger.info("Finding strong spots")
    logger.info("*" * 80)

    strong_refls = find_spots(params.spotfinder, imported_expts)

    logger.info("Saving strong spots to %s", params.laue_output.strong_filename)
    strong_refls.as_file(filename=params.laue_output.strong_filename)

    logger.info("")
    logger.info("Time Taken Spotfinding = %f seconds", time.time() - spotfinding_time)

    # Index at peak wavelength
    index_time = time.time()

    logger.info("")
    logger.info("*" * 80)
    logger.info("Indexing images")
    logger.info("*" * 80)

    indexed_expts, indexed_refls = initial_index(
        params.indexer, imported_expts, strong_refls
    )

    logger.info(
        "Saving indexed experiments to %s", params.laue_output.indexed.experiments
    )
    indexed_expts.as_file(params.laue_output.indexed.experiments)

    logger.info(
        "Saving indexed reflections to %s", params.laue_output.indexed.reflections
    )
    indexed_refls.as_file(filename=params.laue_output.indexed.reflections)

    logger.info("")
    logger.info("Time Taken Indexing = %f seconds", time.time() - index_time)

    # Final logs
    logger.info("")
    logger.info(
        "Time Taken for Total Processing = %f seconds", time.time() - start_time
    )


if __name__ == "__main__":
    run()
