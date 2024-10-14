#!/usr/bin/env python
"""
This script performs monochromatic indexing and optional scan-varying refinement.
"""
import logging
import sys
import time

import libtbx.phil
from dials.util import show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

from laue_dials.algorithms.monochromatic import initial_index, scan_varying_refine
from laue_dials.utils.version import laue_version

# Print laue-dials + DIALS versions
laue_version()

logger = logging.getLogger("laue-dials.command_line.index")

help_message = """
Perform monochromatic indexing and optional scan-varying refinement.

This program takes a DIALS imported experiment list (generated with
dials.import) and a strong reflection table and generates an initial
monochromatic indexing solution to feed into the remainder of the pipeline.
The outputs are a pair of files (monochromatic.expt, monochromatic.refl)
that constitute a monochromatic estimate of a geometric solution for the
experiment.

Examples:

    laue.index [options] imported.expt strong.refl
"""

# Set the phil scope
main_phil = libtbx.phil.parse(
    """
laue_output {

  index_only = True
    .type = bool
    .help = "Whether to only index or also refine the data."

  final_output {
      experiments = 'monochromatic.expt'
        .type = str
        .help = "The final output experiment list filename."
      reflections = 'monochromatic.refl'
        .type = str
        .help = "The final output reflection table filename."
  }

  indexed {
    experiments = 'indexed.expt'
      .type = str
      .help = "The output indexed experiment list filename."

    reflections = 'indexed.refl'
      .type = str
      .help = "The output indexed reflection table filename."
  }

  refined {
      experiments = 'refined.expt'
        .type = str
        .help = "The output refined experiment list stills filename."
      reflections = 'refined.refl'
        .type = str
        .help = "The output refined reflection table stills filename."
  }

  log = 'laue.index.log'
    .type = str
    .help = "The log filename."
}

indexer {
  include scope dials.command_line.index.phil_scope
}

refiner {
  include scope dials.command_line.refine.phil_scope
}
""",
    process_includes=True,
)

indexer_phil = libtbx.phil.parse(
    """
indexer {
  indexing {
    refinement_procotol {
      mode = repredict_only
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

refiner_phil = libtbx.phil.parse(
    """
refiner {
  refinement {
    parameterisation {
      goniometer {
        fix = None
      }
      beam {
        fix = all
      }
      crystal {
        fix = cell
      }
      detector {
        fix = orientation
      }
      scan_varying = True
    }
    reflections {
      outlier {
        algorithm = tukey
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

working_phil = main_phil.fetch(sources=[indexer_phil, refiner_phil])


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    """
    Run the indexing script with the specified command-line arguments.

    Args:
        args (list): Command-line arguments.
        phil: The phil scope for the program.

    Returns:
        None
    """
    # Parse arguments
    usage = "laue.index [options] imported.expt strong.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=True,
        epilog=help_message,
    )
    params, options = parser.parse_args(
        args=args, show_diff_phil=True, ignore_unhandled=True
    )

    # Configure logging
    console = logging.StreamHandler(sys.stdout)
    fh = logging.FileHandler(params.laue_output.log, mode="w", encoding="utf-8")
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

    # Get data
    strong_refls, imported_expts = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    # Get initial time for process
    start_time = time.time()

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

    # In case indexing is last step
    final_expts = indexed_expts
    final_refls = indexed_refls

    if not params.laue_output.index_only:
        # Perform scan-varying refinement
        refine_time = time.time()

        logger.info("")
        logger.info("*" * 80)
        logger.info("Performing geometric refinement")
        logger.info("*" * 80)

        refined_expts, refined_refls = scan_varying_refine(
            params.refiner, indexed_expts, indexed_refls
        )

        logger.info(
            "Saving refined experiments to %s", params.laue_output.refined.experiments
        )
        refined_expts.as_file(params.laue_output.refined.experiments)

        logger.info(
            "Saving refined reflections to %s", params.laue_output.refined.reflections
        )
        refined_refls.as_file(filename=params.laue_output.refined.reflections)

        # In case refinement is last step
        final_expts = refined_expts
        final_refls = refined_refls

        logger.info("")
        logger.info("Time Taken Refining = %f seconds", time.time() - refine_time)

    # Final logs
    logger.info(
        "Saving final experiments to %s", params.laue_output.final_output.experiments
    )
    final_expts.as_file(params.laue_output.final_output.experiments)

    logger.info(
        "Saving final reflections to %s", params.laue_output.final_output.reflections
    )
    final_refls.as_file(filename=params.laue_output.final_output.reflections)

    logger.info("")
    logger.info(
        "Time Taken for Total Processing = %f seconds", time.time() - start_time
    )


if __name__ == "__main__":
    run()
