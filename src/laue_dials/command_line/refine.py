#!/usr/bin/env python
"""
This script handles polychromatic geometry refinement.
"""
import logging
import time

import libtbx.phil
from dials.array_family.flex import reflection_table
from dials.command_line.refine import run_dials_refine
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser
from dxtbx.model.experiment_list import ExperimentListFactory

from laue_dials.algorithms.laue import (gen_beam_models, remove_beam_models,
                                        store_wavelengths)

logger = logging.getLogger("laue-dials.command_line.refine")

help_message = """

This program takes an indexed DIALS experiment list and reflection table
(with wavelengths) and refines the experiment geometry. The outputs are a pair of
files (poly_refined.expt, poly_refined.refl) that contain an experiment geometry
suitable for prediction and integration.

Examples:

    laue.refine [options] optimized.expt optimized.refl
"""

# Set the phil scope
master_phil = libtbx.phil.parse(
    """
include scope dials.command_line.refine.phil_scope

store_beams = False
  .type = bool
  .help = "Whether to store beam objects in expt file"
""",
    process_includes=True,
)

refiner_phil = libtbx.phil.parse(
    """
refinery {
  engine = SparseLevMar
}

refinement {
  reflections {
    weighting_strategy {
      override = stills
    }

    outlier {
      minimum_number_of_reflections = 1

      algorithm = mcd

      separate_images = True
    }
  }

  parameterisation {
    beam {
      fix = *in_spindle_plane *out_spindle_plane
    }

    crystal {
      unit_cell {
        fix_list = real_space_a
      }
    }

    detector {
      fix = distance
    }

    auto_reduction {
      action = remove

      min_nref_per_parameter = 1
    }
  }
}

output {
  experiments = poly_refined.expt

  reflections = poly_refined.refl
}
"""
)

working_phil = master_phil.fetch(sources=[refiner_phil])


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    # Parse arguments
    usage = "laue.refine [options] optimized.expt optimized.refl"

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

    # Load files
    expts = ExperimentListFactory.from_json_file(
        params.input.experiments[0].filename, check_format=True
    )
    refls = reflection_table.from_file(params.input.reflections[0].filename)

    # Get initial time for process
    start_time = time.time()

    # Generate beam models
    logger.info("")
    logger.info("*" * 80)
    logger.info("Generating beam models")
    logger.info("*" * 80)

    multi_expts, multi_refls = gen_beam_models(expts, refls)

    # Perform scan-varying refinement
    time.time()

    logger.info("")
    logger.info("*" * 80)
    logger.info("Performing geometric refinement")
    logger.info("*" * 80)

    refined_expts, refined_refls, _, _ = run_dials_refine(
        multi_expts, multi_refls, params
    )

    logger.info("")
    logger.info("*" * 80)
    logger.info("Storing refined wavelengths")
    logger.info("*" * 80)

    refined_refls = store_wavelengths(refined_expts, refined_refls)

    # Strip beam objects and reset reflection IDs
    if not params.store_beams:
        logger.info("")
        logger.info("*" * 80)
        logger.info("Removing beam objects")
        logger.info("*" * 80)

        refined_refls["id"] = refls["id"]
        refined_expts = remove_beam_models(refined_expts)

    logger.info("Saving refined experiments to %s", params.output.experiments)
    refined_expts.as_file(params.output.experiments)

    logger.info("Saving refined reflections to %s", params.output.reflections)
    refined_refls.as_file(filename=params.output.reflections)

    # Final logs
    logger.info("")
    logger.info("Time Taken Refining = %f seconds", time.time() - start_time)


if __name__ == "__main__":
    run()
