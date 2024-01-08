#!/usr/bin/env python
"""
This script handles polychromatic geometry refinement.
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
from dials.command_line.refine import run_dials_refine
from dials.util import show_mail_handle_errors
from dials.util.options import ArgumentParser
from dxtbx.model import ExperimentList
from dxtbx.model.experiment_list import ExperimentListFactory

from laue_dials.algorithms.laue import (
    gen_beam_models,
    remove_beam_models,
    store_wavelengths,
)
from laue_dials.utils.version import laue_version

# Print laue-dials + DIALS versions
laue_version()

logger = logging.getLogger("laue-dials.command_line.refine")

help_message = """
This script handles polychromatic geometry refinement.

This program takes an indexed DIALS experiment list and reflection table
(with wavelengths) and refines the experiment geometry. The outputs are a pair of
files (poly_refined.expt, poly_refined.refl) that contain an experiment geometry
suitable for prediction and integration.

Examples:

    laue.refine [options] optimized.expt optimized.refl
"""

# Set the phil scope
main_phil = libtbx.phil.parse(
    """
include scope dials.command_line.refine.working_phil

nproc = 1
  .type = int
  .help = Number of parallel processes to run
""",
    process_includes=True,
)

refiner_phil = libtbx.phil.parse(
    """

refinement {
  refinery {
    engine = SparseLevMar
  }

  reflections {
    weighting_strategy {
      override = stills
    }

    outlier {
      nproc = 1

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
      action = fix

      min_nref_per_parameter = 1
    }

    spherical_relp_model = True
  }
}

output {
  experiments = poly_refined.expt

  reflections = poly_refined.refl

  log = laue.poly_refined.log
}
"""
)

working_phil = main_phil.fetch(sources=[refiner_phil])


def refine_image(params, expts, refls):
    """
    Refine image given parameters, experiments, and reflections.

    Args:
        params (libtbx.phil.scope_extract): Refinement parameters.
        expts (dxtbx.model.ExperimentList): Experiment list.
        refls (flex.reflection_table): Reflection table.

    Returns:
        refined_expts (dxtbx.model.experiment_list.ExperimentList): The refined experiment list with updated geometry.
        refined_refls (dials.array_family.flex.reflection_table): The refined reflection table with updated wavelength and centroid data.
    """
    img_num = refls["id"][0]
    original_ids = refls["id"]
    refls["id"] = flex.int([0] * len(refls))
    refls["imageset_id"] = flex.int([0] * len(refls))
    # Generate beam models
    multi_expts, multi_refls = gen_beam_models(expts, refls)

    # Perform refinement
    try:
        refined_expts, refined_refls, _, _ = run_dials_refine(
            multi_expts, multi_refls, params
        )
    except:
        logger.warning(
            f"WARNING: Experiment {img_num} could not be refined. Skipping image."
        )
        return ExperimentList(), reflection_table()  # Return empty

    # Write wavelengths and centroid data
    refined_refls = store_wavelengths(refined_expts, refined_refls)
    refined_refls.map_centroids_to_reciprocal_space(refined_expts)

    # Strip beam objects and reset reflection IDs
    refined_expts = remove_beam_models(refined_expts, original_ids[0])
    refined_refls["id"] = original_ids

    # Return refined data
    return refined_expts, refined_refls


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    """
    Run the refinement script.

    Args:
        args (list): Command-line arguments.
        phil: Working phil scope.
    """
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

    # Load files
    input_expts = ExperimentListFactory.from_json_file(
        params.input.experiments[0].filename, check_format=False
    )
    input_refls = reflection_table.from_file(params.input.reflections[0].filename)
    input_refls = input_refls.select(
        input_refls["wavelength"] != 0
    )  # Remove unindexed reflections

    # Get initial time for process
    start_time = time.time()

    # Prepare parallel input
    ids = list(np.unique(input_refls["id"]).astype(np.int32))
    expts_arr = []
    refls_arr = []
    for i in ids:
        expts_arr.append(ExperimentList([input_expts[i]]))
        refls_arr.append(input_refls.select(flex.bool(input_refls["id"] == i)))
    inputs = list(zip(repeat(params), expts_arr, refls_arr))

    # Refine data
    num_processes = params.nproc
    with Pool(processes=num_processes) as pool:
        output = pool.starmap(refine_image, inputs)

    # Initialize arrays for final results
    total_refined_expts = ExperimentList()
    total_refined_refls = reflection_table()

    # Convert refined data to DIALS objects
    for i in ids:
        total_refined_expts.extend(output[i][0])
        total_refined_refls.extend(output[i][1])

    logger.info("Saving refined experiments to %s", params.output.experiments)
    total_refined_expts.as_file(params.output.experiments)

    logger.info("Saving refined reflections to %s", params.output.reflections)
    total_refined_refls.as_file(filename=params.output.reflections)

    # Final logs
    logger.info("")
    logger.info("Time Taken Refining = %f seconds", time.time() - start_time)


if __name__ == "__main__":
    run()
