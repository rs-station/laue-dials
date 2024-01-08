#!/usr/bin/env python
"""
This script generates integrated MTZ files from refined data with predictions
"""

import logging
import sys
import time
from functools import partial
from itertools import repeat
from multiprocessing import Pool

import gemmi
import libtbx.phil
import numpy as np
import reciprocalspaceship as rs
from cctbx import sgtbx
from dials.array_family import flex
from dials.util import show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

from laue_dials.algorithms.integration import SegmentedImage
from laue_dials.utils.version import laue_version

# Print laue-dials + DIALS versions
laue_version()

logger = logging.getLogger("laue-dials.command_line.integrate")

help_message = """
This script generates integrated MTZ files from refined data with predictions.

The program takes a refined geometry experiment file along with a predicted
reflection table, and uses those to integrate intensities in the data set.
The output is an MTZ file containing integrated intensities suitable for
merging and scaling.

Examples:

    laue.integrate [options] poly_refined.expt predicted.refl
"""

# Set the phil scope
phil_scope = libtbx.phil.parse(
    """
output {
  filename = 'integrated.mtz'
    .type = str
    .help = "The output MTZ filename."

  log = 'laue.integrate.log'
    .type = str
    .help = "The log filename."
  }

nproc = 1
  .type = int
  .help = "Number of parallel integrations to do"

isigi_cutoff = 2.0
  .type = float
  .help = "I/SIGI threshold to use for marking strong spots."
""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


def get_refls_image(refls, img_id):
    """
    Get the set of reflections lying on a particular image.

    Args:
        refls (dials.array_family.flex.reflection_table): Reflection table.
        img_id (int): Image ID.

    Returns:
        refls (dials.array_family.flex.reflection_table): Reflection table for the specified image.
    """
    return refls.select(refls["id"] == img_id)


def integrate_image(img_set, refls, isigi_cutoff):
    """
    Integrate predicted spots on an image.

    Args:
        img_set (dxtbx_imageset_ext.Imageset): Image set.
        refls (dials.array_family.flex.reflection_table): Reflection table.
        isigi_cutoff (float): I/SIGI threshold.

    Returns:
        flex.reflection_table: Updated reflection table.
    """
    img_num = refls["id"][0]
    logger.info(f"Integrating image {img_num}.")
    proctime = time.time()

    # Make SegmentedImage
    all_spots = refls["xyzcal.px"].as_numpy_array()[:, :2].astype("float32")
    pixels = img_set.get_raw_data(0)[0].as_numpy_array().astype("float32")
    sim = SegmentedImage(pixels, all_spots)

    # Get integrated reflections only
    refls = refls.select(flex.bool(sim.used_reflections))

    # Integrate reflections
    sim.integrate(isigi_cutoff)

    # Update reflection data
    i = np.zeros(len(refls))
    sigi = np.zeros(len(refls))
    bg = np.zeros(len(refls))
    sigbg = np.zeros(len(refls))
    profiles = sim.profiles.to_list()
    for j in range(len(refls)):
        prof = profiles[j]
        if prof.success:
            i[j] = prof.I
            sigi[j] = prof.SigI
            bg[j] = np.maximum((prof.background * prof.bg_mask), 0.0).sum()
            sigbg[j] = np.sqrt(np.maximum((prof.background * prof.bg_mask), 0.0)).sum()
    refls["intensity.sum.value"] = flex.double(i)
    refls["intensity.sum.variance"] = flex.double(sigi**2)
    refls["background.sum.value"] = flex.double(bg)
    refls["background.sum.variance"] = flex.double(sigbg**2)
    refls = refls.select(refls["intensity.sum.value"] != 0)
    logger.info(f"Image {img_num} took {time.time() - proctime} seconds.")
    return refls  # Updated reflection table


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    """
    Run the integration script with the specified command-line arguments.

    Args:
        args (list): Command-line arguments.
        phil: The phil scope for the program.

    Returns:
        None
    """
    # Parse arguments
    usage = "laue.integrate [options] poly_refined.expt predicted.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=True,
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

    # Load data
    reflections, expts = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    preds = reflections[0]  # Get predictions

    # Sanity checks
    if len(expts) == 0:
        parser.print_help()
        return

    # Get reflections and image data
    imagesets = expts.imagesets()
    ids = list(np.unique(preds["id"]).astype(np.int32))
    get_refls = partial(get_refls_image, preds)
    tables = list(map(get_refls, ids))
    inputs = list(zip(imagesets, tables, repeat(params.isigi_cutoff)))

    # Get initial time for process
    start_time = time.time()

    # Multiprocess integration
    num_processes = params.nproc
    logger.info("Starting integration.")
    with Pool(processes=num_processes) as pool:
        refls_arr = pool.starmap(integrate_image, inputs, chunksize=1)
    logger.info("Integration finished.")

    # Construct an integrated reflection table
    logger.info("Constructing reflection table")
    final_refls = flex.reflection_table()
    for refls in refls_arr:
        final_refls.extend(refls)
    refls = final_refls

    # Get data needed for MTZ file
    logger.info("Converting to MTZ format.")
    hkl = refls["miller_index"].as_vec3_double()
    cell = np.zeros(6)
    for crystal in expts.crystals():
        cell += np.array(crystal.get_unit_cell().parameters()) / len(expts.crystals())
    cell = gemmi.UnitCell(*cell)
    sginfo = expts.crystals()[0].get_space_group().info()
    symbol = sgtbx.space_group_symbols(sginfo.symbol_and_number().split("(")[0])
    spacegroup = gemmi.SpaceGroup(symbol.universal_hermann_mauguin())

    # Generate rs.DataSet to write to MTZ
    data = rs.DataSet(
        {
            "H": hkl.as_numpy_array()[:, 0].astype(np.int32),
            "K": hkl.as_numpy_array()[:, 1].astype(np.int32),
            "L": hkl.as_numpy_array()[:, 2].astype(np.int32),
            "BATCH": refls["id"].as_numpy_array() + 1,
            "I": refls["intensity.sum.value"].as_numpy_array(),
            "SIGI": refls["intensity.sum.variance"].as_numpy_array() ** 0.5,
            "xcal": refls["xyzcal.px"].as_numpy_array()[:, 0],
            "ycal": refls["xyzcal.px"].as_numpy_array()[:, 1],
            "wavelength": refls["wavelength"].as_numpy_array(),
            "BG": refls["background.sum.value"].as_numpy_array(),
            "SIGBG": refls["background.sum.variance"].as_numpy_array() ** 0.5,
        },
        cell=cell,
        spacegroup=spacegroup,
    ).infer_mtz_dtypes()

    # Save reflections
    logger.info("Saving integrated reflections to %s", params.output.filename)
    data.write_mtz(params.output.filename, skip_problem_mtztypes=True)

    # Final logs
    logger.info("")
    logger.info(
        "Time Taken for Total Processing = %f seconds", time.time() - start_time
    )


if __name__ == "__main__":
    run()
