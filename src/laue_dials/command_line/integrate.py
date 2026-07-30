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

from laue_dials.algorithms.integration import (
    Integrator,
    estimate_integration_radius,
    unmasked_prediction_selection,
)
from laue_dials.utils.version import laue_version

logger = logging.getLogger("laue-dials.command_line.integrate")

help_message = """
This script generates integrated MTZ files from refined data with predictions.

The program takes a refined geometry experiment file along with a predicted
reflection table, and uses those to integrate intensities in the data set.
The output is an MTZ file containing integrated intensities suitable for
merging and scaling.

The algorithm applied here is a variable elliptical profile fitting
algorithm inspired by the VariableElliptical mode in Precognition. Each
predicted centroid is given a circular window of pixels, and the counts in
that window are modeled as an elliptical two-dimensional Gaussian profile
on a flat background, assuming Poisson noise.

Profile shapes are estimated jointly with the background and the
intensities, over at most maxiter iterations. The shape for every
reflection is pooled from the pixels of its knn nearest strong spots, so
weak reflections inherit a well-determined profile from their neighbors
rather than fitting noise. Intensities and their uncertainties are then
obtained by profile fitting, weighting each pixel by its expected
contribution, rather than by summing counts inside a mask.

Unless integration_radius is set, the window radius is estimated from the
spacing of the predicted centroids. That same radius is used to dilate the
detector mask, so predictions whose window would overlap a bad pixel are
discarded before integration.

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

  reflections = None
    .type = str
    .help = "The output reflection table filename. None will output no reflection table."

  log = 'laue.integrate.log'
    .type = str
    .help = "The log filename."
  }

nproc = 1
  .type = int
  .help = "Number of parallel integrations to do"

isigi_cutoff = 3.0
  .type = float
  .help = "I/SIGI threshold to use for marking strong spots."

integration_radius = None
  .type = int(value_min=0)
  .help = "Radius in pixels used both for the integration window around each predicted centroid and for dilating the detector mask when discarding predictions that fall in masked regions. Defaults to a dynamically-computed radius (0.5 * the 20th percentile of nearest-neighbor centroid distances)."

knn = 5
  .type = int(value_min=1)
  .help = "Number of nearest strong spots whose pixels are pooled to estimate the elliptical profile of each reflection. Larger values give steadier profiles but blur genuine variation in spot shape across the detector. Reduced automatically if an image has fewer strong spots than this."

maxiter = 2
  .type = int(value_min=1)
  .help = "Maximum number of profile-fitting iterations. Each iteration re-estimates the background, profiles, and intensities, then re-marks strong spots. Fitting stops early if the Poisson log-likelihood stops improving."
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


def integrate_image(img_set, refls, isigi_cutoff, integration_radius, knn, maxiter):
    """
    Integrate predicted spots on an image.

    The integration radius is estimated once from the full predicted set (or
    taken from ``integration_radius`` if supplied) and reused both to dilate the
    detector mask -- discarding predictions whose integration window would
    overlap masked pixels -- and as the integration window itself. Using a
    single radius for both keeps the dilation and the window self-consistent, so
    no surviving centroid's window ever reaches a masked pixel.

    Args:
        img_set (dxtbx_imageset_ext.Imageset): Image set.
        refls (dials.array_family.flex.reflection_table): Reflection table.
        isigi_cutoff (float): I/SIGI threshold.
        integration_radius (int): Radius in pixels for the integration window
            and mask dilation. If None, a radius is estimated from the spacing
            of the centroids.
        knn (int): Number of nearest strong spots pooled to estimate each
            reflection's profile.
        maxiter (int): Maximum number of profile-fitting iterations.

    Returns:
        flex.reflection_table: Updated reflection table.
    """
    img_num = refls["id"][0]
    logger.info(f"Integrating image {img_num}.")
    proctime = time.time()

    all_spots = refls["xyzcal.px"].as_numpy_array()[:, :2].astype("float32")
    pixels = img_set.get_raw_data(0)[0].as_numpy_array().astype("float32")

    # Estimate the radius once, from the full predicted set, then reuse it for
    # both the mask dilation and the integration window.
    radius = integration_radius
    if radius is None:
        radius = estimate_integration_radius(all_spots)
    logger.info(f"Image {img_num}: computed integration radius = {radius} px.")

    # Discard predictions whose integration window would overlap the detector
    # mask, dilated by that same radius.
    try:
        mask = np.array(img_set.get_mask(0)[0]).reshape(pixels.shape)
        x = np.floor(all_spots[:, 0]).astype(int)
        y = np.floor(all_spots[:, 1]).astype(int)
        sel = unmasked_prediction_selection(mask, x, y, radius, pixels.shape[1])
        refls = refls.select(flex.bool(sel))
        all_spots = all_spots[sel]
    except Exception as e:
        logger.warning(
            "Image %s: could not apply detector mask (%s); "
            "integrating all predictions.",
            img_num,
            e,
        )

    if len(all_spots) == 0:
        logger.warning(
            "Image %s: no predictions remain after masking. Skipping.", img_num
        )
        return flex.reflection_table()

    integrator = Integrator(
        pixels, all_spots, radius=radius, k=knn, isigi_cutoff=isigi_cutoff
    )
    try:
        integrator.fit(maxiter=maxiter)
    except RuntimeError as e:
        logger.warning("Image %s: %s Skipping.", img_num, e)
        return flex.reflection_table()

    # Update reflection data
    refls["intensity.sum.value"] = flex.double(integrator.intensity)
    refls["intensity.sum.variance"] = flex.double(np.square(integrator.uncertainty))
    refls["background.sum.value"] = flex.double(integrator.background.squeeze())
    # For Poisson noise, variance equals the mean background count
    refls["background.sum.variance"] = flex.double(integrator.background.squeeze())
    refls = refls.select(refls["intensity.sum.value"] != 0)
    refls = refls.select(refls["intensity.sum.variance"] > 0)
    logger.info(
        f"Image {img_num} took {time.time() - proctime} seconds to integrate {len(refls)} reflections."
    )
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

    # Print version information
    logger.info(laue_version())

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

    # Remove duplicate expt + refl data
    params.input.experiments = None
    params.input.reflections = None

    # Sanity checks
    if len(expts) == 0:
        parser.print_help()
        return

    # Get reflections and image data
    imagesets = expts.imagesets()
    ids = list(np.unique(preds["id"]).astype(np.int32))
    get_refls = partial(get_refls_image, preds)
    tables = list(map(get_refls, ids))
    if len(imagesets) != len(tables):
        logger.error(
            "Number of imagesets (%d) does not match the number of images with "
            "predictions (%d). Check that the experiment and reflection files "
            "correspond to the same dataset.",
            len(imagesets),
            len(tables),
        )
        return
    inputs = list(
        zip(
            imagesets,
            tables,
            repeat(params.isigi_cutoff),
            repeat(params.integration_radius),
            repeat(params.knn),
            repeat(params.maxiter),
        )
    )

    # Get initial time for process
    start_time = time.time()

    # Multiprocess integration
    num_processes = params.nproc
    logger.info("Starting integration.")
    if num_processes == 1:
        refls_arr = [integrate_image(*i) for i in inputs]
    else:
        with Pool(processes=num_processes) as pool:
            refls_arr = pool.starmap(integrate_image, inputs, chunksize=1)
    logger.info("Integration finished.")

    # Construct an integrated reflection table
    logger.info("Constructing reflection table")
    final_refls = flex.reflection_table()
    for refls in refls_arr:
        final_refls.extend(refls)
    refls = final_refls
    if len(refls) == 0:
        logger.error("No reflections were successfully integrated. Exiting.")
        return
    if params.output.reflections != None:
        logger.info(
            "Saving integrated reflection table to %s", params.output.reflections
        )
        refls.as_file(params.output.reflections)

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
