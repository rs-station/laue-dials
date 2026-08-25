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
    detector_global_pixels,
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

# Minimum number of usable pixels a window must retain -- after the panel edge
# and the detector mask have been accounted for -- to be worth fitting.
MIN_WINDOW_PIXELS = 6


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


def panel_ids_for(refls):
    """
    Panel index of every reflection, defaulting to panel 0.

    Args:
        refls (dials.array_family.flex.reflection_table): Reflection table.

    Returns:
        np.ndarray: Integer panel index per reflection.
    """
    if "panel" in refls:
        return refls["panel"].as_numpy_array().astype(int)
    return np.zeros(len(refls), dtype=int)


def shared_frame_centroids(img_set, panel_ids, spots):
    """
    Express predicted centroids in a frame shared by every panel, in pixels.

    ``xyzcal.px`` is panel-local: every panel restarts at (0, 0). Distances
    between reflections on different panels are therefore meaningless, which
    matters twice over -- for the automatic integration radius, and for finding
    each reflection's nearest strong neighbours when pooling profiles. Mapping
    each centroid through its own panel into the laboratory frame and dividing
    by the pixel size restores a comparable set of coordinates.

    Args:
        img_set (dxtbx_imageset_ext.Imageset): Image set.
        panel_ids (np.ndarray): Panel index per reflection.
        spots (np.ndarray): (n, 2) panel-local centroids in pixels.

    Returns:
        np.ndarray or None: (n, 3) lab-frame centroids in pixel units, or None
        if the detector model is unavailable or the detector has one panel, in
        which case the panel-local coordinates are already shared.
    """
    try:
        detector = img_set.get_detector(0)
    except Exception:
        return None
    if detector is None or len(detector) < 2:
        return None
    try:
        pixel_size = float(np.mean(detector[0].get_pixel_size()))
        lab = np.array(
            [
                detector[int(p)].get_pixel_lab_coord((float(x), float(y)))
                for p, (x, y) in zip(panel_ids, spots)
            ]
        )
    except Exception as e:
        logger.warning(
            "Could not map centroids into the laboratory frame (%s); "
            "falling back to panel-local coordinates.",
            e,
        )
        return None
    return lab / pixel_size


def unmasked_selection_per_panel(panel_masks, panel_ids, spots, radius):
    """
    Apply the mask dilation of :func:`unmasked_prediction_selection` panel by
    panel, using each panel's own mask and row stride.

    Panels whose mask marks no bad pixels are left alone. The dilation's
    no-mask fallback buffers every detector edge by ``radius``, which on a
    segmented detector -- the LADI drum is 85 wedges only 50 px wide -- would
    throw away most of the predictions on grounds of geometry rather than data
    quality. Windows that overhang a panel edge are handled by the Integrator,
    which marks the pixels beyond the edge invalid.

    Args:
        panel_masks (sequence): One boolean mask per panel, True for good pixels.
        panel_ids (np.ndarray): Panel index per reflection.
        spots (np.ndarray): (n, 2) panel-local centroids in pixels.
        radius (int): Dilation radius in pixels.

    Returns:
        np.ndarray: Boolean array, True for predictions to keep.
    """
    keep = np.ones(len(spots), dtype=bool)
    x = np.floor(spots[:, 0]).astype(int)
    y = np.floor(spots[:, 1]).astype(int)
    for panel_id in np.unique(panel_ids):
        mask = np.asarray(panel_masks[int(panel_id)], dtype=bool)
        if mask.all():
            continue
        on_panel = panel_ids == panel_id
        keep[on_panel] = unmasked_prediction_selection(
            mask, x[on_panel], y[on_panel], radius, mask.shape[1]
        )
    return keep


def mtz_centroids(expts, refls):
    """
    Predicted centroids for the MTZ, on one grid per experiment.

    ``xyzcal.px`` is panel-local, so on a multi-panel detector every panel is
    written on top of every other one: metadata keys like ``xcal,ycal`` then
    describe 48 superimposed wedges rather than a position on the detector.
    Where the detector has more than one panel the centroids are laid out on a
    detector-wide grid instead. Single-panel detectors are untouched.

    Args:
        expts (dxtbx_model_ext.ExperimentList): Experiments.
        refls (dials.array_family.flex.reflection_table): Integrated reflections.

    Returns:
        tuple: (xcal, ycal), each an array of length len(refls), in pixels.
    """
    xyz = refls["xyzcal.px"].as_numpy_array()
    xcal, ycal = xyz[:, 0].copy(), xyz[:, 1].copy()

    ids = refls["id"].as_numpy_array()
    if "panel" in refls:
        panel_ids = refls["panel"].as_numpy_array().astype(int)
    else:
        panel_ids = np.zeros(len(refls), dtype=int)

    for expt_id in np.unique(ids):
        try:
            detector = expts[int(expt_id)].detector
        except (IndexError, AttributeError):
            logger.warning(
                "No detector model for experiment %s; leaving its centroids "
                "panel-local in the MTZ.",
                expt_id,
            )
            continue
        if detector is None or len(detector) < 2:
            continue
        on = ids == expt_id
        global_px = detector_global_pixels(detector, panel_ids[on], xyz[on, :2])
        if global_px is None:
            logger.warning(
                "Experiment %s: the panels have no common slow direction, so "
                "there is no detector-wide grid to write; leaving its "
                "centroids panel-local in the MTZ.",
                expt_id,
            )
            continue
        xcal[on] = global_px[:, 0]
        ycal[on] = global_px[:, 1]
        logger.info(
            "Experiment %s: wrote xcal/ycal on the detector-wide grid " "(%d panels).",
            expt_id,
            len(detector),
        )
    return xcal, ycal


def integrate_image(img_set, refls, isigi_cutoff, integration_radius, knn, maxiter):
    """
    Integrate predicted spots on an image.

    Every panel of the image is read, and each reflection is integrated from
    the panel it lies on. ``xyzcal.px`` is panel-local, so integrating the
    whole image against panel 0 -- as this function used to -- reads unrelated
    pixels for every reflection that is not on panel 0.

    The integration radius is estimated once from the full predicted set (or
    taken from ``integration_radius`` if supplied) in a frame shared by all
    panels, and reused both to dilate the detector mask -- discarding
    predictions whose integration window would overlap masked pixels -- and as
    the integration window itself.

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
    panel_ids = panel_ids_for(refls)

    # Pixel data for every panel, not just panel 0.
    panel_pixels = [
        panel.as_numpy_array().astype("float32") for panel in img_set.get_raw_data(0)
    ]
    if panel_ids.size and panel_ids.max() >= len(panel_pixels):
        logger.error(
            "Image %s: reflections reference panel %d but the image has only "
            "%d panel(s). Skipping.",
            img_num,
            int(panel_ids.max()),
            len(panel_pixels),
        )
        return flex.reflection_table()

    try:
        panel_masks = [
            np.array(mask).reshape(pixels.shape)
            for mask, pixels in zip(img_set.get_mask(0), panel_pixels)
        ]
    except Exception as e:
        logger.warning(
            "Image %s: could not read the detector mask (%s); "
            "treating every pixel as valid.",
            img_num,
            e,
        )
        panel_masks = None

    shared_coords = shared_frame_centroids(img_set, panel_ids, all_spots)
    neighbor_coords = all_spots if shared_coords is None else shared_coords

    # Estimate the radius once, from the full predicted set, then reuse it for
    # both the mask dilation and the integration window.
    radius = integration_radius
    if radius is None:
        radius = estimate_integration_radius(neighbor_coords)
    logger.info(f"Image {img_num}: computed integration radius = {radius} px.")

    # Discard predictions whose integration window would overlap the detector
    # mask, dilated by that same radius.
    if panel_masks is not None:
        try:
            sel = unmasked_selection_per_panel(
                panel_masks, panel_ids, all_spots, radius
            )
            refls = refls.select(flex.bool(sel.tolist()))
            all_spots = all_spots[sel]
            panel_ids = panel_ids[sel]
            neighbor_coords = neighbor_coords[sel]
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
        panel_pixels,
        all_spots,
        panel_ids=panel_ids,
        panel_masks=panel_masks,
        neighbor_coords=neighbor_coords,
        radius=radius,
        k=knn,
        isigi_cutoff=isigi_cutoff,
    )

    # A window can be left with too few usable pixels -- against a panel edge,
    # or over a masked region the dilation did not catch. Fitting a background,
    # an ellipse and an intensity to a handful of pixels is not meaningful, so
    # those reflections are dropped and the rest re-integrated without them.
    window_valid = getattr(integrator, "window_valid", None)
    if window_valid is not None:
        enough = window_valid.sum(-1) >= max(MIN_WINDOW_PIXELS, 0.25 * integrator.m)
        if not enough.all():
            logger.info(
                "Image %s: dropping %d prediction(s) with too few usable pixels.",
                img_num,
                int((~enough).sum()),
            )
            if not enough.any():
                return flex.reflection_table()
            refls = refls.select(flex.bool(enough.tolist()))
            integrator = Integrator(
                panel_pixels,
                all_spots[enough],
                panel_ids=panel_ids[enough],
                panel_masks=panel_masks,
                neighbor_coords=neighbor_coords[enough],
                radius=radius,
                k=knn,
                isigi_cutoff=isigi_cutoff,
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

    xcal, ycal = mtz_centroids(expts, refls)

    # Generate rs.DataSet to write to MTZ
    data = rs.DataSet(
        {
            "H": hkl.as_numpy_array()[:, 0].astype(np.int32),
            "K": hkl.as_numpy_array()[:, 1].astype(np.int32),
            "L": hkl.as_numpy_array()[:, 2].astype(np.int32),
            "BATCH": refls["id"].as_numpy_array() + 1,
            "I": refls["intensity.sum.value"].as_numpy_array(),
            "SIGI": refls["intensity.sum.variance"].as_numpy_array() ** 0.5,
            "xcal": xcal,
            "ycal": ycal,
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
