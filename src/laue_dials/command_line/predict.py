#!/usr/bin/env python
"""
This script predicts reflections for integration
"""

import logging
import sys
import time
from itertools import repeat
from multiprocessing import Pool

import gemmi
import libtbx.phil
import numpy as np
from dials.algorithms.spot_prediction import ray_intersection
from dials.array_family import flex
from dials.array_family.flex import reflection_table
from dials.util import show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dxtbx.model import ExperimentList

from laue_dials.algorithms.outliers import gen_kde
from laue_dials.utils.version import laue_version

# Print laue-dials + DIALS versions
laue_version()

logger = logging.getLogger("laue-dials.command_line.predict")

help_message = """
This script predicts reflections for integration using a refined geometry experiment and reflection file.

This program takes a refined geometry experiment and reflection file, builds a
DIALS experiment list and reflection table, and predicts the feasible set of
reflections on the detector for integration using the refined geometry.

The output is a predicted reflection table (predicted.refl), which contains
the necessary information for the integrator to locate predicted spots
and integrate them.

Examples::

    laue.predict [options] poly_refined.expt poly_refined.refl
"""

# Set the phil scope
phil_scope = libtbx.phil.parse(
    """
output {
  reflections = 'predicted.refl'
    .type = str
    .help = "The output reflection table filename."

  log = 'laue.predict.log'
    .type = str
    .help = "The log filename."
  }

nproc = 1
  .type = int
  .help = Number of parallel processes to run

wavelengths {
  lam_min = None
    .type = float(value_min=0.1)
    .help = "Minimum wavelength for beam spectrum"
  lam_max = None
    .type = float(value_min=0.2)
    .help = "Maximum wavelength for beam spectrum"
  }

reciprocal_grid {
  d_min = None
    .type = float(value_min=0.1)
    .help = "Minimum d-spacing for reflecting planes"
  }

cutoff_log_probability = 0.
  .type = float
  .help = "The cutoff threshold for removing unlikely reflections"
""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


def predict_spots(lam_min, lam_max, d_min, refls, expts):
    """
    Predict spots given a geometry.

    Args:
        lam_min (float): Minimum wavelength for the beam spectrum.
        lam_max (float): Maximum wavelength for the beam spectrum.
        d_min (float): Minimum d-spacing for reflecting planes.
        refls (dials.array_family.flex.reflection_table): The reflection table.
        expts (dxtbx.model.experiment_list.ExperimentList): The experiment list.

    Returns:
        final_preds (dials.array_family.flex.reflection_table): Predicted reflection table.
    """
    from laue_dials.algorithms.laue import LauePredictor

    img_num = refls["id"][0]

    # Remove outliers
    refls = refls.select(refls.get_flags(refls.flags.used_in_refinement))

    # Set up reflection table to store valid predictions
    final_preds = reflection_table()

    try:
        # Get experiment data from experiment objects
        experiment = expts[0]
        cryst = experiment.crystal
        spacegroup = gemmi.SpaceGroup(
            cryst.get_space_group().type().universal_hermann_mauguin_symbol()
        )

        # Get mask
        mask = experiment.imageset.get_mask(0)[0]

        # Get beam vector
        s0 = np.array(experiment.beam.get_s0())

        # Get unit cell params
        cell_params = cryst.get_unit_cell().parameters()
        cell = gemmi.UnitCell(*cell_params)

        # Get U matrix
        U = np.asarray(cryst.get_U()).reshape(3, 3)

        # Get observed centroids
        sub_refls = refls.select(refls["id"] == img_num)

        # Generate predictor object
        logger.info(f"Predicting spots on image {img_num}.")
        la = LauePredictor(
            s0,
            cell,
            U,
            lam_min,
            lam_max,
            d_min,
            spacegroup=spacegroup,
        )

        # Predict spots
        s1, new_lams, q_vecs, millers = la.predict_s1()

        # Build new reflection table for predictions
        preds = reflection_table.empty_standard(len(s1))
        del preds["intensity.prf.value"]
        del preds["intensity.prf.variance"]
        del preds["lp"]
        del preds["profile_correlation"]

        # Populate needed columns
        preds["id"] = flex.int([int(experiment.identifier)] * len(preds))
        preds["imageset_id"] = flex.int([sub_refls[0]["imageset_id"]] * len(preds))
        preds["s1"] = flex.vec3_double(s1)
        preds["phi"] = flex.double(np.zeros(len(s1)))  # Data are stills
        preds["wavelength"] = flex.double(new_lams)
        preds["rlp"] = flex.vec3_double(q_vecs)
        preds["miller_index"] = flex.miller_index(millers.astype("int").tolist())

        # Get which reflections intersect detector
        intersects = ray_intersection(experiment.detector, preds)
        preds = preds.select(intersects)
        new_lams = new_lams[intersects]

        # Get predicted centroids
        x, y, _ = preds["xyzcal.mm"].parts()

        # Convert to pixel units
        px_size = experiment.detector.to_dict()["panels"][0]["pixel_size"]
        x = x / px_size[0]
        y = y / px_size[1]

        # Convert centroids to integer pixels
        x = np.asarray(flex.floor(x).iround())
        y = np.asarray(flex.floor(y).iround())

        # Remove predictions in masked areas
        img_row_size = experiment.detector.to_dict()["panels"][0]["image_size"][1]
        sel = np.full(len(x), True)
        for i in range(len(preds)):
            if not mask[x[i] + img_row_size * y[i]]:
                sel[i] = False
        preds = preds.select(flex.bool(sel))
        new_lams = new_lams[sel]
    except:
        logger.warning(
            f"WARNING: Could not predict reflections for experiment {img_num}. Image skipped."
        )
        return reflection_table()  # Return empty on failure

    # Append image predictions to dataset
    final_preds.extend(preds)

    # Return predicted refls
    return final_preds


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    """
    Run the prediction script.

    Args:
        args (list): Command-line arguments.
        phil: Working phil scope.
    """
    # Parse arguments
    usage = "laue.predict [options] poly_refined.expt poly_refined.refl"

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

    # Check for valid parameter values
    if params.reciprocal_grid.d_min == None:
        logger.info("Please provide a d_min.")
        return
    elif params.wavelengths.lam_min == None or params.wavelengths.lam_max == None:
        logger.info(
            "Please provide upper and lower boundaries for the wavelength spectrum."
        )
        return
    elif params.wavelengths.lam_min > params.wavelengths.lam_max:
        logger.info("Minimum wavelength cannot be greater than maximum wavelength.")
        return

    # Get initial time for process
    start_time = time.time()

    # Load data
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    reflections = reflections[0]  # Get table out of list
    reflections = reflections.select(
        reflections.get_flags(reflections.flags.used_in_refinement)
    )

    # Sanity checks
    if len(experiments) == 0:
        parser.print_help()
        return

    # Prepare parallel input
    ids = list(np.unique(reflections["id"]).astype(np.int32))
    expts_arr = []
    refls_arr = []
    for i in range(len(ids)):  # Split DIALS objects into lists
        expts_arr.append(ExperimentList([experiments[i]]))
        refls_arr.append(reflections.select(reflections["id"] == ids[i]))
    inputs = list(
        zip(
            repeat(params.wavelengths.lam_min),
            repeat(params.wavelengths.lam_max),
            repeat(params.reciprocal_grid.d_min),
            refls_arr,
            expts_arr,
        )
    )

    # Predict reflections
    logger.info(f"Predicting reflections")
    num_processes = params.nproc
    with Pool(processes=num_processes) as pool:
        output = pool.starmap(predict_spots, inputs, chunksize=1)
    logger.info(f"Finished predicting feasible spots.")

    # Convert output to single reflection table
    predicted_reflections = reflection_table()
    for table in output:
        predicted_reflections.extend(table)

    # Generate a KDE
    logger.info("Training KDE for resolution-dependent bandwidth.")
    _, _, kde = gen_kde(experiments, reflections)

    # Get probability densities for predictions:
    logger.info(f"Calculating prediction probabilities.")
    rlps = predicted_reflections["rlp"].as_numpy_array()
    norms = (np.linalg.norm(rlps, axis=1)) ** 2
    lams = predicted_reflections["wavelength"].as_numpy_array()
    pred_data = np.vstack([lams, norms])

    # Split array into chunks
    inputs = np.array_split(pred_data, num_processes, axis=1)

    # Multiprocess PDF estimation
    with Pool(processes=num_processes) as pool:
        prob_list = pool.map(kde.pdf, inputs)
    probs = np.concatenate(prob_list)

    # Cut off using log probabilities
    logger.info(f"Removing improbable reflections.")
    cutoff_log = params.cutoff_log_probability
    sel = np.log(probs) >= cutoff_log
    final_predictions = predicted_reflections.select(flex.bool(sel))

    # Mark strong spots
    logger.info("Marking strong predictions")
    idpred, idstrong = final_predictions.match_by_hkle(reflections)
    strongs = np.zeros(len(final_predictions), dtype=int)
    strongs[idpred] = 1
    final_predictions["strong"] = flex.int(strongs)

    logger.info(f"Assigning intensities")
    for i in range(len(idstrong)):
        final_predictions["intensity.sum.value"][idpred[i]] = reflections[
            "intensity.sum.value"
        ][idstrong[i]]
        final_predictions["intensity.sum.variance"][idpred[i]] = reflections[
            "intensity.sum.variance"
        ][idstrong[i]]
        final_predictions["xyzobs.mm.value"][idpred[i]] = reflections[
            "xyzobs.mm.value"
        ][idstrong[i]]
        final_predictions["xyzobs.mm.variance"][idpred[i]] = reflections[
            "xyzobs.mm.variance"
        ][idstrong[i]]
        final_predictions["xyzobs.px.value"][idpred[i]] = reflections[
            "xyzobs.px.value"
        ][idstrong[i]]
        final_predictions["xyzobs.px.variance"][idpred[i]] = reflections[
            "xyzobs.px.variance"
        ][idstrong[i]]

    # Populate 'px' variety of predicted centroids
    # Based on flat rectangular detector
    x, y, z = final_predictions["xyzcal.mm"].parts()
    expt = experiments[0]  # assuming shared detector models
    x = x / expt.detector.to_dict()["panels"][0]["pixel_size"][0]
    y = y / expt.detector.to_dict()["panels"][0]["pixel_size"][1]
    final_predictions["xyzcal.px"] = flex.vec3_double(x, y, z)

    # Save reflections
    logger.info("Saving predicted reflections to %s", params.output.reflections)
    final_predictions.as_file(filename=params.output.reflections)

    # Final logs
    logger.info("")
    logger.info(
        "Time Taken for Total Processing = %f seconds", time.time() - start_time
    )


if __name__ == "__main__":
    run()
