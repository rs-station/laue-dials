#!/usr/bin/env python
"""
This script predicts reflections for integration
"""
import logging

import gemmi
import libtbx.phil
import numpy as np
from dials.algorithms.spot_prediction import ray_intersection
from dials.array_family import flex
from dials.array_family.flex import reflection_table
from dials.util import log, show_mail_handle_errors
from dials.util.options import (ArgumentParser,
                                reflections_and_experiments_from_files)
from tqdm import trange

from laue_dials.algorithms.outliers import gen_kde

logger = logging.getLogger("laue-dials.command_line.predict")

help_message = """

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

wavelengths {
  lam_min = 0.95
    .type = float(value_min=0.1)
    .help = "Minimum wavelength for beam spectrum"
  lam_max = 1.15
    .type = float(value_min=0.2)
    .help = "Maximum wavelength for beam spectrum"
  }

reciprocal_grid {
  d_min = 1.4
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


def predict_spots(params, refls, expts):
    """
    A function for predicting spots given a geometry
    """
    from laue_dials.algorithms.laue import LauePredictor

    # Remove outliers
    print("Removing outliers")
    refls = refls.select(refls.get_flags(refls.flags.used_in_refinement))

    # Set up reflection table to store valid predictions
    final_preds = reflection_table()

    # Get experiment data from experiment objects
    print("Making predictions per image")
    for img_num in trange(len(expts.imagesets())):
        i = 0
        img = expts.imagesets()[img_num]
        experiment = expts[0]
        while True:  # Get first expt for this image
            experiment = expts[i]
            if experiment.imageset == img:
                break
            i = i + 1
        cryst = experiment.crystal
        spacegroup = gemmi.SpaceGroup(
            cryst.get_space_group().type().universal_hermann_mauguin_symbol()
        )

        # Get beam vector
        s0 = np.array(experiment.beam.get_s0())

        # Get unit cell params
        cell_params = cryst.get_unit_cell().parameters()
        cell = gemmi.UnitCell(*cell_params)

        # Get U matrix
        U = np.asarray(cryst.get_U()).reshape(3, 3)

        # Get observed centroids
        sub_refls = refls.select(refls["imageset_id"] == img_num)
        sub_refls["xyzobs.mm.value"].parts()[0].as_numpy_array()
        sub_refls["xyzobs.mm.value"].parts()[1].as_numpy_array()

        # Wavelengths per spot
        sub_refls["wavelength"].as_numpy_array()

        # Generate predictor object
        la = LauePredictor(
            s0,
            cell,
            U,
            params.wavelengths.lam_min,
            params.wavelengths.lam_max,
            params.reciprocal_grid.d_min,
            spacegroup,
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
        preds["imageset_id"] = flex.int([img_num] * len(preds))
        preds["s1"] = flex.vec3_double(s1)
        preds["phi"] = flex.double(np.zeros(len(s1)))  # Data are stills
        preds["wavelength"] = flex.double(new_lams)
        preds["rlp"] = flex.vec3_double(q_vecs)
        preds["miller_index"] = flex.miller_index(millers.astype("int").tolist())

        # Get which reflections intersect detector
        intersects = ray_intersection(experiment.detector, preds)
        preds = preds.select(intersects)
        new_lams = new_lams[intersects]

        # Generate a KDE
        _, _, kde = gen_kde(expts, refls)

        # Get predicted centroids
        x = preds["xyzcal.mm"].parts()[0].as_numpy_array()
        y = preds["xyzcal.mm"].parts()[1].as_numpy_array()

        # Get probability densities for predictions:
        rlps = preds["rlp"].as_numpy_array()
        norms = (np.linalg.norm(rlps, axis=1)) ** 2
        pred_data = [norms, new_lams]
        probs = kde.pdf(pred_data)

        # Cut off using log probabilities
        cutoff_log = params.cutoff_log_probability
        sel = np.log(probs) >= cutoff_log
        x[sel]
        y[sel]
        probs[sel]
        preds = preds.select(flex.bool(sel))

        # Append image predictions to dataset
        final_preds.extend(preds)

    # Return predicted refls
    return final_preds


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
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

    params, options = parser.parse_args(args=args, show_diff_phil=True)

    # Configure logging
    log.config(verbosity=options.verbose, logfile=params.output.log)

    # Log diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Load data
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    reflections = reflections[0]  # Get table out of list

    # Sanity checks
    if len(experiments) == 0:
        parser.print_help()
        return

    # Loop over experiments
    # Reindex data
    print(f"Predicting reflections")
    predicted_reflections = predict_spots(params, reflections, experiments)

    # Mark strong spots
    idpred, idstrong = predicted_reflections.match_by_hkle(reflections)
    strongs = np.zeros(len(predicted_reflections), dtype=int)
    strongs[idpred] = 1
    predicted_reflections["strong"] = flex.int(strongs)

    print(f"Assigning intensities")
    for i in trange(len(idstrong)):
        predicted_reflections["intensity.sum.value"][idpred[i]] = reflections[
            "intensity.sum.value"
        ][idstrong[i]]
        predicted_reflections["intensity.sum.variance"][idpred[i]] = reflections[
            "intensity.sum.variance"
        ][idstrong[i]]
        predicted_reflections["xyzobs.mm.value"][idpred[i]] = reflections[
            "xyzobs.mm.value"
        ][idstrong[i]]
        predicted_reflections["xyzobs.mm.variance"][idpred[i]] = reflections[
            "xyzobs.mm.variance"
        ][idstrong[i]]
        predicted_reflections["xyzobs.px.value"][idpred[i]] = reflections[
            "xyzobs.px.value"
        ][idstrong[i]]
        predicted_reflections["xyzobs.px.variance"][idpred[i]] = reflections[
            "xyzobs.px.variance"
        ][idstrong[i]]

    # Populate 'px' variety of predicted centroids
    # Based on flat rectangular detector
    x, y, z = predicted_reflections["xyzcal.mm"].parts()
    expt = experiments[0]  # assuming shared detector models
    x = x / expt.detector.to_dict()["panels"][0]["pixel_size"][0]
    y = y / expt.detector.to_dict()["panels"][0]["pixel_size"][1]
    predicted_reflections["xyzcal.px"] = flex.vec3_double(x, y, z)

    # Save reflections
    logger.info("Saving predicted reflections to %s", params.output.reflections)
    predicted_reflections.as_file(filename=params.output.reflections)


if __name__ == "__main__":
    run()
