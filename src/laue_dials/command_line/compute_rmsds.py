#!/usr/bin/env python
"""
This script computes and plots RMSDs for a pair of DIALS experiment/reflection files.
"""

import logging
import sys

import gemmi
import libtbx.phil
import numpy as np
import pandas as pd
import reciprocalspaceship as rs
from cctbx import sgtbx
from dials.util import show_mail_handle_errors
from dials.util.options import (ArgumentParser,
                                reflections_and_experiments_from_files)
from matplotlib import pyplot as plt

from laue_dials.utils.version import laue_version

# Print laue-dials + DIALS versions
laue_version()

logger = logging.getLogger("laue-dials.command_line.compute_rmsds")

help_message = """
This program computes and plots the RMSDs between observed and predicted centroids in a reflection table.

Examples::

    laue.compute_rmsds [options] filename.expt filename.refl
"""

# Set the phil scope
phil_scope = libtbx.phil.parse(
    """
  show = True
    .type = bool
    .help = "Show the plot of centroid RMSDs per image."

  save = False
    .type = bool
    .help = "Save the plot of centroid RMSDs per image to a PNG file."

  output = "residuals.png"
    .type = str
    .help = "The filename for the generated plot."

  refined_only = False
    .type = bool
    .help = "Only compute refined spot RMSDs."

  log = 'laue.compute_rmsds.log'
    .type = str
    .help = "The log filename."
""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    """
    Compute and plot RMSDs for a pair of DIALS experiment/reflection files.

    Args:
        args (list): Command-line arguments.
        phil: The phil scope for the program.

    Returns:
        None
    """
    # Parse arguments
    usage = "laue.compute_rmsds [options] filename.expt filename.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=working_phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure logging
    console = logging.StreamHandler(sys.stdout)
    fh = logging.FileHandler(params.log, mode="w", encoding="utf-8")
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

    # Print help if no input
    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        exit()

    # Log diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Load data
    refls, expts = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    refls = refls[0]

    if params.refined_only:
        refls = refls.select(refls.get_flags(refls.flags.used_in_refinement))

    if len(refls) == 0:
        logger.info("No reflections in table after filtering.")
        return

    # Get data from reflection table
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
            "image": refls["id"].as_numpy_array() + 1,
            "xobs": refls["xyzobs.px.value"].as_numpy_array()[:, 0],
            "yobs": refls["xyzobs.px.value"].as_numpy_array()[:, 1],
            "xcal": refls["xyzcal.px"].as_numpy_array()[:, 0],
            "ycal": refls["xyzcal.px"].as_numpy_array()[:, 1],
        },
        cell=cell,
        spacegroup=spacegroup,
    ).infer_mtz_dtypes()

    logger.info(f"Total Number of Spots: {len(data)}.")

    # Calculate image residuals
    images = np.unique(data["image"])
    x_resids = data["xcal"] - data["xobs"]
    y_resids = data["ycal"] - data["yobs"]
    sqr_resids = x_resids**2 + y_resids**2
    mean_resids = np.zeros(len(images))
    for i, img_num in enumerate(images):
        sel = data["image"] == img_num
        mean_resids[i] = np.mean(sqr_resids[sel])
    rmsds = np.sqrt(mean_resids)

    resid_data = pd.DataFrame({"Image": images, "RMSD (px)": rmsds})

    logger.info(f"RMSDs per image: \n{resid_data}")

    # Get pixel size (assume square)
    # Not sure if this will be needed but I never remember
    # this incantation so leaving it here
    expts.detectors()[0].to_dict()["panels"][0]["pixel_size"][0]

    # Plot residuals
    fig = plt.figure()
    plt.scatter(images, rmsds)
    plt.title("Image RMSDs")
    plt.xlabel("Image #")
    plt.ylabel("RMSD (px)")
    if params.save:
        fig.savefig(params.output, format="png")
    if params.show:
        plt.show()


if __name__ == "__main__":
    run()
