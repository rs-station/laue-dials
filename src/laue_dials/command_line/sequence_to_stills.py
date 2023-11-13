"""
Split a sequence into a set of stills.

Example:

  laue.sequence_to_stills monochromatic.expt monochromatic.refl
"""

import logging
import sys
import time

import numpy as np
from dials.array_family import flex
from dials.array_family.flex import reflection_table
from dials.util import show_mail_handle_errors
from dials.util.options import (ArgumentParser,
                                reflections_and_experiments_from_files)
from dxtbx.model import MosaicCrystalSauter2014
from dxtbx.model.experiment_list import Experiment, ExperimentList
from libtbx.phil import parse
from scitbx import matrix
from tqdm import trange

logger = logging.getLogger("laue.command_line.sequence_to_stills")

# The phil scope
phil_scope = parse(
    """
output {
  experiments = stills.expt
    .type = str
    .help = "Filename for the experimental models that have been converted to stills"
  reflections = stills.refl
    .type = str
    .help = "Filename for the reflection tables with split shoeboxes (3D to 2D)"
  log = laue.sequence_to_stills.log
    .type = str
    .help = "Filename for the reflection tables with split shoeboxes (3D to 2D)"
  domain_size_ang = None
    .type = float
    .help = "Override for domain size. If None, use the crystal's domain size, if"
            "available"
  half_mosaicity_deg = None
    .type = float
    .help = "Override for mosaic angle. If None, use the crystal's mosaic angle, if"
            "available"
}
max_scan_points = None
  .type = int
  .expert_level = 2
  .help = Limit number of scan points
"""
)


def sequence_to_stills(experiments, reflections, params):
    assert len(reflections) == 1
    reflections = reflections[0]

    new_experiments = []
    new_reflections = []

    first_loop = True
    crystals = []
    imagesets = []

    print("Building experiment lists.")
    for expt_id, experiment in enumerate(experiments):
        # Generate crystals and imagesets for each scan point on first loop
        if first_loop:
            for i_scan_point in range(*experiment.scan.get_array_range()):
                # Get the goniometer setting matrix
                goniometer_setting_matrix = matrix.sqr(
                    experiment.goniometer.get_setting_rotation()
                )
                goniometer_axis = matrix.col(experiment.goniometer.get_rotation_axis())
                step = experiment.scan.get_oscillation()[1]

                # The A matrix is the goniometer setting matrix for this scan point
                # times the scan varying A matrix at this scan point. Note, the
                # goniometer setting matrix for scan point zero will be the identity
                # matrix and represents the beginning of the oscillation.
                # For stills, the A matrix needs to be positioned in the midpoint of an
                # oscillation step. Hence, here the goniometer setting matrixis rotated
                # by a further half oscillation step.
                A = (
                    goniometer_axis.axis_and_angle_as_r3_rotation_matrix(
                        angle=experiment.scan.get_angle_from_array_index(i_scan_point)
                        + (step / 2),
                        deg=True,
                    )
                    * goniometer_setting_matrix
                    * matrix.sqr(experiment.crystal.get_A())
                )
                crystal = MosaicCrystalSauter2014(experiment.crystal)
                crystal.set_A(A)

                # Copy in mosaic parameters if available
                if params.output.domain_size_ang is None and hasattr(
                    experiment.crystal, "get_domain_size_ang"
                ):
                    crystal.set_domain_size_ang(
                        experiment.crystal.get_domain_size_ang()
                    )
                elif params.output.domain_size_ang is not None:
                    crystal.set_domain_size_ang(params.output.domain_size_ang)

                if params.output.half_mosaicity_deg is None and hasattr(
                    experiment.crystal, "get_half_mosaicity_deg"
                ):
                    crystal.set_half_mosaicity_deg(
                        experiment.crystal.get_half_mosaicity_deg()
                    )
                elif params.output.half_mosaicity_deg is not None:
                    crystal.set_half_mosaicity_deg(params.output.half_mosaicity_deg)

                # Add to list of crystals
                crystals.append(crystal)

                # Add to list of imagesets
                imagesets.append(
                    experiment.imageset.as_imageset()[i_scan_point : i_scan_point + 1]
                )

                # Don't regenerate this list
                first_loop = False

        # Split experiment by scan points
        for i_scan_point in range(*experiment.scan.get_array_range()):
            new_experiment = Experiment(
                identifier=str(i_scan_point),
                detector=experiment.detector,
                beam=experiment.beam,
                crystal=crystals[i_scan_point],
                imageset=imagesets[i_scan_point],
            )
            new_experiments.append(new_experiment)

    # ----------------EXPERIMENTS CREATED---------------------------------
    print("Building reflection table.")
    for i_scan_point in trange(len(crystals)):
        # Get subset of reflections on this image
        _, _, _, _, z1, z2 = reflections["bbox"].parts()
        subrefls = reflections.select((i_scan_point >= z1) & (i_scan_point < z2))
        new_refls = subrefls.copy()
        centroids = np.asarray(subrefls["xyzobs.px.value"] - [0.0, 0.0, 0.5])
        centroids[:, 2] = i_scan_point
        new_refls["xyzobs.px.value"] = flex.vec3_double(centroids)
        new_refls["imageset_id"] = flex.int(len(new_refls), i_scan_point)
        x, y, _ = subrefls["xyzobs.mm.value"].parts()
        new_refls["xyzobs.mm.value"] = flex.vec3_double(
            x, y, flex.double(len(new_refls), 0)
        )
        new_refls["id"] = flex.int(len(new_refls), i_scan_point)
        new_reflections.append(new_refls)
    return (new_experiments, new_reflections)


@show_mail_handle_errors()
def run(args=None, phil=phil_scope):
    """
    Validate the arguments and load experiments/reflections for sequence_to_stills

    Arguments:
        args: The command line arguments to use. Defaults to sys.argv[1:]
        phil: The phil_scope. Defaults to the master phil_scope for this program
    """
    # The script usage
    usage = "usage: laue.sequence_to_stills [options] [param.phil] models.expt reflections.refl"

    # Create the parser
    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=__doc__,
    )
    params, options = parser.parse_args(args=args, show_diff_phil=True)

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

    # Try to load the models and data
    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        return

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    # Get initial time for process
    start_time = time.time()

    (new_experiments, new_reflections) = sequence_to_stills(
        experiments, reflections, params
    )
    # Write out the output experiments, reflections
    logger.info("Writing output data.")
    total_reflections = reflection_table()
    for i in trange(len(new_experiments)):
        total_reflections.extend(new_reflections[i])

    # Final logs
    logger.info("")
    logger.info(
        "Time Taken to Split into Stills = %f seconds", time.time() - start_time
    )
    ExperimentList(new_experiments).as_file("stills.expt")
    total_reflections.as_file("stills.refl")


if __name__ == "__main__":
    run()
