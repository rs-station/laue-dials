# This is a script to compare A matrices between DIALS and Precognition files

# NOTE: Python 3.5 or higher required. Matrix multiplication operator not defined in previous versions.

# NOTE: Build some functions to run this across entire datasets and histogram it. Throw it up on Github.

# In reciprocal space if not noted otherwise


from cctbx import crystal_orientation
from scitbx.matrix import sqr
from dxtbx.model import experiment_list
from dxtbx.model import crystal
import numpy as np
from gemmi import UnitCell, Fractional
import math
import glob
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator


def align_mats(crystal1, crystal2):  # Unfinished
    """
    Return A matrices in same alignment for both crystals

    Parameters
    ----------
    crystal1 : MosaicCrystalSauter2014, or similar object with a get_A method that returns a tuple with length 9
    crystal2 : MosaicCrystalSauter2014, or similar object with a get_A method that returns a tuple with length 9

    Returns
    -------
    U1 : numpy (3,3) matrix
        Aligned U matrix for first crystal
    U2 : numpy (3,3) matrix
        Aligned U matrix for second crystal
    """
    ref_orientation = crystal_orientation.crystal_orientation(
        crystal1.get_A(), crystal_orientation.basis_type.reciprocal
    )
    other_orientation = crystal_orientation.crystal_orientation(
        crystal2.get_A(), crystal_orientation.basis_type.reciprocal
    )
    # aligned_ori = crystal_orientation.crystal_orientation(other.get_A(), crystal_orientation.basis_type.reciprocal)
    try:
        op_align = sqr(
            other_orientation.best_similarity_transformation(
                ref_orientation, math.inf, 1
            )
        )
    except:
        print("Transformation past tolerance.")
        return np.reshape(np.asarray(crystal1.get_A()), (3, 3)), np.reshape(
            np.asarray(crystal2.get_A()), (3, 3)
        )
    aligned_ori = other_orientation.change_basis(op_align)

    U1 = np.reshape(np.asarray(ref_orientation.get_U_as_sqr()), (3, 3))
    U2 = np.reshape(np.asarray(aligned_ori.get_U_as_sqr()), (3, 3))
    return U1, U2


# Kevin's function for finding rotation angle between two A matrices
def rot_angle(crystal, other):
    """
    Calculate the angular rotation between two crystal objects.

    Parameters
    ----------
    crystal : MosaicCrystalSauter2014, or similar object with a get_A method that returns a tuple with length 9
    other   : MosaicCrystalSauter2014, or similar object with a get_A method that returns a tuple with length 9

    Returns
    -------
    angle : float
        The rotation angle between the two A matrices
    axis : float scitbx.matrix
        Axis of rotation between the two A-matrices. Exists as a 3x1 matrix.
    """
    ref_orientation = crystal_orientation.crystal_orientation(
        crystal.get_A(), crystal_orientation.basis_type.reciprocal
    )
    other_orientation = crystal_orientation.crystal_orientation(
        other.get_A(), crystal_orientation.basis_type.reciprocal
    )
    # aligned_ori = crystal_orientation.crystal_orientation(other.get_A(), crystal_orientation.basis_type.reciprocal)
    try:
        op_align = sqr(
            other_orientation.best_similarity_transformation(
                ref_orientation, math.inf, 1
            )
        )
    except:
        print("Transformation past tolerance.")
        return 0, 0
    aligned_ori = other_orientation.change_basis(op_align)

    U_ref = ref_orientation.get_U_as_sqr()
    U_other = aligned_ori.get_U_as_sqr()
    missetting_rot = U_other * U_ref.inverse()
    (
        angle,
        axis,
    ) = missetting_rot.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(
        deg=True
    )
    return angle, axis


def diff_cell_params(crystal, other):
    """
    Calculate the magnitude of the difference in unit cell
    parameters between two crystal objects.

    Parameters
    ----------

    Returns
    -------
    da : float

    db : float

    dc : float

    dalpha : float

    dbeta : float

    dgamma : float

    """

    cell1 = np.asarray(crystal.get_unit_cell().parameters(), dtype=float)
    cell2 = np.asarray(other.get_unit_cell().parameters(), dtype=float)
    cell_diff = abs(cell1 - cell2)
    da, db, dc, dalpha, dbeta, dgamma = cell_diff
    return da, db, dc, dalpha, dbeta, dgamma


def quot_cell_params(crystal, other):
    """
    Calculate the magnitude of the quotient in unit cell
    parameters between two crystal objects.

    Parameters
    ----------

    Returns
    -------
    qa : float

    qb : float

    qc : float

    qalpha : float

    qbeta : float

    qgamma : float

    """

    cell1 = np.asarray(crystal.get_unit_cell().parameters(), dtype=float)
    cell2 = np.asarray(other.get_unit_cell().parameters(), dtype=float)
    cell_quot = abs(cell1 / cell2)
    qa, qb, qc, qalpha, qbeta, qgamma = cell_quot
    return qa, qb, qc


# Filenames for DIALS and Precognition Inputs
# dials_filename_expt = "processing/monochromatic.expt"
dials_filename_expt = "processing/optimized.expt"
# dials_filename_expt = "processing/poly_refined.expt"
precog_filenames = sorted(
    glob.glob("precognition_files/input/e080_*.mccd.inp")
)  # "e080_001.mccd.inp"
plot_diffs = True

# Initialize  arrays
angle_diffs = np.zeros(len(precog_filenames))
angle_axes = np.ndarray(shape=(len(precog_filenames), 3), dtype=float)
da = np.zeros(len(precog_filenames))
db = np.zeros(len(precog_filenames))
dc = np.zeros(len(precog_filenames))
dalpha = np.zeros(len(precog_filenames))
dbeta = np.zeros(len(precog_filenames))
dgamma = np.zeros(len(precog_filenames))
qa = np.zeros(len(precog_filenames))
qb = np.zeros(len(precog_filenames))
qc = np.zeros(len(precog_filenames))

# DIALS part
# Load experiment list from JSON
# This should ONLY be ran on stills -- convert from sequence to stills
dials_experiment_model = experiment_list.ExperimentListFactory.from_json_file(
    dials_filename_expt
)
dials_crystal_model = dials_experiment_model.crystals()  # Plural from stills file


def get_rotation_matrix(axis, angle):
    u = axis
    sin, cos = np.sin(angle), np.cos(angle)
    return (
        cos * np.eye(3) + sin * np.cross(u, -np.eye(3)) + (1.0 - cos) * np.outer(u, u)
    )


for i in np.arange(len(precog_filenames) - 1):
    # Units for the following Precognition Parameters:
    #
    # Unit Cell: Angstroms and Degrees
    # Space Group: Integer denoting space group
    # Missetting Matrix: ???
    # Omega Polar Orientation: Degrees
    # Goniometer (omega, chi, phi): Degrees
    # Format: Detector Name
    # Detector Distance/Uncertainty: mm/??
    # Beam Center/Uncertainty: pixel/pixel
    # Pixel Size/Uncertainty: mm/mm
    # Swing Angle/Uncertainty: Degrees/??
    # Tilt Angle/Uncertainty: Degrees/??
    # Bulge Correction/Uncertainty: ??/??
    # Resolution Range: Angstroms
    # Wavelength Range: Angstroms

    # Precognition part
    # Parse Precognition file for unit cell params
    # NOTE: These conventions using Precognition conventions as listed in User Guide Release 5.0
    # Please check for correspondence with conventions used in DIALS
    # Do we need detector format? e.g. "RayonixMX340" ??
    # What is a bulge correction?
    # What is a swing angle?

    for line in open(precog_filenames[i]):
        rec = line.strip()
        if rec.startswith("Crystal"):
            unit_cell = rec.split()[1:7]  # a,b,c,alpha,beta,gamma
            space_group = rec.split()[7]
        if rec.startswith("Matrix"):
            missetting_matrix = rec.split()[1:10]
        if rec.startswith("Omega"):
            omega_polar_orientation = rec.split()[1:3]
        if rec.startswith("Goniometer"):
            goniometer = rec.split()[1:4]
        if rec.startswith("Format"):
            detector_format = rec.split()[1]
        if rec.startswith("Distance"):
            detector_distance = rec.split()[1]
            detector_distance_uncertainty = rec.split()[2]
        if rec.startswith("Center"):
            beam_center = rec.split()[1:3]
            beam_center_uncertainty = rec.split()[3:5]
        if rec.startswith("Pixel"):
            pixel_size = rec.split()[1:3]
            pixel_size_uncertainty = rec.split()[3:5]
        if rec.startswith("Swing"):
            swing_angles = rec.split()[1:3]
            swing_angles_uncertainty = rec.split()[
                3:5
            ]  # I'm guessing? User guide isn't explicit...
        if rec.startswith("Tilt"):
            detector_tilt_angles = rec.split()[1:3]
            detector_tilt_angles_uncertainty = rec.split()[
                3:5
            ]  # Also guessing this is uncertainty
        if rec.startswith("Bulge"):
            detector_bulge_correction = rec.split()[1:3]
            detector_bulge_correction_uncertainty = rec.split()[3:5]  # More guessing...
        #    if rec.startswith('Image'): rec.split()[1] # I think this only means something to Precognition's indexer
        if rec.startswith("Resolution"):
            resolution_range = rec.split()[1:3]
        if rec.startswith("Wavelength"):
            wavelength_range = rec.split()[1:3]

    # Convert angles to radians and prepare matrices for calculation
    # o1 \in [0,pi)
    # o2 \in [0,2pi)
    M = np.array(missetting_matrix, dtype=float).reshape((3, 3))
    o1 = np.deg2rad(float(omega_polar_orientation[0]))
    o2 = np.deg2rad(float(omega_polar_orientation[1]))
    gonio_phi = np.deg2rad(float(goniometer[2]))
    cell = UnitCell(*[float(i) for i in unit_cell])  # i is an iterating variable

    # Get rotation matrix from B to U
    R = get_rotation_matrix(np.array([0.0, 0.0, -1.0]), o1)
    temp = get_rotation_matrix(np.array([0.0, 1.0, 0.0]), o2)
    R = temp @ R
    temp = get_rotation_matrix(
        (R @ np.array([0.0, 1.0, 0.0])[:, None])[:, 0], gonio_phi
    )
    R = temp @ R

    # Create transposed orthogonalization matrix
    O = np.vstack(
        (
            cell.orthogonalize(Fractional(1.0, 0.0, 0.0)).tolist(),
            cell.orthogonalize(Fractional(0.0, 1.0, 0.0)).tolist(),
            cell.orthogonalize(Fractional(0.0, 0.0, 1.0)).tolist(),
        )
    ).T

    # Compute U, B, A matrices
    precog2mosflm = np.array(
        [[0, 0, 1], [0, -1, 0], [1, 0, 0]]
    )  # Change from Precognition to MOSFLM convention (this is a left operator)
    precog_A = precog2mosflm @ (R @ M @ np.linalg.inv(O))

    # precog_A = precog_U@precog_B # This is a lie
    # U is a properly oriented real-space crystallographic basis for frame in lab coordinate system
    # So is B

    # Create crystal object
    precog_crystal_model = crystal.CrystalFactory.from_mosflm_matrix(precog_A.flatten())

    #    # Plot DIALS + Precog A matrices per frame
    #    if plot_diffs:
    ##        dials_matrix = np.reshape(dials_crystal_model[i].get_A(), (3,3))
    ##        precog_matrix = np.reshape(precog_crystal_model.get_A(), (3,3))
    #        precog_matrix, dials_matrix = align_mats(precog_crystal_model, dials_crystal_model[i])
    #        fig, axs = plt.subplots(1,2)
    #        datamin = min(np.min(dials_matrix), np.min(precog_matrix))
    #        datamax = max(np.max(dials_matrix), np.max(precog_matrix))
    #        im1 = axs[0].matshow(dials_matrix,vmin=datamin,vmax=datamax)
    #        im2 = axs[1].matshow(precog_matrix,vmin=datamin,vmax=datamax)
    #        plt.suptitle(f'A Matrices For Frame {i}')
    #        axs[0].set_xlabel('DIALS')
    #        axs[1].set_xlabel('Precognition')
    #        fig.subplots_adjust(right=0.85)
    #        fig.colorbar(im2)
    #        plt.savefig(f'frame_A_matrix_diffs/frame_{i}')
    #        plt.close()

    # Compare A matrices and print information
    angle_diffs[i], angle_axes[i] = rot_angle(
        dials_crystal_model[i], precog_crystal_model
    )
    da[i], db[i], dc[i], dalpha[i], dbeta[i], dgamma[i] = diff_cell_params(
        dials_crystal_model[i], precog_crystal_model
    )
    qa[i], qb[i], qc[i] = quot_cell_params(dials_crystal_model[i], precog_crystal_model)

# -----------------------------------------------------------------------------
# Plot a histogram of the angle differences
n_bins = 30

frames = np.arange(len(precog_filenames))
ax = plt.figure().gca()
plt.plot(frames, angle_diffs)
plt.grid(axis="y", alpha=1)
plt.xlabel("Frames")
plt.ylabel("Angular differences (degrees)")
plt.title("DIALS vs Precognition Crystal Model Angular Offsets by Frame")
# plt.ylim([0,10])
# ax.yaxis.set_major_locator(MaxNLocator(integer=True)) # Crystals come in integers units
plt.savefig("diff_angles_per_frame.png")

ax = plt.figure().gca()
plt.plot(frames, angle_axes[:, 0])
plt.grid(axis="y", alpha=1)
plt.xlabel("Frames")
plt.ylabel("Rotation Axis x")
plt.title("DIALS vs Precognition Crystal Model Rotation Axis X-Components by Frame")
plt.savefig("axis_x_per_frame.png")

ax = plt.figure().gca()
plt.plot(frames, angle_axes[:, 1])
plt.grid(axis="y", alpha=1)
plt.xlabel("Frames")
plt.ylabel("Rotation Axis y")
plt.title("DIALS vs Precognition Crystal Model Rotation Axis Y-Components by Frame")
plt.savefig("axis_y_per_frame.png")

ax = plt.figure().gca()
plt.plot(frames, angle_axes[:, 2])
plt.grid(axis="y", alpha=1)
plt.xlabel("Frames")
plt.ylabel("Rotation Axis z")
plt.title("DIALS vs Precognition Crystal Model Rotation Axis Z-Components by Frame")
plt.savefig("axis_z_per_frame.png")


nrows = 2
ncols = 3
fig, axs = plt.subplots(2, 3, sharex=True)
fig.suptitle("Cell Parameter Differences by Frame")
fig.set_size_inches(12, 8)
axs[0, 0].plot(frames, da)
axs[0, 0].set_title(r"$\Delta a$")
axs[0, 1].plot(frames, db)
axs[0, 1].set_title(r"$\Delta b$")
axs[0, 2].plot(frames, dc)
axs[0, 2].set_title(r"$\Delta c$")
axs[1, 0].plot(frames, dalpha)
axs[1, 0].set_title(r"$\Delta \alpha$")
axs[1, 1].plot(frames, dbeta)
axs[1, 1].set_title(r"$\Delta \beta$")
axs[1, 2].plot(frames, dgamma)
axs[1, 2].set_title(r"$\Delta \gamma$")
for ax in axs.flat:
    ax.set(xlabel="Frames")
for ax in axs[0,]:
    ax.set(ylabel="Angstroms")
for ax in axs[1,]:
    ax.set(ylabel="Degrees")
# for ax in axs.flat:
#    ax.label_outer()
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig("cell_param_diffs_per_frame.png")

fig, axs = plt.subplots(1, 3, sharex=True)
fig.suptitle("Cell Parameter Quotients by Frame (DIALS/Precognition)")
fig.set_size_inches(12, 8)
axs[0].plot(frames, qa)
axs[0].set_title("a1/a2")
axs[1].plot(frames, qb)
axs[1].set_title("b1/b2")
axs[2].plot(frames, qc)
axs[2].set_title("c1/c2")
for ax in axs.flat:
    ax.set(xlabel="Frames")
# for ax in axs.flat:
#    ax.label_outer()
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig("cell_param_quots_per_frame.png")

## Plot changes in DIALS A matrix between frames
# for i in np.arange(len(precog_filenames)-2):
#    diff_A_matrix = np.reshape(np.asarray(dials_crystal_model[i+1].get_A()) - np.asarray(dials_crystal_model[i].get_A()), (3,3))
#    plt.matshow(diff_A_matrix, vmin=-0.01, vmax=0.01)
#    plt.title(f'DIALS A Matrix Difference Between Frame {i}/{i+1}')
#    plt.colorbar()
#    plt.savefig(f'frame_A_matrix_diffs/startframe_{i}')
#    plt.close()

# Plot angular changes in DIALS A matrix between frames
dials_angles = np.zeros(len(precog_filenames) - 2)
dials_axes = np.ndarray(shape=(len(precog_filenames) - 2, 3), dtype=float)
for i in np.arange(len(precog_filenames) - 2):
    dials_angles[i], dials_axes[i] = rot_angle(
        dials_crystal_model[i], dials_crystal_model[i + 1]
    )
ax = plt.figure().gca()
plt.plot(np.arange(len(dials_angles)), dials_angles)
plt.grid(axis="y", alpha=1)
plt.xlabel("Start Frames")
plt.ylabel("Angular difference to next frame (degrees)")
plt.title("DIALS Crystal Model Angular Changes by Frame")
plt.savefig("dials_diff_angles_per_frame.png")

# Plot angular changes in DIALS A matrix between frames
dials_angles = np.zeros(len(precog_filenames) - 1)
dials_axes = np.ndarray(shape=(len(precog_filenames) - 1, 3), dtype=float)
for i in np.arange(len(precog_filenames) - 1):
    dials_angles[i], dials_axes[i] = rot_angle(
        dials_crystal_model[0], dials_crystal_model[i]
    )
ax = plt.figure().gca()
plt.plot(np.arange(len(dials_angles)), dials_angles)
plt.grid(axis="y", alpha=1)
plt.xlabel("Frames")
plt.ylabel("Angular difference to from start (degrees)")
plt.title("DIALS Crystal Model Angular Changes by Frame Relative to First")
plt.savefig("dials_diff_angles_per_frame_from_first.png")
