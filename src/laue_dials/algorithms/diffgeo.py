"""
This file contains useful functions for diffraction geometry calculations
"""

import numpy as np
import reciprocalspaceship as rs


def get_UB_matrices(A):
    """
    Convert the reciprocal space indexing matrix, A into the product of
    an orthogonal matrix U and a lower triangular matrix B

    Parameters
    ----------
    A : np.array
        A 3x3 indexing matrix such that the scattering vector `S1-S0=Q=A@h`

    Returns
    -------
    U : np.array
        An orthogonal 3x3 matrix
    B : np.array
        A 3x3 lower triangular matrix
    """
    from scipy.linalg import rq

    # R is a "right" or upper triangular matrix
    # Q is an orthogonal (rotation) matrix
    # A = Q.T@R.T
    R, Q = rq(A.T)
    return Q.T, R.T


def normalize(A):
    """
    Normalize the last dimension of an array by dividing by its L2 norm

    Parameters
    ----------
    A : np.array
        A non-normal matrix

    Returns
    -------
    A : np.array
        A normalized matrix
    """
    normalized_A = A / np.linalg.norm(A, axis=-1)[..., None]

    # In case of magnitude-0 rows:
    normalized_A[np.isnan(normalized_A)] = 0
    return normalized_A


def hkl2ray(hkl, wavelength=None):
    """
    Convert a miller index to the shortest member of its central ray.
    Optionally, adjust its wavelength accordingly.

    Parameters
    ----------
    hkl : array
        `n x 3` array of miller indices. the dtype must be interpretable as an integer
    wavelength : array (optional)
        length `n` array of wavelengths corresponding to each miller index

    Returns
    -------
    reduced_hkl : array
        The miller index of the shortest vector on the same central ray as hkl
    reduced_wavelength : array (optional)
        The wavelengths corresponding to reduced_hkl
    """
    gcd = np.gcd.reduce(hkl.astype(int), axis=-1)
    if wavelength is not None:
        return hkl / gcd[..., None], wavelength * gcd
    else:
        return hkl / gcd[..., None]


def is_ray_equivalent(hkl1, hkl2):
    """
    Test for equivalency between two miller indices in a Laue experiment. Returns a boolean array for each of the `n` hkls in `hkl{1,2}`.

    Parameters
    ----------
    hkl1 : array
        `n x 3` array of miller indices. the dtype must be interpretable as an integer
    hkl2 : array
        `n x 3` array of miller indices. the dtype must be interpretable as an integer

    Returns
    -------
    equal_hkl : array
        Boolean array for hkl equivalency by index from original arrays
    """
    return np.all(np.isclose(hkl2ray(hkl1), hkl2ray(hkl2)), axis=1)


def align_hkls(reference, target, spacegroup, anomalous=True):
    """
    Use the point group operators in `spacegroup` to align target hkls to reference.

    Parameters
    ----------
    reference : array
        n x 3 array of miller indices
    target : array
        n x 3 array of miller indices
    spacegroup : gemmi.Spacegroup
        The space group of reference/target.
    anomalous : bool(optional)
        If true, test Friedel symmetries too.

    Returns
    -------
    aligned : array
        n x 3 array of miller indices equivalent to target
    """
    aligned = target
    cc = -1.0
    for op in spacegroup.operations():
        aligned_ = rs.utils.apply_to_hkl(target, op)
        cc_mat = np.corrcoef(aligned_.T, reference.T)
        cc_ = np.trace(cc_mat[:3, 3:])
        if cc_ > cc:
            aligned = aligned_
            cc = cc_
    if anomalous:
        for op in spacegroup.operations():
            aligned_ = -rs.utils.apply_to_hkl(target, op)
            cc_mat = np.corrcoef(aligned_.T, reference.T)
            cc_ = np.trace(cc_mat[:3, 3:])
            if cc_ > cc:
                aligned = aligned_
                cc = cc_
    return aligned


def orthogonalization(a, b, c, alpha, beta, gamma):
    """
    Compute the orthogonalization matrix from cell params

    Parameters
    ----------
    a,b,c : floats
        Floats representing the magnitudes of the three unit cell axes
    alpha, beta, gamma : floats
        Floats representing the three unit cell angles

    Returns
    -------
    orthogonalization_matrix : array
        3 x 3 orthogonalization matrix for unit cell
    """
    alpha, beta, gamma = (
        np.pi * alpha / 180.0,
        np.pi * beta / 180.0,
        np.pi * gamma / 180.0,
    )
    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    V = (
        a
        * b
        * c
        * np.sqrt(
            1.0 - cosa * cosa - cosb * cosb - cosg * cosg + 2 * cosa * cosb * cosg
        )
    )
    return np.array(
        [
            [a, b * cosg, c * cosb],
            [0.0, b * sing, c * (cosa - cosb * cosg) / sing],
            [0.0, 0.0, V / a / b / sing],
        ]
    )


def mat_to_rot_xyz(R, deg=True):
    """
    Decompose a rotation matrix into euler angles

    Parameters
    ----------
    R : array
        3 x 3 rotation matrix

    Returns
    -------
    rot_x, rot_y, rot_z : floats
        Euler angles to generate rotation matrix
    """
    if R[2, 0] < 1:
        if R[2, 0] > -1:
            rot_y = np.arcsin(-R[2, 0])
            rot_z = np.arctan2(R[1, 0], R[0, 0])
            rot_x = np.arctan2(R[2, 1], R[2, 2])
        else:
            rot_y = np.pi / 2.0
            rot_z = -np.arctan2(-R[1, 2], R[1, 1])
            rot_x = 0.0
    else:
        rot_y = np.pi / 2.0
        rot_z = np.arctan2(-R[1, 2], R[1, 1])
        rot_x = 0.0
    if deg:
        rot_x = np.rad2deg(rot_x)
        rot_y = np.rad2deg(rot_y)
        rot_z = np.rad2deg(rot_z)
    return rot_x, rot_y, rot_z


def rot_xyz_to_mat(rot_x, rot_y, rot_z, deg=True):
    """
    Convert euler angles into a rotation matrix

    Parameters
    ----------
    rot_x, rot_y, rot_z : floats
        Euler angles to generate rotation matrix

    Returns
    -------
    R : array
        3 x 3 rotation matrix
    """
    if deg:
        rot_x = np.deg2rad(rot_x)
        rot_y = np.deg2rad(rot_y)
        rot_z = np.deg2rad(rot_z)
    cx, sx = np.cos(rot_x), np.sin(rot_x)
    cy, sy = np.cos(rot_y), np.sin(rot_y)
    cz, sz = np.cos(rot_z), np.sin(rot_z)

    return np.array(
        [
            [cy * cz, cz * sx * sy - cx * sz, cx * cz * sy + sx * sz],
            [cy * sz, cx * cz + sx * sy * sz, -cz * sx + cx * sy * sz],
            [-sy, cy * sx, cx * cy],
        ]
    )
