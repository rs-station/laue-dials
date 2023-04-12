import gemmi
import numpy as np

from laue_dials.algorithms.diffgeo import *


def test_get_UB_matrices():
    # Define an example indexing matrix
    A = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])

    # Get the U and B matrices
    U, B = get_UB_matrices(A)

    # Verify that U is orthogonal
    assert np.allclose(np.dot(U, U.T), np.eye(3))

    # Verify that B is lower triangular
    assert np.allclose(np.tril(B), B)

    # Verify that A can be reconstructed from U and B
    assert np.allclose(A, U @ B)

    # Define another example indexing matrix
    A = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])

    # Get the U and B matrices
    U, B = get_UB_matrices(A)

    # Verify that U is orthogonal
    assert np.allclose(np.dot(U, U.T), np.eye(3))

    # Verify that B is lower triangular
    assert np.allclose(np.tril(B), B)

    # Verify that A can be reconstructed from U and B
    assert np.allclose(A, U @ B)


def test_normalize():
    # Define an example array
    A = np.array([[1, 2, 3], [4, 5, 6]])

    # Normalize the array
    normalized_A = normalize(A)

    # Verify that the L2 norm of each row is 1
    assert np.allclose(
        np.linalg.norm(normalized_A, axis=-1), np.ones(normalized_A.shape[0])
    )

    # Define another example array
    A = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])

    # Normalize the array
    normalized_A = normalize(A)

    # Verify that the L2 norm of each row is 1
    assert np.allclose(np.linalg.norm(normalized_A, axis=-1), np.array([0, 1, 1]))

    # Verify that the normalized array is equal to the original array divided by its L2 norm
    assert np.allclose(normalized_A, A)


def test_hkl2ray():
    # Define an example array of Miller indices
    hkl = np.array([[1, 2, 3], [4, 5, 6], [2, 4, 6]])

    # Compute the shortest member of the central ray for the Miller indices
    reduced_hkl = hkl2ray(hkl)

    # Verify that the reduced Miller indices are correct
    expected_reduced_hkl = np.array([[1, 2, 3], [4, 5, 6], [1, 2, 3]])
    assert np.allclose(reduced_hkl, expected_reduced_hkl)

    # Define an example array of wavelengths
    wavelength = np.array([1.0, 2.0, 3.0])

    # Compute the shortest member of the central ray and adjust the wavelengths accordingly
    reduced_hkl, reduced_wavelength = hkl2ray(hkl, wavelength=wavelength)

    # Verify that the reduced Miller indices are correct
    expected_reduced_hkl = np.array([[1, 2, 3], [4, 5, 6], [1, 2, 3]])
    assert np.allclose(reduced_hkl, expected_reduced_hkl)

    # Verify that the reduced wavelengths are correct
    expected_reduced_wavelength = np.array([1.0, 2.0, 6.0])
    assert np.allclose(reduced_wavelength, expected_reduced_wavelength)


def test_is_ray_equivalent():
    # Define two sets of Miller indices
    hkl1 = np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]])
    hkl2 = np.array([[2, 4, 8], [1, 2, 3], [6, 12, 18]])

    # Define expected output
    expected_output = np.array([False, True, True])

    # Check that output is accurate
    output = is_ray_equivalent(hkl1, hkl2)
    np.testing.assert_array_equal(output, expected_output)


def test_align_hkls():
    # Test data
    reference = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    target = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])
    spacegroup = gemmi.SpaceGroup(19)

    # Expected output
    expected_output = np.array([[1, -1, 0], [1, 0, 1], [0, -1, 1]])

    # Test the function
    aligned = align_hkls(reference, target, spacegroup)
    assert np.allclose(aligned, expected_output)


def test_orthogonalization():
    # Define unit cell
    a, b, c = 10, 15, 20
    alpha, beta, gamma = 90, 90, 90

    # Define expected orthogonalization matrix
    expected_output = np.array([[10, 0, 0], [0, 15, 0], [0, 0, 20]])

    # Test that generated matrix matches to rounding error
    np.testing.assert_array_almost_equal(
        orthogonalization(a, b, c, alpha, beta, gamma), expected_output
    )


def test_mat_to_rot_xyz():
    # Test a rotation matrix representing a 90 degree rotation about the z-axis
    R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    expected_angles = (0, 0, 90)

    # Test with degrees
    angles_deg = mat_to_rot_xyz(R, deg=True)
    np.testing.assert_almost_equal(angles_deg, expected_angles)

    # Test with radians
    angles_rad = mat_to_rot_xyz(R, deg=False)
    expected_angles_rad = np.radians(expected_angles)
    np.testing.assert_almost_equal(angles_rad, expected_angles_rad)


def test_rot_xyz_to_mat():
    # Test conversion of Euler angles to rotation matrix
    rot_x = np.pi / 2.0
    rot_y = np.pi / 3.0
    rot_z = np.pi / 4.0
    R_expected = np.array(
        [[0.354, 0.612, 0.707], [0.354, 0.612, -0.707], [-0.866, 0.500, 0.000]]
    )
    R_actual = rot_xyz_to_mat(rot_x, rot_y, rot_z, deg=False)
    np.testing.assert_allclose(R_actual, R_expected, atol=1e-3)

    # Test conversion of Euler angles in degrees to rotation matrix
    rot_x = 45.0
    rot_y = 60.0
    rot_z = 90.0
    R_expected = np.array(
        [[0.000, -0.707, 0.707], [0.500, 0.612, 0.612], [-0.866, 0.354, 0.353]]
    )
    R_actual = rot_xyz_to_mat(rot_x, rot_y, rot_z, deg=True)
    np.testing.assert_allclose(R_actual, R_expected, atol=1e-3)
