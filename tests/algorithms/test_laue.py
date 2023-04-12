import unittest

import gemmi
import numpy as np
import reciprocalspaceship as rs

from laue_dials.algorithms.laue import *


class TestLaueBase(unittest.TestCase):
    def setUp(self):
        # Set up a test instance of LaueBase
        self.s0 = np.array([1, 0, 0])
        self.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
        self.R = np.eye(3)
        self.lam_min = 0.5
        self.lam_max = 1.5
        self.dmin = 2
        self.spacegroup = "P1"
        self.LaueBase_obj = LaueBase(
            self.s0,
            self.cell,
            self.R,
            self.lam_min,
            self.lam_max,
            self.dmin,
            self.spacegroup,
        )

    def test_cell(self):
        self.assertEqual(self.LaueBase_obj.cell.a, 50)
        self.assertEqual(self.LaueBase_obj.cell.b, 50)
        self.assertEqual(self.LaueBase_obj.cell.c, 50)
        self.assertEqual(self.LaueBase_obj.cell.alpha, 90)
        self.assertEqual(self.LaueBase_obj.cell.beta, 90)
        self.assertEqual(self.LaueBase_obj.cell.gamma, 90)

    def test_spacegroup(self):
        self.assertEqual(self.LaueBase_obj.spacegroup.number, 1)

    def test_RB(self):
        RB_expected = np.array(self.LaueBase_obj.cell.fractionalization_matrix).T
        np.testing.assert_array_almost_equal(
            self.LaueBase_obj.RB, RB_expected, decimal=3
        )

    def test_s0(self):
        s0_expected = np.array([1, 0, 0])
        np.testing.assert_array_almost_equal(
            self.LaueBase_obj.s0, s0_expected, decimal=3
        )

    def test_Hall(self):
        # Check that all reflections have a d-value above the minimum
        d = self.LaueBase_obj.cell.calculate_d_array(self.LaueBase_obj.Hall)
        self.assertTrue(np.all(d >= self.dmin))

        # Check that all systematic absences are removed
        self.assertTrue(
            np.all(
                ~rs.utils.is_absent(
                    self.LaueBase_obj.Hall, self.LaueBase_obj.spacegroup
                )
            )
        )


class TestLaueAssigner(unittest.TestCase):
    def setUp(self):
        # Set up a test instance of LaueAssigner
        self.s0 = np.array([1, 0, 0])
        self.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
        self.R = np.eye(3)
        self.lam_min = 0.5
        self.lam_max = 1.5
        self.dmin = 2
        self.spacegroup = "P1"
        self.laue_assigner = LaueAssigner(
            self.s0,
            self.cell,
            self.R,
            self.lam_min,
            self.lam_max,
            self.dmin,
            self.spacegroup,
        )
