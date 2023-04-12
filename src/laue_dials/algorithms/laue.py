"""
Classes and functions for Laue-specific processes.
"""

import numpy as np
import reciprocalspaceship as rs
import gemmi
from laue_dials.algorithms.diffgeo import hkl2ray


class LaueBase:
    """
    A base class to be extended for Laue procedures
    """

    def __init__(self, s0, cell, R, lam_min, lam_max, dmin, spacegroup="1"):
        """
        Parameters
        ----------
        s0 : array
            a 3 vector indicating the direction of the incoming beam wavevector.
            This can be any length, it will be unit normalized in the constructor.
        cell : iterable or gemmi.UnitCell
            A tuple or list of unit cell params (a, b, c, alpha, beta, gamma) or a gemmi.UnitCell object
        R : array
            A 3x3 rotation matrix corresponding to the crystal orientation for the frame.
        lam_min : float
            The lower end of the wavelength range of the beam.
        lam_max : float
            The upper end of the wavelength range of the beam.
        dmin : float
            The maximum resolution of the model
        spacegroup : gemmi.SpaceGroup (optional)
            Anything that the gemmi.SpaceGroup constructor understands.
        """
        if not isinstance(cell, gemmi.UnitCell):
            cell = gemmi.UnitCell(*cell)
        self.cell = cell

        if not isinstance(spacegroup, gemmi.SpaceGroup):
            spacegroup = gemmi.SpaceGroup(spacegroup)
        self.spacegroup = spacegroup

        self.R = R
        self.lam_min = lam_min
        self.lam_max = lam_max
        self.dmin = dmin

        self.s0 = s0 / np.linalg.norm(s0)
        self.B = np.array(self.cell.fractionalization_matrix).T

        # Initialize the full reciprocal grid
        hmax, kmax, lmax = self.cell.get_hkl_limits(dmin)
        Hall = (
            np.mgrid[
                -hmax : hmax + 1 : 1.0,
                -kmax : kmax + 1 : 1.0,
                -lmax : lmax + 1 : 1.0,
            ]
            .reshape((3, -1))
            .T
        )
        Hall = Hall[np.any(Hall != 0, axis=1)]
        d = cell.calculate_d_array(Hall)
        Hall = Hall[d >= dmin]

        # Remove any systematic absences in the space group
        Hall = Hall[~rs.utils.is_absent(Hall, self.spacegroup)]
        self.Hall = Hall

    @property
    def RB(self):
        return self.R @ self.B


class LaueAssigner(LaueBase):
    """
    An object to assign miller indices to a laue still
    """

    def __init__(self, s0, s1, cell, R, lam_min, lam_max, dmin, spacegroup="1"):
        """
        Parameters
        ----------
        s1 : array
            n x 3 array indicating the direction of the scatterd beam wavevector.
            This can be any length, it will be unit normalized in the constructor.
        """
        super().__init__(s0, cell, R, lam_min, lam_max, dmin, spacegroup="1")

        self.ewald_offset = None

        self._s1 = s1 / np.linalg.norm(s1, axis=-1)[:, None]
        self._qobs = self._s1 - self.s0
        self._qpred = np.zeros_like(self._s1)
        self._H = np.zeros_like(self._s1)
        self._wav = np.zeros(len(self._H))
        self._harmonics = np.zeros(len(self._s1), dtype=bool)
        self._inliers = np.ones(len(self._s1), dtype=bool)

    @property
    def s1(self):
        return self._s1[self._inliers]

    @property
    def qobs(self):
        return self._qobs[self._inliers]

    @property
    def qpred(self):
        return self._qpred[self._inliers]

    @property
    def H(self):
        return self._H[self._inliers]

    @property
    def wav(self):
        return self._wav[self._inliers]

    @property
    def harmonics(self):
        return self._harmonics[self._inliers]

    # <-- setters that operate on the currently inlying set
    def set_qpred(self, qpred):
        self._qpred[self._inliers] = qpred

    def set_H(self, H):
        self._H[self._inliers] = H
        self._H[~self._inliers] = 0.0

    def set_wav(self, wav):
        self._wav[self._inliers] = wav
        self._wav[~self._inliers] = 0.0

    def set_inliers(self, inliers):
        self._inliers[self._inliers] = inliers

    def set_harmonics(self, harmonics):
        self._harmonics[self._inliers] = harmonics

    # --> setters that operate on the currently inlying set

    def reset_inliers(self):
        self._inliers = np.ones(len(self._inliers), dtype=bool)

    def reject_outliers(self, nstd=10.0):
        """update the list of inliers"""
        from sklearn.covariance import MinCovDet

        X = np.concatenate((self.qobs, self.qpred * self.wav[:, None]), axis=-1)
        dist = MinCovDet().fit(X).dist_
        self.set_inliers(dist <= nstd**2.0)

    def assign(self):
        """
        Assign miller indices to the inlier reflections
        This method will update:
            self.H     -- miller indices
            self.wav   -- wavelengths
            self.qpred -- predicted scattering vector
        """
        # Generate the feasible set of reflections from the current geometry
        Hall = self.Hall
        qall = (self.RB @ Hall.T).T
        feasible = (
            np.linalg.norm(qall + self.s0 / self.lam_min, axis=-1) < 1 / self.lam_min
        ) & (np.linalg.norm(qall + self.s0 / self.lam_max, axis=-1) > 1 / self.lam_max)
        Hall = Hall[feasible]
        qall = qall[feasible]

        # Keep track of harmonics in the feasible set
        Raypred = hkl2ray(Hall)
        _, idx, counts = np.unique(
            Raypred, return_index=True, return_counts=True, axis=0
        )
        harmonics = counts > 1

        # Remove harmonics from the feasible set
        Hall = Hall[idx]
        qall = qall[idx]

        dmat = rs.utils.angle_between(self.qobs[..., None, :], qall[None, ..., :])
        cost = dmat

        from scipy.optimize import linear_sum_assignment

        _, idx = linear_sum_assignment(cost)
        H = Hall[idx]
        qpred = qall[idx]
        harmonics = harmonics[idx]

        # Set all attributes to match the current assignment
        self.set_H(H)
        self.set_qpred(qpred)
        self.set_harmonics(harmonics)

        # wav_pred = -2.*(self.s0 * qpred).sum(-1) / (qpred*qpred).sum(-1)
        wav_obs = np.linalg.norm(self.qobs, axis=-1) / np.linalg.norm(
            self.qpred, axis=-1
        )
        self.set_wav(wav_obs)

    def update_rotation(self):
        """Update the rotation matrix (self.R) based on the inlying refls"""
        from scipy.linalg import orthogonal_procrustes

        misset, _ = orthogonal_procrustes(
            self.qobs,
            self.qpred * self.wav[:, None],
        )
        self.R = misset @ self.R


class LauePredictor(LaueBase):
    """
    An object to predict spots given a Laue experiment.
    """

    def __init__(self, s0, cell, R, lam_min, lam_max, dmin, spacegroup="1"):
        """
        Calls parent LaueBase init
        """
        super().__init__(s0, cell, R, lam_min, lam_max, dmin, spacegroup="1")

    def predict_s1(self, delete_harmonics=False):
        """
        Predicts all s1 vectors for all feasible spots given some resolution-dependent bandwidth
        This method provides:
            s1_pred -- predicted feasible s1 vectors
            lams -- the wavelengths (in Angstroms) associated with these s1 vectors
            qall -- the q vectors associated with these s1 vectors
        """
        # Generate the feasible set of reflections from the current geometry
        Hall = self.Hall
        qall = (self.RB @ Hall.T).T
        feasible = (
            np.linalg.norm(qall + self.s0 / self.lam_min, axis=-1) < 1 / self.lam_min
        ) & (np.linalg.norm(qall + self.s0 / self.lam_max, axis=-1) > 1 / self.lam_max)
        Hall = Hall[feasible]
        qall = qall[feasible]

        # Remove harmonics from the feasible set
        Raypred = hkl2ray(Hall)
        _, idx, counts = np.unique(
            Raypred, return_index=True, return_counts=True, axis=0
        )
        Hall = Hall[idx]  # Remove duplicates
        qall = qall[idx]
        if delete_harmonics:
            idx = counts > 1  # Remove last harmonic
            Hall = Hall[idx]
            qall = qall[idx]

        # For each q, find the wavelength of the Ewald sphere it lies on
        lams = -2.0 * (self.s0 * qall).sum(-1) / (qall * qall).sum(-1)

        # Using this wavelength per q, generate s1 vectors
        s0 = self.s0[None, :] / lams[:, None]
        s1_pred = qall + s0

        # Write s1 predictions
        return s1_pred, lams, qall, Hall
