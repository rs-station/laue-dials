"""
Classes and functions for Laue-specific processes.
"""

import gemmi
import numpy as np
import reciprocalspaceship as rs
from dials.array_family import flex
from dxtbx.model import ExperimentList
from tqdm import tqdm, trange

from laue_dials.algorithms.diffgeo import hkl2ray


def store_wavelengths(expts, refls):
    """
    A function for storing wavelength data in beam objects in the reflection table.

    Parameters:
        expts (dxtbx.model.ExperimentList): ExperimentList to retrieve wavelengths from.
        refls (dials.array_family.flex.reflection_table): A reflection_table to store wavelengths in.

    Returns:
        new_refls (dials.array_family.flex.reflection_table): A reflection_table with updated wavelengths.
    """
    # Make new reflection_table
    new_refls = refls.copy()

    # Get experiment ID per reflection
    idx = new_refls["id"].as_numpy_array()

    # Initialize wavelength array
    lams = np.zeros(len(idx))

    # For each reflection, get corresponding beam wavelength
    for i in trange(len(idx)):
        lams[i] = expts[idx[i]].beam.get_wavelength()

    # Store wavelengths into reflection_table
    new_refls["wavelength"] = flex.double(lams)
    return new_refls


def remove_beam_models(expts, new_id):
    """
    A function for removing beam models no longer needed after refinement.

    Parameters:
        expts (dxtbx.model.ExperimentList): ExperimentList to remove beams from.
        new_id (int): ID to assign to the final experiment.

    Returns:
        new_expts (dxtbx.model.ExperimentList): An ExperimentList with an experiment per image, all sharing a beam.
    """
    # Instantiate new ExperimentList/reflection_table
    new_expts = ExperimentList()

    # Copy first experiment per image to new ExperimentList
    for img_num in trange(len(expts.imagesets())):
        i = 0
        img = expts.imagesets()[img_num]
        expt = expts[0]
        while True:
            expt = expts[i]
            if expt.imageset == img:
                break
            i = i + 1
        expt = expts[i]
        expt.identifier = str(new_id)  # Reset identifier to match image
        new_expts.append(expt)
    return new_expts


def gen_beam_models(expts, refls):
    """
    A function for generating beam models according to wavelengths in a reflection table.

    Parameters:
        expts (ExperimentList): ExperimentList to insert beams into.
        refls (dials.array_family.flex.reflection_table): A reflection_table containing wavelengths for beam models.

    Returns:
        new_expts (dxtbx.model.ExperimentList): ExperimentList with additional beam models and adjusted identifiers.
        refls (dials.array_family.flex.reflection_table): A reflection table with updated identifiers.
    """
    # Imports
    from copy import deepcopy

    from dials.algorithms.refinement.prediction.managed_predictors import \
        ExperimentsPredictorFactory

    # Instantiate new ExperimentList/reflection_table
    new_expts = ExperimentList()
    new_refls = refls.copy()

    # Unassign all reflections from experiments
    new_refls["id"] = flex.int([-1] * len(new_refls))

    # Generate beams per reflection
    exp_id = -1
    for i, refl in tqdm(enumerate(refls.rows())):
        try:
            # New beam per reflection
            expt = expts[refl["id"]]
            new_expt = expt
            new_expt.beam = deepcopy(expt.beam)
            new_expt.beam.set_wavelength(refl["wavelength"])  # Also updates s0
            exp_id = exp_id + 1  # Increment experiment ID
            new_expt.identifier = str(exp_id)
            new_expts.append(new_expt)
            new_refls["id"][i] = exp_id
        except:
            print("Warning: Unindexed strong spot has wavelength 0.")
            continue

    # Repredict centroids
    predictor = ExperimentsPredictorFactory.from_experiments(new_expts)
    new_refls = predictor(new_refls)
    return new_expts, new_refls


class LaueBase:
    """
    A base class to be extended for Laue procedures
    """

    def __init__(self, s0, cell, R, lam_min, lam_max, dmin, spacegroup="1"):
        """
        Parameters:
            s0 (np.ndarray): A 3-vector indicating the direction of the incoming beam wavevector.
            cell (Union[Tuple, gemmi.UnitCell]): A tuple or list of unit cell params (a, b, c, alpha, beta, gamma) or a gemmi.UnitCell object.
            R (np.ndarray): A 3x3 rotation matrix corresponding to the crystal orientation for the frame.
            lam_min (float): The lower end of the wavelength range of the beam.
            lam_max (float): The upper end of the wavelength range of the beam.
            dmin (float): The maximum resolution of the model.
            spacegroup (str): Anything that the gemmi.SpaceGroup constructor understands.
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
        """
        Calculates the product RB.

        Returns:
            np.ndarray: The product of the rotation matrix (R) and the fractionalization matrix (B).
        """
        return self.R @ self.B


class LaueAssigner(LaueBase):
    """
    An object to assign miller indices to a laue still
    """

    def __init__(self, s0, s1, cell, R, lam_min, lam_max, dmin, spacegroup="1"):
        """
        Parameters:
            s1 (np.ndarray): n x 3 array indicating the direction of the scattered beam wavevector.
        """
        super().__init__(s0, cell, R, lam_min, lam_max, dmin, spacegroup)

        self._s1 = s1 / np.linalg.norm(s1, axis=-1)[:, None]
        self._qobs = self._s1 - self.s0
        self._qpred = np.zeros_like(self._s1)
        self._H = np.zeros_like(self._s1)
        self._wav = np.zeros(len(self._H))
        self._harmonics = np.zeros(len(self._s1), dtype=bool)
        self._inliers = np.ones(len(self._s1), dtype=bool)

    @property
    def s1(self):
        """
        np.ndarray: The normalized direction of the scattered beam wavevector for inlying reflections.
        """
        return self._s1[self._inliers]

    @property
    def qobs(self):
        """
        np.ndarray: The observed q vectors for inlying reflections.
        """
        return self._qobs[self._inliers]

    @property
    def qpred(self):
        """
        np.ndarray: The predicted q vectors for inlying reflections.
        """
        return self._qpred[self._inliers]

    @property
    def H(self):
        """
        np.ndarray: The Miller indices for inlying reflections.
        """
        return self._H[self._inliers]

    @property
    def wav(self):
        """
        np.ndarray: The wavelengths associated with inlying reflections.
        """
        return self._wav[self._inliers]

    @property
    def harmonics(self):
        """
        np.ndarray: Boolean array indicating whether inlying reflections correspond to harmonic reflections.
        """
        return self._harmonics[self._inliers]

    # <-- setters that operate on the currently inlying set
    def set_qpred(self, qpred):
        """
        Set the predicted q vectors for inlying reflections.

        Parameters:
            qpred (np.ndarray): Predicted q vectors.
        """
        self._qpred[self._inliers] = qpred

    def set_H(self, H):
        """
        Set the Miller indices for inlying reflections.

        Parameters:
            H (np.ndarray): Miller indices.
        """
        self._H[self._inliers] = H
        self._H[~self._inliers] = 0.0

    def set_wav(self, wav):
        """
        Set the wavelengths for inlying reflections.

        Parameters:
            wav (np.ndarray): Wavelengths.
        """
        self._wav[self._inliers] = wav
        self._wav[~self._inliers] = 0.0

    def set_inliers(self, inliers):
        """
        Set the inliers for the current set of reflections.

        Parameters:
            inliers (np.ndarray): Boolean array indicating inliers.
        """
        self._inliers[self._inliers] = inliers

    def set_harmonics(self, harmonics):
        """
        Set whether inlying reflections correspond to harmonics.

        Parameters:
            harmonics (np.ndarray): Boolean array indicating harmonics.
        """
        self._harmonics[self._inliers] = harmonics

    # --> setters that operate on the currently inlying set

    def reset_inliers(self):
        """Reset all reflections as inliers."""
        self._inliers = np.ones(len(self._inliers), dtype=bool)

    def reject_outliers(self, nstd=10.0):
        """
        Update the list of inliers based on robust statistical measures.

        Parameters:
            nstd (float): Number of standard deviations from the median to consider as inliers.

        This method uses Minimum Covariance Determinant (MCD) to robustly estimate the covariance matrix
        of the concatenated observed and predicted q vectors. It then considers inliers as the points
        within a specified number of standard deviations from the median Mahalanobis distance.

        The list of inliers is updated accordingly.
        """
        from sklearn.covariance import MinCovDet

        X = np.concatenate((self.qobs, self.qpred * self.wav[:, None]), axis=-1)
        dist = MinCovDet().fit(X).dist_
        self.set_inliers(dist <= nstd**2.0)

    def assign(self):
        """
        Assign miller indices to the inlier reflections.

        This method updates the following attributes:
        - self.H: Miller indices associated with the assigned reflections.
        - self.wav: Wavelengths associated with the assigned reflections.
        - self.qpred: Predicted scattering vectors associated with the assigned reflections.
        - self.harmonics: Boolean array indicating whether the assigned reflections are harmonics.

        The assignment is performed by solving the linear sum assignment problem using scipy.optimize.linear_sum_assignment.
        The cost matrix is computed based on the angular distance between observed and predicted scattering vectors.

        The feasible set of reflections is determined from the reciprocal lattice points within the specified geometry.
        Reflections with duplicated scattering vectors are removed, and then the assignment is performed.

        This method is essential for updating the information of assigned reflections, allowing subsequent steps in
        Laue indexing procedures to utilize the assigned miller indices and associated parameters.
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

        ido, idx = linear_sum_assignment(cost)

        # Update appropriate variables
        H = Hall[idx]
        qpred = qall[idx]
        harmonics = harmonics[idx]

        # Set all attributes to match the current assignment
        self.set_H(H)
        self.set_qpred(qpred)
        self.set_harmonics(harmonics)

        # wav_pred = -2.*(self.s0 * qpred).sum(-1) / (qpred*qpred).sum(-1)
        with np.errstate(divide="ignore"):
            wav_obs = np.linalg.norm(self.qobs, axis=-1) / np.linalg.norm(
                self.qpred, axis=-1
            )
        self.set_wav(wav_obs)

    def update_rotation(self):
        """Update the rotation matrix (self.R) based on the inlying reflections."""
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
        super().__init__(s0, cell, R, lam_min, lam_max, dmin, spacegroup=spacegroup)

    def predict_s1(self, delete_harmonics=False):
        """
        Predicts scattering vectors and associated parameters for feasible Laue spots.

        Args:
            delete_harmonics (Optional[bool]): If True, removes harmonic reflections from the prediction.

        Returns:
            s1_pred (np.ndarray): Predicted feasible scattering vectors.
            lams (np.ndarray): Wavelengths associated with the predicted scattering vectors.
            qall (np.ndarray): Q vectors associated with the predicted scattering vectors.
            Hall (np.ndarray): Miller indices associated with the predicted scattering vectors.

        This method predicts scattering vectors for feasible spots given the current Laue experiment's geometry.
        The feasible set of reflections is determined within the specified geometry, and the scattering vectors are predicted.

        If `delete_harmonics` is True, harmonic reflections are removed from the prediction.
        Harmonics are identified based on duplicated scattering vectors in the reciprocal lattice.

        The returned arrays provide essential information about the predicted scattering vectors and associated parameters,
        allowing further analysis and utilization in Laue indexing procedures.
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
