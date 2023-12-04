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

def norm2(array, axis=-1, keepdims=False):
    """Faster version of np.linalg.norm for 3d coords and L2 norm"""
    a2 = np.square(array)
    out = 0.
    for vec in np.split(a2, a2.shape[axis], axis=axis):
        out += vec
    out = np.sqrt(out)
    if not keepdims:
        out = np.squeeze(out, axis=axis)
    return out

def store_wavelengths(expts, refls):
    """
    A function for storing wavelength data in beam objects in the reflection table
        Parameters
        ----------
        expts : ExperimentList
            ExperimentList to retrieve wavelengths from
        refls : reflection_table
            A reflection_table to store wavelengths in

        Returns
        -------
        new_refls : reflection_table
            A reflection_table with updated wavelengths
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
    A function for removing beam models no longer needed after refinement
        Parameters
        ----------
        expts : ExperimentList
            ExperimentList to remove beams from

        new_id : Int
            ID to assign to final experiment

        Returns
        -------
        new_expts : ExperimentList
            An ExperimentList with an experiment per image all sharing a beam
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
    A function for generating beam models according to wavelengths in a reflection table
        Parameters
        ----------
        expts : ExperimentList
            ExperimentList to insert beams into
        refls : reflection_table
            A reflection_table containing wavelengths for beam models

        Returns
        -------
        new_expts : ExperimentList
            ExperimentList with additional beam models and adjusted identifiers
        refls : reflection_table
            A reflection_table with updated identifiers
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

    def predict_s1(self, delete_harmonics=False):
        """
        Predicts all s1 vectors for all feasible spots given some resolution-dependent bandwidth
        This method provides:
            s1_pred -- predicted feasible s1 vectors
            lams -- the wavelengths (in Angstroms) associated with these s1 vectors
            qall -- the q vectors associated with these s1 vectors
            Hall -- the miller indices associated with s1 vectors
        """
        # Generate the feasible set of reflections from the current geometry
        Hall = self.Hall
        qall = (self.RB @ Hall.T).T
        feasible = (
            norm2(qall + self.s0 / self.lam_min, axis=-1) < 1 / self.lam_min
        ) & (norm2(qall + self.s0 / self.lam_max, axis=-1) > 1 / self.lam_max)
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
        super().__init__(s0, cell, R, lam_min, lam_max, dmin, spacegroup)

        self.s1 = s1 / norm2(s1, axis=-1, keepdims=True)
        self.qobs = self.s1 - self.s0
        self.qpred = np.zeros_like(self.s1)
        self.H = np.zeros_like(self.s1)
        self.wav = np.zeros(len(self.H))
        self.harmonics = np.zeros(len(self.s1), dtype=bool)

    def plot_preds(self):
        from matplotlib import pyplot as plt
        xyobs = self.s1[:,:2] / self.wav[:,None]
        xycal = self.predict_s1()[0][:,:2]
        plt.plot(xyobs[:,0], xyobs[:,1], 'k.', label='Observed')
        plt.scatter(xycal[:,0], xycal[:,1], s=80, facecolors='none', edgecolors='r', label='Calculated')
        plt.legend()

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
            norm2(qall + self.s0 / self.lam_min, axis=-1) < 1 / self.lam_min
        ) & (norm2(qall + self.s0 / self.lam_max, axis=-1) > 1 / self.lam_max)
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
        self.harmonics = harmonics

        # Set all attributes to match the current assignment
        self.H = H
        self.qpred = qpred

        # wav_pred = -2.*(self.s0 * qpred).sum(-1) / (qpred*qpred).sum(-1)
        with np.errstate(divide="ignore"):
            wav_obs = norm2(self.qobs, axis=-1) / norm2(
                self.qpred, axis=-1
            )
        self.wav = wav_obs

    def update_rotation(self):
        """Update the rotation matrix (self.R) based on the inlying refls"""
        from scipy.linalg import orthogonal_procrustes

        misset, _ = orthogonal_procrustes(
            self.qobs,
            self.qpred * self.wav[:, None],
        )
        self.R = misset @ self.R

    def refine_cell(self):
        """ Refine the cell axes keeping b fixed """
        from scipy.optimize import minimize
        def loss_fun(params):
            cell = gemmi.UnitCell(*params)
            cell.fractionalization_matrix
            B = np.array(self.cell.fractionalization_matrix).T

    def index_pink(self, max_size=50):
        """
        Infer the setting rotation, self.R using the PinkIndexer algorithm 
        (http://scripts.iucr.org/cgi-bin/paper?S2053273319015559). This implementation
        will use at most max_size reflections
        """
        q = self.s1 - self.s0[None,:]
        q_len = norm2(q, keepdims=True)
        qhat = q / q_len

        # Possible resolution range for each observation
        res_min = self.lam_min / q_len.squeeze(-1)
        res_max = self.lam_max / q_len.squeeze(-1)

        if len(res_min) > max_size:
            dmin = np.sort(res_min)[-max_size]
        else:
            dmin = res_min.min()

        # Generate the feasible set of reflections from the current geometry
        dall = self.cell.calculate_d_array(self.Hall)
        in_res_range = dall >= dmin
        Hall = self.Hall[in_res_range]
        h = (self.B @ Hall.T).T
        hhat = h / norm2(h, keepdims=True)
        dall = dall[in_res_range]

        #Remove entries that incompatible with resolution range of observation
        mask = (res_max[:,None] >= dall) & (res_min[:,None] <= dall)
        i,j = np.where(mask)

        hhat = hhat[j]
        qhat = qhat[i]

        mhat = hhat + qhat
        mhat = mhat / norm2(mhat, keepdims=True)

        phi_grid_size = 200
        phi = np.linspace(-2*np.pi, 0., phi_grid_size)
        c1 = np.sin(phi / 2.)
        c2 = -np.cos(phi / 2.)
        d1 = np.einsum('ad,ad->a', mhat, qhat)
        d2 = np.cross(qhat, mhat)
        theta = 2 * np.arccos(-d1[...,None] * c1)
        num = (mhat[...,None,:] * c2[...,:,None] + d2[...,None,:] * c1[...,:,None])
        denom = np.sin(theta[...,None] / 2.)
        ehat = num / denom
        v = np.arctan(theta / 4)[...,None] * ehat

        #Discretize
        v_grid_size = 100
        rad = np.arctan(np.pi / 4.) #Radius of ball in which rotograms live
        _discretized = np.round(v_grid_size * v / rad).astype('int8')
        discretized = _discretized.reshape((-1, 3))

        # This block is optimized for speed. There are lot of nice, legible ways to 
        # do this in numpy but they can be rather slow. Here's one example alternative:
        # idx,counts = np.unique(discretized, axis=0, return_counts=True)
        # max_idx = np.argmax(counts)
        # grid_max = idx[max_idx] 
        w = 2 * v_grid_size + 1
        to_id = np.arange(w * w * w).reshape((w, w, w))
        flat_ids = to_id[*discretized.T]
        counts = np.bincount(flat_ids)
        max_idx = np.argmax(counts)
        brick = np.ones((w, w, w), dtype='int8')
        vrange = np.arange(-v_grid_size, v_grid_size + 1, dtype='int8')
        vrange = np.roll(vrange, v_grid_size + 1)
        grid_max = np.array([
            (vrange[:,None,None]*brick).flatten()[max_idx],
            (vrange[None,:,None]*brick).flatten()[max_idx],
            (vrange[None,None,:]*brick).flatten()[max_idx],
        ])
        v_max = rad * grid_max / v_grid_size

        #This complicated indexing stuff will figure out the assignment
        #corresponding to the best rotogram voxel.
        flat_assignment = (_discretized == grid_max).all(-1).any(-1)
        Hidx = i[flat_assignment]
        Hobs = Hall[j[flat_assignment]]

        from scipy.linalg import orthogonal_procrustes
        qh = qhat[flat_assignment]
        hh = hhat[flat_assignment]

        a,b,c,alpha,beta,gamma = self.cell.parameters
        def loss_fn(params, return_cell_R_B=False):
            cell = gemmi.UnitCell(a, *params)
            B = np.array(cell.fractionalization_matrix).T
            hh = (B @ Hobs.T).T
            R, _ = orthogonal_procrustes(hh, qh)
            R = R.T
            angles = rs.utils.angle_between(
                (R@hh.T).T,
                qh,
            )
            loss = angles.sum()
            if return_cell_R_B:
                return loss, cell, R, B
            return loss

        from scipy.optimize import minimize
        guess = [b, c, alpha, beta, gamma]
        loss_fn(guess)
        result = minimize(loss_fn, guess)
        #assert result.success
        _, self.cell, self.R, self.B = loss_fn(result.x, True)

class LauePredictor(LaueBase):
    """
    An object to predict spots given a Laue experiment.
    """

    def __init__(self, s0, cell, R, lam_min, lam_max, dmin, spacegroup="1"):
        """
        Calls parent LaueBase init
        """
        super().__init__(s0, cell, R, lam_min, lam_max, dmin, spacegroup=spacegroup)

