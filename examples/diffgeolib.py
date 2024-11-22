import numpy as np
import reciprocalspaceship as rs
import gemmi


"""
This is just a place to stash useful functions and classes for diffraction geometry calculations
"""

def indexing_rotation_tril(A):
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
    #R is a "right" or upper triangular matrix
    #Q is an orthogonal (rotation) matrix
    # A = Q.T@R.T
    R,Q = rq(A.T)
    return Q.T, R.T

def normalize(A):
    """
    Normalize the last dimension of an array by dividing by its L2 norm
    """
    return A / np.linalg.norm(A, axis=-1)[...,None]

def hkl2ray(hkl, wavelength=None):
    """ 
    Convert a miller index to the shortest member of its central ray. 
    Optionally, adjust its wavelength accordingly.

    Parameters
    ----------
    hkl : array
        `n x 3` array of miller indices. the dtype must be interpretable as an integeger
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
        return hkl/gcd[...,None], wavelength*gcd
    else:
        return hkl/gcd[...,None]

def is_ray_equivalent(hkl1, hkl2):
    """
    Test for equivalency between two miller indices in a Laue experiment. Returns a boolean array for each of the `n` hkls in `hkl{1,2}`.
    """
    return np.all(np.isclose(hkl2ray(hkl1),  hkl2ray(hkl2)), axis=1)

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
    cc = -1.
    for op in spacegroup.operations():
        aligned_ = rs.utils.apply_to_hkl(target, op)
        cc_mat = np.corrcoef(aligned_.T, reference.T)
        cc_ = np.trace(cc_mat[:3,3:])
        if cc_ > cc:
            aligned = aligned_
            cc = cc_
    if anomalous:
        for op in spacegroup.operations():
            aligned_ = -rs.utils.apply_to_hkl(target, op)
            cc_mat = np.corrcoef(aligned_.T, reference.T)
            cc_ = np.trace(cc_mat[:3,3:])
            if cc_ > cc:
                aligned = aligned_
                cc = cc_
    return aligned


class Detector():
    def __init__(self, ax1, ax2, ori):
        """ 
        This is a simple class to represent a detector panel. 
        """
        super().__init__()
        self.D = np.vstack((
            ax1,
            ax2,
            ori,
        ))

    @classmethod
    def from_inp_file(cls, filename, **kwargs):
        """ doesn't account for tilt yet """
        for line in  open(filename):
            if 'Distance' in line:
                ddist = float(line.split()[1])
            if 'Center' in line:
                beam = np.array([float(line.split()[1]), float(line.split()[1])])
            if 'Pixel' in line:
                size = np.array([float(line.split()[1]), float(line.split()[1])])

        ax1 = np.array([size[0], 0., 0.])
        ax2 = np.array([0., size[1], 0.])
        ori = np.array([
                -beam[0] * size[0], 
                -beam[1] * size[1], 
                ddist
            ])
        return cls(ax1, ax2, ori, **kwargs)

    @classmethod
    def from_expt_file(cls, filename, **kwargs):
        """
        Return a generator of Detector instances for each panel described in `filename`
        """
        from dxtbx.model.experiment_list import ExperimentListFactory
        elist = ExperimentListFactory.from_json_file(filename)
        for detector in elist.detectors():
            for panel in detector.iter_panels():
                ori = np.array(panel.get_origin())
                size = panel.get_pixel_size()
                ax1 = np.array(panel.get_fast_axis())*size[0]
                ax2 = np.array(panel.get_slow_axis())*size[1]
                yield cls(ax1, ax2, ori, **kwargs)

    def pix2lab(self, xy, D=None):
        """ Map xy to 3d coordinates in the lab frame in mm"""
        if D is None:
            D = self.D
        xy = np.pad(xy, ((0,0), (0,1)), constant_values=1.)
        S1 = xy@D
        return S1

    def lab2pix(self, xyz, D=None):
        """ Map 3d coordinates in the lab frame in mm to pixels """
        if D is None:
            D = self.D
        Dinv = np.linalg.inv(D)
        return (xyz@Dinv)[:,:2]

    def s12pix(self, s1, D=None):
        """ Project a scattered beam wavevector into pixel coordinates. """
        if D is None:
            D = self.D
        Dinv = np.linalg.inv(D)
        xya = s1@Dinv
        xy = xya[:,:2]/xya[:,2,None]
        return xy

def orthogonalization(a, b, c, alpha, beta, gamma):
    """ Compute the orthogonalization matrix from cell params """
    alpha,beta,gamma = np.pi*alpha/180.,np.pi*beta/180.,np.pi*gamma/180.
    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    V = a*b*c*np.sqrt(1. - cosa*cosa-cosb*cosb-cosg*cosg + 2*cosa*cosb*cosg)
    return np.array([
        [a, b*cosg, c*cosb],
        [0., b*sing, c*(cosa-cosb*cosg)/sing],
        [0., 0., V/a/b/sing],
    ], device=a.device)

def mat_to_rot_xyz(R, deg=True):
    """ Decompose a rotation matrix into euler angles """
    if R[2, 0] < 1:
        if R[2,0] > -1:
            rot_y = np.asin(-R[2, 0])
            rot_z = np.atan2(R[1, 0], R[0, 0])
            rot_x = np.atan2(R[2, 1], R[2, 2])
        else:
            rot_y = np.pi/2.
            rot_z = -np.atan2(-R[1,2], R[1,1])
            rot_x = 0.
    else:
        rot_y = np.pi/2.
        rot_z = np.atan2(-R[1,2], R[1,1])
        rot_x = 0.
    if deg:
        rot_x = np.rad2deg(rot_x) 
        rot_y = np.rad2deg(rot_y)
        rot_z = np.rad2deg(rot_z)
    return rot_x, rot_y, rot_z

def rot_xyz_to_mat(rot_x, rot_y, rot_z, deg=True):
    """ Convert euler angles into a rotation matrix """
    if deg:
        rot_x = np.deg2rad(rot_x)
        rot_y = np.deg2rad(rot_y)
        rot_z = np.deg2rad(rot_z)
    cx,sx = np.cos(rot_x),np.sin(rot_x)
    cy,sy = np.cos(rot_y),np.sin(rot_y)
    cz,sz = np.cos(rot_z),np.sin(rot_z)

    return np.array([
        [cy*cz, cz*sx*sy-cx*sz, cx*cz*sy+sx*sz],
        [cy*sz, cx*cz+sx*sy*sz, -cz*sx+cx*sy*sz],
        [-sy, cy*sx, cx*cy],
    ], device=rot_x.device)




class LaueAssigner():
    """
    An object to assign miller indices to a laue still. 
    """
    def __init__(self, s0, s1, cell, R, lam_min, lam_max, dmin, spacegroup='1'):
        """
        Parameters
        ----------
        s0 : array
            a 3 vector indicating the direction of the incoming beam wavevector.
            This can be any length, it will be unit normalized in the constructor.
        s1 : array
            n x 3 array indicating the direction of the scatterd beam wavevector.
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
        self.B = np.array(self.cell.fractionalization_matrix).T
        self.ewald_offset = None

        # self.s{0,1} are dynamically masked by their outlier status
        self.s0 = s0 / np.linalg.norm(s0)
        self._s1 = s1 / np.linalg.norm(s1, axis=-1)[:,None]
        self._qobs = self._s1 - self.s0
        self._qpred = np.zeros_like(self._s1)
        self._H = np.zeros_like(self._s1)
        self._wav = np.zeros(len(self._H))
        self._harmonics = np.zeros(len(self._s1), dtype=bool)
        self._inliers = np.ones(len(self._s1), dtype=bool)

        #Initialize the full reciprocal grid
        hmax,kmax,lmax = self.cell.get_hkl_limits(dmin)
        Hall = np.mgrid[
            -hmax:hmax+1:1.,
            -kmax:kmax+1:1.,
            -lmax:lmax+1:1.,
        ].reshape((3, -1)).T
        Hall = Hall[np.any(Hall != 0, axis=1)]
        d = cell.calculate_d_array(Hall)
        Hall = Hall[d >= dmin]

        #Just remove any systematic absences in the space group
        Hall = Hall[~rs.utils.is_absent(Hall, self.spacegroup)]
        self.Hall = Hall

    @property
    def RB(self):
        return self.R@self.B

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

    #<-- setters that operate on the currently inlying set
    def set_qpred(self, qpred):
        self._qpred[self._inliers] = qpred

    def set_H(self, H):
        self._H[self._inliers]  = H
        self._H[~self._inliers] = 0.

    def set_wav(self, wav):
        self._wav[self._inliers]  = wav
        self._wav[~self._inliers] = 0.

    def set_inliers(self, inliers):
        self._inliers[self._inliers] = inliers

    def set_harmonics(self, harmonics):
        self._harmonics[self._inliers] = harmonics
    #--> setters that operate on the currently inlying set

    def reset_inliers(self):
        self._inliers = np.ones(len(self._inliers), dtype=bool)

    def reject_outliers(self, nstd=10.):
        """ update the list of inliers """
        from sklearn.covariance import MinCovDet
        X = np.concatenate((self.qobs, self.qpred * self.wav[:,None]), axis=-1)
        dist = MinCovDet().fit(X).dist_
        self.set_inliers(dist <= nstd**2.)

    def assign(self):
        """ 
        Assign miller indices to the inlier reflections 

        This method will update: 
            self.H     -- miller indices
            self.wav   -- wavelengths
            self.qpred -- predicted scattering vector
        """
        #Generate the feasible set of reflections from the current geometry
        Hall = self.Hall
        qall = (self.RB@Hall.T).T
        feasible = (
            (np.linalg.norm(qall + self.s0/self.lam_min, axis=-1) < 1/self.lam_min) & 
            (np.linalg.norm(qall + self.s0/self.lam_max, axis=-1) > 1/self.lam_max)
        ) 
        Hall = Hall[feasible]
        qall = qall[feasible]

        # Keep track of harmonics in the feasible set TODO
        Raypred = hkl2ray(Hall)
        _,idx, counts = np.unique(Raypred, return_index=True, return_counts=True, axis=0)
        harmonics = (counts > 1)

        #Remove harmonics from the feasible set
        Hall = Hall[idx]
        qall = qall[idx]


        dmat = rs.utils.angle_between(self.qobs[...,None,:], qall[None,...,:])
        cost = dmat

        from scipy.optimize import linear_sum_assignment
        _,idx = linear_sum_assignment(cost)
        H   = Hall[idx]
        qpred = qall[idx]
        harmonics = harmonics[idx]

        # Set all attributes to match the current assignment
        self.set_H(H)
        self.set_qpred(qpred)
        self.set_harmonics(harmonics)

        #wav_pred = -2.*(self.s0 * qpred).sum(-1) / (qpred*qpred).sum(-1)
        wav_obs = np.linalg.norm(self.qobs, axis=-1) / np.linalg.norm(self.qpred, axis=-1)
        self.set_wav(wav_obs)

    def update_rotation(self):
        """ Update the rotation matrix (self.R) based on the inlying refls """
        from scipy.linalg import orthogonal_procrustes
        misset,_ = orthogonal_procrustes(
            self.qobs, 
            self.qpred * self.wav[:,None],
        )
        self.R = misset@self.R

class LauePredictor():
    """
    An object to predict spots given a Laue experiment.
    """
    def __init__(self, s0, cell, R, lam_min, lam_max, dmin, spacegroup='1'):
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
        self.B = np.array(self.cell.fractionalization_matrix).T

        # self.s{0,1} are dynamically masked by their outlier status
        self.s0 = s0 / np.linalg.norm(s0)

        # Initialize the full reciprocal grid
        hmax,kmax,lmax = self.cell.get_hkl_limits(dmin)
        Hall = np.mgrid[
            -hmax:hmax+1:1.,
            -kmax:kmax+1:1.,
            -lmax:lmax+1:1.,
        ].reshape((3, -1)).T
        Hall = Hall[np.any(Hall != 0, axis=1)]
        d = cell.calculate_d_array(Hall)
        Hall = Hall[d >= dmin]
        self.Hall = Hall

    @property
    def RB(self):
        return self.R@self.B

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
        qall = (self.RB@Hall.T).T
        feasible = (
            (np.linalg.norm(qall + self.s0/self.lam_min, axis=-1) < 1/self.lam_min) & 
            (np.linalg.norm(qall + self.s0/self.lam_max, axis=-1) > 1/self.lam_max)
        ) 
        Hall = Hall[feasible]
        qall = qall[feasible]

        # Remove harmonics from the feasible set
        Raypred = hkl2ray(Hall)
        _,idx,counts   = np.unique(Raypred, return_index=True, return_counts=True, axis=0)
        Hall = Hall[idx] # Remove duplicates
        qall = qall[idx]
        if(delete_harmonics):
            idx = (counts > 1) # Remove last harmonic
            Hall = Hall[idx]
            qall = qall[idx]

        # Filter reflections which do not satisfy the resolution-dependent bandwidth
        # TODO: Skip for now and implement this filtration later just to see -- we'll overpredict at high resolution

        # For each q, find the wavelength of the Ewald sphere it lies on
        lams = -2.*(self.s0 * qall).sum(-1) / (qall*qall).sum(-1)

        # Using this wavelength per q, generate s1 vectors
        s0 = self.s0[None,:] / lams[:,None]
        s1_pred = qall + s0

        # Write s1 predictions
        return s1_pred, lams, qall, Hall
