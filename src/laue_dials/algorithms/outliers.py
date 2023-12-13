"""
Outlier rejection functions.
"""
import numpy as np
import scipy
from dials.array_family import flex


def gen_kde(elist, refls):
    """
    Train a Gaussian Kernel Density Estimator (KDE) on 1/d^2 and wavelengths of submitted strong spots.

    Args:
        elist (dxtbx.model.ExperimentList): The list of experiment objects.
        refls (dials.array_family.flex.reflection_table): The reflection table containing strong spots.

    Returns:
        normalized_resolution (np.ndarray): The normalized resolutions (1/d^2) of the strong spots.
        lams (np.ndarray): The wavelengths of the strong spots.
        kde (scipy.stats.gaussian_kde): The trained Gaussian KDE model.

    Note:
        Harmonic reflections are removed from consideration before training the KDE.
    """
    # Firstly remove all harmonic reflections from consideration
    if "harmonics" in refls:
        harmonics = refls["harmonics"].as_numpy_array()
        refls = refls.select(flex.bool(~harmonics))

    # Get rlps and normalize
    refls.map_centroids_to_reciprocal_space(elist)
    rlps = refls["rlp"].as_numpy_array()
    norms = np.linalg.norm(rlps, axis=1)

    # Get wavelength of rlp
    lams = refls["wavelength"].as_numpy_array()

    # Fit with kernel density estimator
    normalized_resolution = norms**2
    train_data = np.vstack([lams, normalized_resolution])
    kde = scipy.stats.gaussian_kde(train_data)
    return normalized_resolution, lams, kde
