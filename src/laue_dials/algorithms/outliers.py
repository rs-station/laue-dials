"""
Outlier rejection functions.
"""
import numpy as np
import scipy
from dials.array_family import flex


def gen_kde(elist, refls):
    """This function trains a Gaussian KDE on 1/d^2 and wavelengths of submitted strong spots"""
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
