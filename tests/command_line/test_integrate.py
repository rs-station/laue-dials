"""
Tests for the laue.integrate command-line layer.

These cover the wiring between the PHIL scope and the ``Integrator``, rather
than the integration maths itself (which lives in
``tests/algorithms/test_integration.py``).
"""

import numpy as np
import pytest
from dials.array_family import flex

from laue_dials.command_line.integrate import integrate_image


class FakeImageSet:
    """Minimal stand-in for a dxtbx imageset with one fully-valid panel."""

    def __init__(self, pixels):
        self._pixels = pixels

    def get_raw_data(self, index):
        return (flex.double(self._pixels.astype(float)),)

    def get_mask(self, index):
        return (flex.bool(np.ones(self._pixels.shape, dtype=bool).flatten()),)


class RecordingIntegrator:
    """Captures the arguments it was called with and returns dummy results."""

    calls = []

    def __init__(
        self, pixels, centroids, radius=None, k=None, isigi_cutoff=None, **kwargs
    ):
        self.record = {
            "radius": radius,
            "k": k,
            "isigi_cutoff": isigi_cutoff,
            "maxiter": None,
            "n": len(centroids),
        }
        RecordingIntegrator.calls.append(self.record)
        self.n = len(centroids)

    def fit(self, maxiter=None):
        self.record["maxiter"] = maxiter

    @property
    def intensity(self):
        return np.full(self.n, 10.0)

    @property
    def uncertainty(self):
        return np.full(self.n, 2.0)

    @property
    def background(self):
        return np.full((self.n, 1), 1.0)


@pytest.fixture
def refls_and_imageset():
    """A small reflection table on a blank image, both on experiment id 0."""
    rng = np.random.default_rng(0)
    pixels = rng.poisson(5.0, size=(60, 60)).astype("float32")

    n = 12
    xy = rng.uniform(15, 45, size=(n, 2))
    refls = flex.reflection_table()
    refls["id"] = flex.int(n, 0)
    refls["xyzcal.px"] = flex.vec3_double(
        flex.double(np.ascontiguousarray(xy[:, 0])),
        flex.double(np.ascontiguousarray(xy[:, 1])),
        flex.double(n, 0.0),
    )
    return refls, FakeImageSet(pixels)


@pytest.mark.parametrize(
    "phil_name, phil_value, recorded_name",
    [
        ("isigi_cutoff", 4.5, "isigi_cutoff"),
        ("integration_radius", 7, "radius"),
        ("knn", 9, "k"),
        ("maxiter", 6, "maxiter"),
    ],
)
def test_phil_parameters_reach_the_integrator(
    refls_and_imageset, monkeypatch, phil_name, phil_value, recorded_name
):
    """
    Every tuning parameter integrate_image accepts must reach the Integrator.

    isigi_cutoff was previously accepted and then dropped at the construction
    site, so the class default silently overrode whatever the user asked for.
    These assertions pin the whole set against that class of regression.
    """
    refls, img_set = refls_and_imageset
    RecordingIntegrator.calls = []
    monkeypatch.setattr(
        "laue_dials.command_line.integrate.Integrator", RecordingIntegrator
    )

    kwargs = {
        "isigi_cutoff": 2.0,
        "integration_radius": 4,
        "knn": 5,
        "maxiter": 2,
        phil_name: phil_value,
    }
    integrate_image(img_set, refls, **kwargs)

    assert len(RecordingIntegrator.calls) == 1
    assert RecordingIntegrator.calls[0][recorded_name] == phil_value
