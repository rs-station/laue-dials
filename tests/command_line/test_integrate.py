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
    """Captures the kwargs it was constructed with and returns dummy results."""

    calls = []

    def __init__(self, pixels, centroids, radius=None, isigi_cutoff=None, **kwargs):
        RecordingIntegrator.calls.append(
            {"radius": radius, "isigi_cutoff": isigi_cutoff, "n": len(centroids)}
        )
        self.n = len(centroids)

    def fit(self):
        pass

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


def test_isigi_cutoff_is_forwarded_to_integrator(refls_and_imageset, monkeypatch):
    """
    The isigi_cutoff PHIL parameter must reach the Integrator.

    It was previously accepted by integrate_image and then dropped, so the
    class default silently overrode whatever the user asked for.
    """
    refls, img_set = refls_and_imageset
    RecordingIntegrator.calls = []
    monkeypatch.setattr(
        "laue_dials.command_line.integrate.Integrator", RecordingIntegrator
    )

    integrate_image(img_set, refls, isigi_cutoff=4.5, integration_radius=4)

    assert len(RecordingIntegrator.calls) == 1
    assert RecordingIntegrator.calls[0]["isigi_cutoff"] == 4.5


def test_integration_radius_is_forwarded_to_integrator(refls_and_imageset, monkeypatch):
    """An explicit integration_radius is used rather than an estimated one."""
    refls, img_set = refls_and_imageset
    RecordingIntegrator.calls = []
    monkeypatch.setattr(
        "laue_dials.command_line.integrate.Integrator", RecordingIntegrator
    )

    integrate_image(img_set, refls, isigi_cutoff=2.0, integration_radius=7)

    assert RecordingIntegrator.calls[0]["radius"] == 7
