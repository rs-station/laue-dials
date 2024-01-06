import numpy as np
import pytest

from laue_dials.algorithms.integrate import Profile, SegmentedImage


@pytest.fixture
def example_profile_data():
    x = np.linspace(1, 1000, num=1000)
    y = np.linspace(1, 1000, num=1000)
    counts = np.array([100.0, 750.0, 60.0])
    return x, y, counts


def test_profile_init(example_profile_data):
    x, y, counts = example_profile_data
    profile = Profile(x, y, counts)
    assert profile is not None
    # Add more assertions based on your expectations


@pytest.fixture
def example_segmented_image_data():
    # Replace this with actual data for creating a SegmentedImage instance
    pixels = np.array([[0, 0, 1], [1, 1, 0], [0, 1, 0]])
    centroids = np.array([[0, 0], [1, 1], [2, 2]])
    return pixels, centroids


def test_segmented_image_init(example_segmented_image_data):
    pixels, centroids = example_segmented_image_data
    segmented_image = SegmentedImage(pixels, centroids)
    assert segmented_image is not None
