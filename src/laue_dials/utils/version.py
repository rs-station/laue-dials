# Get package versioning information
import os
from importlib.metadata import version


def laue_version():
    width = os.get_terminal_size().columns
    dials_version = version("dials")
    laue_dials_version = version("laue-dials")

    print("-" * width)
    print(f"DIALS version " + dials_version + ".")
    print(f"laue-dials version " + laue_dials_version + ".")
    print("-" * width)
