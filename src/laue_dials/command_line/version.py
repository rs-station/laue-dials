#!/usr/bin/env python
"""
This script performs spotfinding on an imported experiment.
"""

from laue_dials.utils.version import laue_version

help_message = """
Prints version information for Laue-DIALS and DIALS.

Examples:

    laue.version
"""


def run():
    """
    Print version information to the terminal and logfile.

    Args:
        None

    Returns:
        None
    """
    # Print version information
    print(laue_version())


if __name__ == "__main__":
    run()
