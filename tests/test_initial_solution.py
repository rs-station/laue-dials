
from laue_dials.command_line.initial_solution import dials_version

__author__ = "PrinceWalnut"
__copyright__ = "PrinceWalnut"
__license__ = "MIT"


def test_dials_version():
    """API Tests"""
    assert dials_version() == True
