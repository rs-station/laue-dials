"""
    Setup file for laue_dials.
    Use setup.cfg to configure.
"""
from setuptools import setup

if __name__ == "__main__":
    try:
        setup(use_scm_version=True, setup_requires=["setuptools_scm"])
    except:  # noqa
        print(
            "\n\nAn error occurred while building the project, "
            "please ensure you have the most updated version of setuptools, "
            "setuptools_scm and wheel with:\n"
            "   pip install -U setuptools setuptools_scm wheel\n\n"
        )
        raise
