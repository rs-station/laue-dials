    .. image:: https://api.cirrus-ci.com/github/rs-station/laue_dials.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/rs-station/laue_dials
    .. image:: https://readthedocs.org/projects/laue_dials/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://laue_dials.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/rs-station/laue_dials/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/rs-station/laue_dials
    .. image:: https://img.shields.io/pypi/v/laue_dials.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/laue_dials/

==========
laue_dials
==========

Data analysis package for Laue crystallography.

``laue_dials`` is an extension to the `DIALS`_ code for analyzing polychromatic crystallographic data.
Building off the ``DIALS`` framework, and including modern tools like ``numpy``, ``scipy``, and
``reciprocalspaceship``, this package allows for analysis of X-ray crystallographic data using
wide-bandwidth light sources. This package is intended to be used in conjunction with `DIALS`_,
`careless`_, and `Phenix`_ in order to generate molecular models from raw data.

============
Installation
============

To install ``laue_dials``, the DIALS code must first be installed, and then using the DIALS python
environment run the following:

.. code:: python

   pip install laue-dials

If any issues occur either with installation or use of the software, please file an issue at `issue tracker`_.

.. _careless: https://github.com/rs-station/careless
.. _DIALS: https://dials.github.io/index.html
.. _issue tracker: https://github.com/rs-station/laue_dials/issues
.. _Phenix: http://www.phenix-online.org
.. _reciprocalspaceship: https://github.com/rs-station/reciprocalspaceship
