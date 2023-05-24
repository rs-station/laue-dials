.. image:: https://github.com/rs-station/laue-dials/actions/workflows/build.yml/badge.svg
   :alt: Build Status
   :target: https://github.com/rs-station/laue-dials/actions/workflows/build.yml

.. image:: https://img.shields.io/pypi/v/laue-dials?color=blue
   :alt: PyPI Release
   :target: https://pypi.org/project/laue-dials/

.. image:: https://codecov.io/gh/rs-station/laue-dials/branch/main/graph/badge.svg
   :alt: Code Coverage
   :target: https://codecov.io/gh/rs-station/laue-dials

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :alt: MIT License
   :target: https://github.com/rs-station/laue-dials/blob/main/LICENSE.txt

==========
laue-dials
==========

Data analysis package for Laue crystallography.

``laue-dials`` is an extension to the `DIALS`_ code for analyzing polychromatic crystallographic data.
Building off the ``DIALS`` framework, and including modern tools like ``numpy``, ``scipy``, and
``reciprocalspaceship``, this package allows for analysis of X-ray crystallographic data using
wide-bandwidth light sources. This package is intended to be used in conjunction with `DIALS`_,
`careless`_, and `Phenix`_ in order to generate molecular models from raw data.

============
Installation
============

To install ``laue-dials``, the DIALS code must first be installed using

.. code:: python

   conda install -c conda-forge dials

Then, using the DIALS python environment, run the following:

.. code:: python

   pip install laue-dials

``laue-dials`` consists of seven command-line scripts for the processing of Laue diffraction data, which are

.. code:: python

   laue.initial_solution
   laue.sequence_to_stills
   laue.optimize_indexing
   laue.refine
   laue.predict
   laue.integrate
   laue.combine_mtzs

Note that you need to import the image data using ``dials.import``. For information on how to use this command, visit https://dials.github.io/documentation/programs/dials_import.html. An example of how to analyze a full dataset lives at https://github.com/rs-station/laue-dials/blob/main/examples/pipeline.sh.

If any issues occur either with installation or use of the software, please file an issue at `issue tracker`_. Any and all questions, concerns, or feature requests are welcome.

.. _careless: https://github.com/rs-station/careless
.. _DIALS: https://dials.github.io/index.html
.. _issue tracker: https://github.com/rs-station/laue-dials/issues
.. _Phenix: http://www.phenix-online.org
.. _reciprocalspaceship: https://github.com/rs-station/reciprocalspaceship
