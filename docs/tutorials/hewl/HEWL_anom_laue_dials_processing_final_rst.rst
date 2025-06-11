Installation
============

The easiest way to install ``laue-dials`` and its dependencies is using
`Anaconda <https://docs.anaconda.com/free/anaconda/install/index.html>`__.
First we update and install the libmamba solver with

::

   conda update -n base conda
   conda install -n base conda-libmamba-solver
   conda config --set solver libmamba

With Anaconda, we can then create and activate a custom environment for
your install by running

::

   conda create --name laue-dials
   conda activate laue-dials

Now we are ready to install the main dependency and framework:
`DIALS <https://dials.github.io>`__. After installing that, we can
install ``laue-dials`` using pip, as below:

::

   conda install -c conda-forge dials
   pip install laue-dials

All other dependencies will then be automatically installed for us, and
we’ll be ready to analyze your first Laue data set! Reopen this notebook
with the appropriate environment activated when ready.

Documentation for ``laue-dials`` can be found at
`here <https://rs-station.github.io/laue-dials/index.html>`__, and
entering a command with no arguments on the command line will also print
a help page!

Introduction
============

In this notebook, we will process an anomalous HEWL dataset. The data is
comprised of 3049 frames that constitute several rotations of a HEWL
crystal.

At the end of processing, we would like an integrated ``.mtz`` file we
can then process with ``careless`` for merging. We must run
``laue-dials`` on the whole dataset for completeness, but we can use
fewer images for tutorial purposes if needed. The data are spaced by one
degree per frame, so 180 frames represents a single full rotation.

Data processing will rely on images found in ``./data`` and scripts
found in ``./scripts``.

Importing Data
==============

We can use ``dials.import`` as a way to import the data files written at
experimental facilities into a format that is friendly to both ``DIALS``
and ``laue-dials``. We provide the peak wavelength of the beam, the
detector pixel measurements (in mm), and the goniometer axes.

.. code:: ipython3

    %%time
    %%bash
    
    # Import data
    dials.import geometry.scan.oscillation=0,1 \
        geometry.goniometer.axes=-1,0,0 \
        geometry.beam.wavelength=1.05 \
        geometry.detector.panel.pixel=0.08854,0.08854 \
        input.template=$(pwd)'/data/HEWL_NaI_3_2_####.mccd' \
        output.experiments=imported_HEWL_anom_3049.expt \
        output.log=dials.import_HEWL_anom_3049.log

Getting an Initial Estimate
===========================

After importing our data, the first thing we need to do is get an
initial estimate for the experimental geometry. Here, we’ll use some
monochromatic algorithms from DIALS to help! This step can be tricky –
failure can be due to several causes. In the event of failure, here are
a few common causes:

1. The spotfinding gain is either too high or too low. Try looking at
   the results of ``dials.image_viewer imported.expt strong.refl`` as
   below and seeing if you have too many (or too few) reflections. Lower
   gain gives you more spots, but also more likely to give false
   positives.
2. Supplying the space group or unit cell during indexing can be
   helpful. When supplying the unit cell, allow for some variation in
   the lengths of the axes, since the monochromatic algorithms may
   result in a slightly scaled unit cell depending on the chosen
   wavelength.
3. You may have intensities that need to be masked. These can come from
   bad panels or extraneous scatter. You can use ``dials.image_viewer``
   (described below) to create a mask file for your data, and then
   provide the ``spotfinder.lookup.mask="pixels.mask"`` command below to
   use that mask during spotfinding.

First we will run spotfinding algorithms, then check using the DIALS
image viewer to ensure we have the gain set to an appropriate level, and
then try our hand at ``laue.index`` to see the quality of the indexed
unit cell.

.. code:: ipython3

    %%time
    %%bash
    
    laue.find_spots imported_HEWL_anom_3049.expt \
        spotfinder.mp.nproc=60 \
        spotfinder.threshold.dispersion.gain=0.15 \
        spotfinder.filter.max_separation=10 \
        output.reflections=strong_HEWL_anom_3049.refl \
        output.log=laue.find_spots_HEWL_anom_3049.log

Viewing Images
==============

Sometimes it’s helpful to be able to see the analysis data overlayed on
the raw data. DIALS has a utility for viewing spot information on the
raw images called ``dials.image_viewer``. For example, the spotfinding
gain parameter can be tuned to capture more spots, but lowering it too
much finds nonexistent spots. To check this, we can use the image viewer
to see what spots were found on images. We need to provide an ``expt``
file and a ``refl`` file – the ``imported.expt`` and ``strong.refl``
files will do for checking spotfinding. This program also has utilities
for generating masks if they are needed. The red dots from the checkbox
“Mark centers of mass” are the spots found by ``laue.find_spots`` (which
in turn makes a call to ``dials.find_spots``). These are best used for
judging whether you need to adjust the gain higher (for fewer spots) or
lower (for more) during spotfinding. You can find more details on the
image viewer in the `DIALS tutorial
here <https://dials.github.io/documentation/tutorials/processing_in_detail_betalactamase.html>`__.

.. code:: ipython3

    %%time
    %%bash
    # I found it helpful to set the brightness to 30
    
    dials.image_viewer imported_HEWL_anom_3049.expt strong_HEWL_anom_3049.refl

.. code:: ipython3

    %%time
    %%bash
    
    laue.index imported_HEWL_anom_3049.expt strong_HEWL_anom_3049.refl \
        indexer.indexing.nproc=60 \
        indexer.indexing.known_symmetry.space_group=96 \
        indexer.indexing.refinement_protocol.mode=refine_shells \
        indexer.refinement.parameterisation.auto_reduction.action=fix \
        laue_output.index_only=False \
        laue_output.indexed.experiments=indexed_HEWL_anom_3049.expt \
        laue_output.indexed.reflections=indexed_HEWL_anom_3049.refl \
        laue_output.refined.experiments=refined_HEWL_anom_3049.expt \
        laue_output.refined.reflections=refined_HEWL_anom_3049.refl \
        laue_output.final_output.experiments=monochromatic_HEWL_anom_3049.expt \
        laue_output.final_output.reflections=monochromatic_HEWL_anom_3049.refl \
        laue_output.log=laue.index_HEWL_anom_3049.log

Making Stills
=============

Here we will now split our monochromatic estimate into a series of
stills to prepare it for the polychromatic pipeline. There is a useful
utility called ``laue.sequence_to_stills`` for this.

NOTE: Do not use ``dials.sequence_to_stills``, as there are data columns
which do not match between the two programs.

.. code:: ipython3

    %%time
    %%bash
    
    laue.sequence_to_stills monochromatic_HEWL_anom_3049.expt \
        monochromatic_HEWL_anom_3049.refl \
        output.experiments=stills_HEWL_anom_3049.expt \
        output.reflections=stills_HEWL_anom_3049.refl \
        output.log=laue.sequence_to_stills_HEWL_anom_3049.log

Polychromatic Analysis
======================

Here we will use four other programs in ``laue-dials`` to create a
polychromatic experimental geometry using our initial monochromatic
estimate. Each of the programs does the following:

``laue.optimize_indexing`` assigns wavelengths to reflections and
refines the crystal orientation jointly.

``laue.refine`` is a polychromatic wrapper for ``dials.refine`` and
allows for refining the experimental geometry overall to one suitable
for spot prediction and integration.

``laue.predict`` takes the refined experimental geometry and predicts
the centroids of all strong and weak reflections on the detector.

``laue.integrate`` then builds spot profiles and integrates intensities
on the detector.

.. code:: ipython3

    %%time
    %%bash
    
    laue.optimize_indexing stills_HEWL_anom_3049.refl \
        stills_HEWL_anom_3049.expt \
        output.experiments=optimized_HEWL_anom_3049.expt \
        output.reflections=optimized_HEWL_anom_3049.refl \
        output.log=laue.optimize_indexing_HEWL_anom_3049.log \
        wavelengths.lam_min=0.97 \
        wavelengths.lam_max=1.25 \
        reciprocal_grid.d_min=1.4 \
        nproc=60

.. code:: ipython3

    %%time
    %%bash
    
    laue.refine optimized_HEWL_anom_3049.expt \
        optimized_HEWL_anom_3049.refl \
        output.experiments=poly_refined_HEWL_anom_3049.expt \
        output.reflections=poly_refined_HEWL_anom_3049.refl \
        output.log=laue.poly_refined_HEWL_anom_3049.log \
        nproc=60

Note: even without maxing out the available cores, jupyterlab has a
tendency to crash/think that the above cell is running indefinitely.
After confirming that the experiment & reflections files had been
successfully written out via terminal, I had to interrupt the kernel,
restart it, and then resume processing below.

Check results in image viewer
-----------------------------

.. code:: ipython3

    %%time
    %%bash
    
    dials.image_viewer monochromatic_HEWL_anom_3049.expt monochromatic_HEWL_anom_3049.refl

Predictions do not look great - many shoeboxes do not have a predicted
spot, and there are also some predicted spots that are off-target or
fully false positives.

.. code:: ipython3

    %%time
    %%bash
    
    dials.image_viewer poly_refined_HEWL_anom_3049.expt poly_refined_HEWL_anom_3049.refl

The polychromatic predictions look much better!

Check wavelength spectrum
-------------------------

There is a utility in ``laue-dials`` called ``laue.plot_wavelengths``.
This command generates a histogram of the assigned wavelength spectrum.
If you know approximately the shape of your beam spectrum, this can be a
useful check to ensure that nothing has gone wrong with wavelength
assignment at this stage before predicting the full set of reflections.

.. code:: ipython3

    %%time
    %%bash
    
    laue.plot_wavelengths poly_refined_HEWL_anom_3049.refl \
        refined_only=True \
        save=True \
        show=False \
        output=wavelengths_HEWL_anom_3049.png \
        log=laue.plot_wavelengths_HEWL_anom_3049.log

.. code:: ipython3

    from IPython.display import Image
    import os
    cwd = os.getcwd()
    Image(filename=cwd+'/wavelengths_HEWL_anom_3049.png') 

Spot prediction
---------------

Since the assigned spectrum looks good, we can move on to predicting the
full set of reflections. If the assigned beam spectrum ends up narrower
than the wavelength limits you provided in ``laue.optimize_indexing``,
you can always narrow down the spectrum here for ``laue.predict``. The
predictor will find the locations of all feasible spots and build
profiles for the weak spots based on the observed strong spots. The
output reflection table can then be fed along with the refined ``expt``
file into ``laue.integrate`` to generate ``mtz`` files suitable for
merging in a program like ``careless``.

.. code:: ipython3

    %%time
    %%bash
    
    laue.predict poly_refined_HEWL_anom_3049.expt \
        poly_refined_HEWL_anom_3049.refl \
        output.reflections=predicted_HEWL_anom_3049.refl \
        output.log=laue.predict_HEWL_anom_3049.log \
        wavelengths.lam_min=0.97 \
        wavelengths.lam_max=1.25 \
        reciprocal_grid.d_min=1.4 \
        nproc=60

Integration
-----------

.. code:: ipython3

    %%time
    %%bash
    
    laue.integrate poly_refined_HEWL_anom_3049.expt \
        predicted_HEWL_anom_3049.refl \
        output.filename=integrated_HEWL_anom_3049.mtz \
        output.log=laue.integrate_HEWL_anom_3049.log \
        nproc=12

Conclusion
==========

At this point, you now have integrated ``mtz`` files that you can pass
to `careless <https://github.com/rs-station/careless>`__ for scaling and
merging. We provide an example SLURM-compatible ``careless`` script,
found at ``scripts/sbatch_careless_varied_frames.sh``. There are also
several other scripts that can be used for further processing that are
described by ``README.txt``.

Note that throughout this pipeline, you can use DIALS utilities like
``dials.image_viewer`` or ``dials.report`` to check progress and ensure
your data is being analyzed properly. We recommend regularly checking
the analysis by looking at the data on images, which can be done by

``dials.image_viewer FILE.expt FILE.refl``.

These files are generally written as pairs with the same base name, with
the exception of combining ``imported.expt`` + ``strong.refl``, or
``poly_refined.expt`` + ``predicted.refl``.

Also note that you can take any program and enter it on the command-line
for further help. For example, writing

``laue.optimize_indexing``

will print a help page for the program. You can see all configurable
parameters by using

``laue.optimize_indexing -c``.

This applies to all ``laue-dials`` command-line programs.

For further processing of these data in programs like ``careless``, the
``README.txt`` file includes instructions for using the programs in
``/scripts/`` (reproduced below).

Congratulations! This tutorial is now over. For further questions, feel
free to consult documentation or email the
`authors <https://pypi.org/project/laue-dials/>`__.

Post-Laue-DIALS processing
==========================

All HEWL anomalous data analysis and figure generation post-laue-dials
was done using the 5 scripts below, in order:

1. HEWL_anom_cut_friedelize_careless.sh

   - Copies the integrated (unmerged) mtz file produced by the
     HEWL_anom_laue_dials_processing_final.ipynb notebook into the
     working directory
   - Calls the cut_unmerged_mtz_by_frames.py utility to create mtzs with
     only a subset of the overall images
   - Calls the friedelize.py utility to split the Friedel mates into two
     mtz files (\*_plus.mtz and \*_minus.mtz)
   - Copies those split mtzs into the appropriate directory
   - Calls the sbatch_careless_varied_frames.sh utility to scale those
     mtzs

2. HEWL_anom_unfriedelize.sh

   - Calls the unfriedelize.py utility to recombine the Friedel mates
     into a single mtz file
   - Moves the resulting mtz to the refinement directory

3. HEWL_anom_refine.sh

   - Copies files with a set of custom refinement parameters for each
     step of refinement in Phenix into the appropriate directory.
     Refinement 1 is a rigid-body refinment only, while Refinement 2
     also refines individual B-factors.
   - Calls the utility sbatch_phenix_Refine.sh to run Phenix refinement

4. HEWL_anom_peak_heights.sh

   - Calls the anomalous_peak_heights.py utility to calculate the
     anomalous peak heights for each I and S atom accross all frame
     number sizes and store the resulting outputs in csv files
   - Calls the concatenate_anomalous_peak_csv.py utility to concatenate
     the resulting 13 csv files into one

5. HEWL_anom_figures.sh

   - Calls the HEWL_anom_peaks.pml utility to generate the PyMOL figure
     showing anomalous density
   - Calls the careless.ccanom and careless.cchalf function to prepare
     data for subsequent plotting in Jupyter notebooks
