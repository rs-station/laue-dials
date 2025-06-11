==========================
Jupyter Notebook Tutorials
==========================

These tutorials can help guide you through analyzing a data set with Laue-DIALS. The data for
these tutorials are stored on SBGrid, and can be downloaded using the following in any bash terminal:

.. code:: bash

   # HEWL
   rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/1118 .
   cd 1118 ; shasum -c files.sha

   # PDZ2
   rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/1116 .
   cd 1116 ; shasum -c files.sha

   # DHFR
   rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/1117 .
   cd 1117 ; shasum -c files.sha

These tutorials are tested most recently with version 3.17.0 of DIALS.
To reproduce a python environment for analyzing these data, you can
install DIALS (see `DIALS`_) and use the YAML file specified
:download:`here <manuscript_env.yaml>`.

All referenced scripts can be downloaded :download:`here for PDZ2 <pdz2/pdz2_scripts.tar.gz>` and :download:`here for HEWL <hewl/hewl_scripts.tar.gz>`.


.. _careless: https://github.com/rs-station/careless
.. _DIALS: https://dials.github.io/index.html
.. _issue tracker: https://github.com/rs-station/laue-dials/issues
.. _Phenix: http://www.phenix-online.org
.. _reciprocalspaceship: https://github.com/rs-station/reciprocalspaceship
