All HEWL anomalous data analysis and figure generation post-laue-dials was done using the 5 scripts below, in order:

1_HEWL_anom_cut_friedelize_careless.sh
    - Copies the integrated (unmerged) mtz file produced by the HEWL_anom_laue_dials_processing_final.ipynb notebook into the working directory
    - Calls the cut_unmerged_mtz_by_frames.py utility to create mtzs with only a subset of the overall images
    - Calls the friedelize.py utility to split the Friedel mates into two mtz files (*_plus.mtz and *_minus.mtz)
    - Copies those split mtzs into the appropriate directory
    - Calls the sbatch_careless_varied_frames.sh utility to scale those mtzs

2_HEWL_anom_unfriedelize.sh
    - Calls the unfriedelize.py utility to recombine the Friedel mates into a single mtz file
    - Moves the resulting mtz to the refinement directory

3_HEWL_anom_refine.sh
    - Copies files with a set of custom refinement parameters for each step of refinement in Phenix into the appropriate directory. Refinement 1 is a rigid-body refinment only, while Refinement 2 also refines individual B-factors.
    - Calls the utility sbatch_phenix_Refine.sh to run Phenix refinement

4_HEWL_anom_peak_heights.sh
    - Calls the anomalous_peak_heights.py utility to calculate the anomalous peak heights for each I and S atom accross all frame number sizes and store the resulting outputs in csv files
    - Calls the concatenate_anomalous_peak_csv.py utility to concatenate the resulting 13 csv files into one

5_HEWL_anom_figures.sh
    - Calls the HEWL_anom_peaks.pml utility to generate the PyMOL figure showing anomalous density
    - Calls the careless.ccanom and careless.cchalf function to prepare data for subsequent plotting in Jupyter notebooks
