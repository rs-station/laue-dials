# An example script for processing a fixed-target rotation-series dataset. The data used
# here can be downloaded at https://zenodo.org/record/6407157. You may need to alter
# option values for other datasets.

#SBATCH --account=<project_id>
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00

FILE_INPUT_TEMPLATE="PATH/TO/DATA/e080_###.mccd"

# Import data
# You may need to change some of these values to match your data set
dials.import geometry.scan.oscillation=0,1 \
    geometry.goniometer.axes=0,1,0 \
    geometry.beam.wavelength=1.04 \
    geometry.detector.panel.pixel=0.08854,0.08854 \
    lookup.mask="pixels.mask" \
    input.template=$FILE_INPUT_TEMPLATE

# Get a monochromatic geometry model
laue.find_spots imported.expt \
    spotfinder.mp.nproc=8 \
    spotfinder.threshold.dispersion.gain=0.35 \
    spotfinder.filter.max_separation=10 \
    lookup.mask="pixels.mask"

# Get a monochromatic geometry model
laue.index imported.expt strong.refl \
    indexer.indexing.known_symmetry.space_group=19 \
    indexer.indexing.refinement_protocol.mode=refine_shells \
    indexer.refinement.parameterisation.auto_reduction.action=fix \
    laue_output.index_only=False

# Split sequence into stills
laue.sequence_to_stills monochromatic.*

# Polychromatic pipeline
N=12 # Max multiprocessing
laue.optimize_indexing stills.* \
    output.experiments="optimized.expt" \
    output.reflections="optimized.refl" \
    output.log="laue.optimize_indexing.log" \
    wavelengths.lam_min=0.95 \
    wavelengths.lam_max=1.15 \
    reciprocal_grid.d_min=1.4 \
    n_proc=$N

N=8 # Max multiprocessing
laue.refine optimized.* \
    output.experiments="poly_refined.expt" \
    output.reflections="poly_refined.refl" \
    output.log="laue.poly_refined.log" \
    n_proc=$N

laue.predict poly_refined.* \
    output.reflections="predicted.refl" \
    output.log="laue.predict.log" \
    wavelengths.lam_min=0.95 \
    wavelengths.lam_max=1.15 \
    reciprocal_grid.d_min=1.4 \
    n_proc=$N

laue.integrate poly_refined.expt predicted.refl \
    output.filename="integrated.mtz" \
    output.log="laue.integrate.log" \
    n_proc=$N

# This is where laue_dials ends. The output file integrated.mtz can be merged in careless and refined in phenix to get a model
