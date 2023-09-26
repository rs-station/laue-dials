# An example script for processing a fixed-target rotation-series dataset. The data used
# here can be downloaded at https://zenodo.org/record/6407157.

FILE_INPUT_TEMPLATE="PATH/TO/DATA/e080_###.mccd"

# Import data into DIALS files
dials.import geometry.scan.oscillation=0,1 \
    geometry.goniometer.axes=0,1,0 \
    geometry.beam.wavelength=1.04 \
    geometry.detector.panel.pixel_size=0.08854,0.08854 \
    lookup.mask="pixels.mask" \
    input.template=$FILE_INPUT_TEMPLATE \

# Get a monochromatic geometry model
laue.initial_solution imported.expt indexer.indexing.known_symmetry.space_group=19 spotfinder.lookup.mask="pixels.mask" indexer.refinement.parameterisation.auto_reduction.action=fix

# Split sequence into stills
laue.sequence_to_stills monochromatic.*

# Polychromatic portion
N=12 # Max multiprocessing
laue.optimize_indexing stills.* output.experiments="optimized.expt" output.reflections="optimized.refl" log="laue.optimize_indexing.log" wavelengths.lam_min=0.95 wavelengths.lam_max=1.15 n_proc=$N
laue.refine optimized.* output.experiments="poly_refined.expt" output.reflections="poly_refined.refl" output.log="laue.poly_refined.log" n_proc=$N
laue.predict poly_refined.* output.reflections="predicted.refl" output.log="laue.predict.log" n_proc=$N
laue.integrate poly_refined.expt predicted.refl output.filename="integrated.mtz" output.log="laue.integrate.log" n_proc=$N

# This is where laue_dials ends. The output file total_integrated.mtz can be merged in careless and refined in phenix to get a model
