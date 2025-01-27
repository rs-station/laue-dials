# The values in this script are for the DHFR data at 
# https://zenodo.org/records/10199220. You may need 
# to change these values if analyzing a different data set.

FILE_INPUT_TEMPLATE="/PATH/TO/DATA"
MASK="/PATH/TO/MASK/FILE"
N=8

# Import data
# You may need to change some of these values to match your data set
dials.import geometry.scan.oscillation=0,1 \
    geometry.goniometer.axes=0,1,0 \
    convert_sequences_to_stills=True \
    geometry.beam.wavelength=1.04 \
    geometry.detector.panel.pixel=0.08854,0.08854 \
    lookup.mask=$MASK \
    input.template=$FILE_INPUT_TEMPLATE

# Get a monochromatic geometry model
dials.find_spots imported.expt \
    output.shoeboxes=False \
    spotfinder.mp.nproc=$N \
    spotfinder.threshold.dispersion.gain=0.15 \
    spotfinder.filter.max_separation=10 \
    lookup.mask=$MASK 

# Get a monochromatic geometry model
laue.index imported.expt strong.refl \
    indexer.indexing.nproc=$N \
    indexer.indexing.method="pink_indexer" \
    indexer.indexing.pink_indexer.wavelength=1.1 \
    indexer.indexing.pink_indexer.percent_bandwidth=15 \
    indexer.indexing.pink_indexer.max_refls=50 \
    indexer.indexing.pink_indexer.min_lattices=5 \
    indexer.indexing.pink_indexer.rotogram_grid_points=360 \
    indexer.indexing.pink_indexer.voxel_grid_points=250 \
    indexer.indexing.known_symmetry.space_group=19 \
    indexer.indexing.known_symmetry.unit_cell=34.297,45.552,99.035,90,90,90 \
    indexer.indexing.refinement_protocol.mode=None \
    indexer.indexing.stills.ewald_proximal_volume_max=0.03 \
    indexer.indexing.stills.rmsd_min_px=20 \
    indexer.indexing.stills.refine_all_candidates=True \
    indexer.indexing.joint_indexing=False \
    indexer.refinement.reflections.outlier.algorithm=None \
    indexer.refinement.parameterisation.auto_reduction.action=fix \
    indexer.refinement.parameterisation.scan_varying=False \
    laue_output.index_only=True

# Polychromatic pipeline
laue.optimize_indexing monochromatic.expt monochromatic.refl \
    output.experiments="optimized.expt" \
    output.reflections="optimized.refl" \
    output.log="laue.optimize_indexing.log" \
    wavelengths.lam_min=1.00 \
    wavelengths.lam_max=1.20 \
    reciprocal_grid.d_min=2.0 \
    geometry.unit_cell=34.297,45.552,99.035,90,90,90 \
    n_macrocycles=5 \
    keep_unindexed=False \
    filter_spectrum=True \
    nproc=$N

laue.refine optimized.* \
    output.experiments="poly_refined.expt" \
    output.reflections="poly_refined.refl" \
    output.log="laue.poly_refined.log" \
    nproc=$N

laue.predict poly_refined.* \
    output.reflections="predicted.refl" \
    output.log="laue.predict.log" \
    wavelengths.lam_min=1.00 \
    wavelengths.lam_max=1.20 \
    reciprocal_grid.d_min=1.4 \
    nproc=$N

laue.integrate poly_refined.expt predicted.refl \
    output.filename="integrated.mtz" \
    output.log="laue.integrate.log" \
    nproc=$N
