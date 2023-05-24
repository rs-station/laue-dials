# An example script for processing a fixed-target rotation-series dataset. The data used
# here can be downloaded at https://zenodo.org/record/6407157.

FILE_INPUT_TEMPLATE="PATH/TO/DATA/e080_###.mccd"

# Make needed directories
mkdir stills
mkdir optimized
mkdir refined
mkdir predicted
mkdir integrated

# Import data into DIALS files
dials.import geometry.scan.oscillation=0,1 \
    geometry.goniometer.axes=0,1,0 \
    geometry.beam.wavelength=1.04 \
    geometry.detector.panel.pixel_size=0.08854,0.08854 \
    input.template=$FILE_INPUT_TEMPLATE \

# Get a monochromatic geometry model
laue.initial_solution imported.expt indexer.indexing.known_symmetry.space_group=19 spotfinder.lookup.mask="pixels.mask" indexer.refinement.parameterisation.auto_reduction.action=fix

# Split sequence into stills
laue.sequence_to_stills monochromatic.*
mv split_image* stills/

# Polychromatic portion
N=12 # Max multiprocessing
for i in $(seq -f "%06g" 0 178) # Can change 178 to smaller number of images to analyze
do
  ((j=j%N)); ((j++==0)) && wait; # Batch processing
  (echo "Analyzing image ${i}.";

  laue.optimize_indexing stills/split_image${i}.* output.experiments="optimized/optimized${i}.expt" output.reflections="optimized/optimized${i}.refl" log="optimized/laue.optimize_indexing${i}.log";

  laue.refine optimized/optimized${i}.* output.experiments="refined/poly_refined${i}.expt" output.reflections="refined/poly_refined${i}.refl" output.log="refined/laue.poly_refined${i}.log";

  laue.predict refined/poly_refined${i}.* output.reflections="predicted/predicted${i}.refl" output.log="predicted/laue.predict${i}.log";

  laue.integrate refined/poly_refined${i}.expt predicted/predicted${i}.refl output.filename="integrated/integrated${i}.mtz" output.log="integrated/laue.integrate${i}.log" n_proc=12
  );
done

laue.combine_mtzs integrated/integrated*.mtz

# This is where laue_dials ends. The output file total_integrated.mtz can be merged in careless and refined in phenix to get a model
