# Integration regression fixtures

These fixtures back `tests/regression/test_integration_regression.py`, which
checks that `laue.integrate` still reproduces Precognition intensities on two
DHFR Laue patterns.

## Data source

The diffraction images come from SBGrid dataset **1117** (DHFR Laue,
DOI `10.15785/SBGRID/1117`). The raw images (`e080_015.mccd`, `e080_016.mccd`,
~29 MB each) are **not** version controlled; the test downloads them at run time
via rsync from `rsync://data.sbgrid.org/10.15785/SBGRID/1117/`.

Images 15 and 16 were selected because both have `laue.compute_rmsds` RMSDs
below half a pixel (0.421 px / 0.436 px) and are adjacent. Their experiment ids
are kept as 14 and 15 so the DIALS `BATCH` (= id + 1) equals the image number
(15, 16), matching the Precognition `BATCH`.

## Files

| File | Description |
| --- | --- |
| `poly_refined_2img.expt` | Refined stills geometry for the two images. Image paths are stored as bare basenames and the mask as `pixels.mask`; the test rewrites both to absolute paths at run time. |
| `predicted_2img.refl` | Predicted reflections for the two images, produced by `laue.predict`. |
| `pixels.mask.gz` | Gzipped detector mask (single 3840×3840 panel, ~6% masked). Shipped compressed because the raw pickle is ~14.7 MB; the test expands it before use. |
| `precog_reference_2img.mtz` | Precognition integrated intensities for the same two images, used as the regression reference. |

## Mask provenance

`laue.predict` reads the detector mask via `experiment.imageset.get_mask(0)` and
**discards predictions that land directly on masked pixels** (radius-independent).
`laue.integrate` then reads the same mask, **dilates it by the integration
radius**, and discards predictions whose integration window would overlap the
mask before summing pixels. The predictions in `predicted_2img.refl` were
generated with `pixels.mask`, and the mask is shipped so both the prediction and
integration steps are reproducible.

To use the mask directly (e.g. to re-run `laue.predict`):

```bash
gunzip -k pixels.mask.gz   # -> pixels.mask
```

## Precognition reference

`precog_reference_2img.mtz` was built from the Precognition `.ii` intensity files
for `e080_015` and `e080_016` (fundamental reflections only, i.e. harmonic
multiplicity == 1). Precognition Miller indices are converted to the DIALS
setting with the reindexing operator **(h, k, l) → (−h, k, −l)**, and the batch
number is taken from the image number in the filename (15, 16).

The current baseline is Pearson CC ≈ 0.977 over ~2290 matched reflections.

## Regenerating the fixtures

`regenerate_fixtures.sh` in this directory reproduces every fixture from the full
pipeline outputs. It:

1. Subsets the refined geometry + reflections (`poly_refined.{expt,refl}`) to
   experiment indices 14 and 15 (images 15 and 16).
2. Runs `laue.compute_rmsds` to confirm both images refine below half a pixel.
3. Runs `laue.predict` (which now removes only on-mask centroids) to produce
   `predicted_2img.refl`.
4. Runs `laue.integrate` (which dilates the mask internally) to a scratch MTZ
   used for the cell/spacegroup and baseline correlation.
5. Builds `precog_reference_2img.mtz` from the matching Precognition `.ii` files
   using the reindexing operator above.
6. Gzips the detector mask to `pixels.mask.gz`.

Adjust the `SRC` / `IMAGES_DIR` / `II_DIR` paths at the top of the script for
your environment, run it from a scratch directory, then copy the resulting
fixtures back here (rewriting the `poly_refined_2img.expt` image paths to bare
basenames and the mask to `pixels.mask`).
