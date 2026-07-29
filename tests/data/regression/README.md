# Integration regression fixtures

These fixtures back `tests/regression/test_integration_regression.py`, which
checks that `laue.integrate` still reproduces Precognition intensities on two
DHFR Laue patterns.

## Data source

The diffraction images come from SBGrid dataset **1117** (DHFR Laue,
DOI `10.15785/SBGRID/1117`). The raw images (`e080_001.mccd`, `e080_002.mccd`,
~29 MB each) are **not** version controlled; the test downloads them at run time
via rsync from `rsync://data.sbgrid.org/10.15785/SBGRID/1117/`.

## Files

| File | Description |
| --- | --- |
| `poly_refined_2img.expt` | Refined stills geometry for the two images. Image paths are stored as bare basenames and the mask as `pixels.mask`; the test rewrites both to absolute paths at run time. |
| `predicted_2img.refl` | Predicted reflections for the two images, produced by `laue.predict`. |
| `pixels.mask.gz` | Gzipped detector mask (single 3840×3840 panel, ~6% masked). Shipped compressed because the raw pickle is ~14.7 MB; the test expands it before use. |
| `precog_reference_2img.mtz` | Precognition integrated intensities for the same two images, used as the regression reference. |

## Mask provenance

`laue.predict` reads the detector mask via `experiment.imageset.get_mask(0)` and
**discards predictions that land in masked (or mask-dilated) regions**. The
predictions in `predicted_2img.refl` were therefore generated with
`pixels.mask`, so the mask is included here to keep the prediction step
reproducible. `laue.integrate` does **not** read the mask — it sums raw pixel
windows around each predicted centroid — so the mask does not affect the
correlation the regression test computes. (Verified: integrating with and
without the mask referenced yields an identical Pearson CC of 0.978 over the
same 2225 matched reflections.)

To use the mask directly (e.g. to re-run `laue.predict`):

```bash
gunzip -k pixels.mask.gz   # -> pixels.mask
```

## Precognition reference

`precog_reference_2img.mtz` was built from the Precognition `.ii` intensity files
for `e080_001` and `e080_002` (fundamental reflections only, i.e. harmonic
multiplicity == 1). Precognition Miller indices are converted to the DIALS
setting with the reindexing operator **(h, k, l) → (−h, k, −l)**, and the batch
number is taken from the image number in the filename.

## Regenerating the fixtures

1. Run the laue-dials pipeline (import → find_spots → index →
   sequence_to_stills → optimize_indexing → refine → predict) on the DHFR
   images, using `pixels.mask` as the lookup mask.
2. Subset the refined experiments and predictions to the first two images
   (experiment ids 0 and 1).
3. Build `precog_reference_2img.mtz` from the matching Precognition `.ii` files
   using the reindexing operator above.
