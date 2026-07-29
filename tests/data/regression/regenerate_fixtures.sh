#!/usr/bin/env bash
#
# Regenerate the integration regression fixtures in tests/data/regression/.
#
# Fixtures produced:
#   poly_refined_2img.expt      - refined stills geometry for images 15 & 16
#   predicted_2img.refl         - laue.predict output (on-mask centroids removed)
#   integrated_2img.mtz         - laue.integrate output (scratch; used for CC/cell)
#   precog_reference_2img.mtz   - Precognition reference intensities, batches 15/16
#   pixels.mask.gz              - gzip of the detector mask
#
# Images 15 & 16 (e080_015/016, experiment indices 14/15) were chosen because
# both have laue.compute_rmsds RMSDs below half a pixel (0.421 / 0.436 px) and
# are adjacent. Experiment ids are kept as 14/15 so DIALS BATCH (= id + 1) equals
# the image number (15/16), matching the Precognition BATCH.
#
# Requires a DIALS/laue-dials environment and the local pipeline outputs under
# $SRC. Run from a scratch working directory; copy the resulting fixtures into
# tests/data/regression/ afterwards.
set -euo pipefail

# --- Inputs (adjust to your environment) -----------------------------------
SRC="${SRC:-$HOME/scratch/test_integrator_branch}"          # full pipeline outputs
IMAGES_DIR="${IMAGES_DIR:-$HOME/data/dhfr_data/thirty_images}"  # raw .mccd images
POLY_EXPT="$SRC/poly_refined.expt"     # refined geometry, all 29 stills
POLY_REFL="$SRC/poly_refined.refl"
MASK="$SRC/pixels.mask"                # detector mask (raw pickle)
II_DIR="$SRC/precog/precognition_files/intensities"  # Precognition .ii files

KEEP_IDS="14 15"                       # experiment indices -> images 15 & 16
LAM_MIN=0.95
LAM_MAX=1.15
D_MIN=1.4

# --- 1. Subset the refined geometry + reflections to the two images --------
python - "$POLY_EXPT" "$POLY_REFL" $KEEP_IDS <<'PY'
import sys
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.model import ExperimentList
from dials.array_family import flex

expt_path, refl_path = sys.argv[1], sys.argv[2]
keep_ids = [int(a) for a in sys.argv[3:]]

el = ExperimentListFactory.from_json_file(expt_path, check_format=False)
refls = flex.reflection_table.from_file(refl_path)

sub_el = ExperimentList([el[i] for i in keep_ids])
sel = flex.bool(len(refls), False)
for i in keep_ids:
    sel |= (refls["id"] == i)
sub_refls = refls.select(sel)

idmap = sub_refls.experiment_identifiers()
for k in list(idmap.keys()):
    if k not in keep_ids:
        del idmap[k]

sub_el.as_json("poly_refined_2img.expt")
sub_refls.as_file("poly_refined_2img.refl")
print(f"Subset: {len(sub_el)} experiments, {len(sub_refls)} reflections, ids {sorted(set(sub_refls['id']))}")
PY

# --- 2. Verify the indexing solution (RMSDs should be < 0.5 px) ------------
laue.compute_rmsds poly_refined_2img.expt poly_refined_2img.refl

# --- 3. Predict reflections (new predict: removes on-mask centroids only) ---
laue.predict poly_refined_2img.expt poly_refined_2img.refl \
    output.reflections="predicted_2img.refl" \
    output.log="laue.predict.log" \
    wavelengths.lam_min=$LAM_MIN wavelengths.lam_max=$LAM_MAX \
    reciprocal_grid.d_min=$D_MIN nproc=1

# --- 4. Integrate (dilates the mask by the integration radius internally) ---
#     Produces a scratch MTZ used below for cell/spacegroup and the baseline CC.
laue.integrate poly_refined_2img.expt predicted_2img.refl \
    output.filename="integrated_2img.mtz" \
    output.log="laue.integrate.log" nproc=1

# --- 5. Build the Precognition reference MTZ (batches 15 & 16) --------------
python - "$II_DIR" <<'PY'
import os, re, sys
import numpy as np, pandas as pd, reciprocalspaceship as rs

II_DIR = sys.argv[1]
REINDEX = np.array([-1, 1, -1], dtype=int)   # Precognition -> DIALS: (-h, k, -l)
IMAGES = ["e080_015", "e080_016"]            # BATCH = image number (15, 16)

ref = rs.read_mtz("integrated_2img.mtz")     # cell / spacegroup only

frames = []
for stem in IMAGES:
    data = np.loadtxt(os.path.join(II_DIR, f"{stem}.mccd.ii"))
    if data.ndim == 1:
        data = data[np.newaxis, :]
    data = data[data[:, 3].astype(int) == 1]  # fundamentals only (mult == 1)
    batch = int(re.search(r"e080_(\d+)", stem).group(1))
    frames.append(pd.DataFrame({
        "H": (REINDEX[0] * data[:, 0]).astype(int),
        "K": (REINDEX[1] * data[:, 1]).astype(int),
        "L": (REINDEX[2] * data[:, 2]).astype(int),
        "BATCH": batch,
        "xcal": data[:, 4].astype(np.float32),
        "ycal": data[:, 5].astype(np.float32),
        "wavelength": data[:, 7].astype(np.float32),
        "I": data[:, 8].astype(np.float32),
        "SIGI": data[:, 9].astype(np.float32),
    }))

combined = pd.concat(frames, ignore_index=True)
mtz = rs.DataSet({
    "BATCH": rs.DataSeries(combined["BATCH"].values, dtype="B"),
    "I": rs.DataSeries(combined["I"].values, dtype="J"),
    "SIGI": rs.DataSeries(combined["SIGI"].values, dtype="Q"),
    "xcal": rs.DataSeries(combined["xcal"].values, dtype="R"),
    "ycal": rs.DataSeries(combined["ycal"].values, dtype="R"),
    "wavelength": rs.DataSeries(combined["wavelength"].values, dtype="R"),
}, cell=ref.cell, spacegroup=ref.spacegroup)
mtz.index = pd.MultiIndex.from_arrays(
    [combined["H"].values, combined["K"].values, combined["L"].values],
    names=["H", "K", "L"])
mtz.write_mtz("precog_reference_2img.mtz")
print(f"Wrote precog_reference_2img.mtz: {len(mtz)} reflections, batches {sorted(set(combined['BATCH']))}")
PY

# --- 6. Compress the detector mask -----------------------------------------
gzip -c "$MASK" > pixels.mask.gz

# --- 7. Report the baseline correlation vs Precognition --------------------
python - <<'PY'
import numpy as np, reciprocalspaceship as rs
dials = rs.read_mtz("integrated_2img.mtz").reset_index()
precog = rs.read_mtz("precog_reference_2img.mtz").reset_index()
m = dials.merge(precog, on=["H", "K", "L", "BATCH"], suffixes=("_d", "_p"))
i_d, i_p = m["I_d"].to_numpy(float), m["I_p"].to_numpy(float)
g = np.isfinite(i_d) & np.isfinite(i_p)
print(f"Baseline: n_matched={g.sum()}, Pearson CC={np.corrcoef(i_d[g], i_p[g])[0,1]:.4f}")
PY

echo
echo "Fixtures to copy into tests/data/regression/:"
echo "  poly_refined_2img.expt  (rewrite image paths to bare basenames + mask='pixels.mask')"
echo "  predicted_2img.refl"
echo "  precog_reference_2img.mtz"
echo "  pixels.mask.gz"
