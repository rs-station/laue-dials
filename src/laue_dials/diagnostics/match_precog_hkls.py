from glob import glob

import gemmi
import numpy as np
import pandas as pd
import reciprocalspaceship as rs
from diffgeolib import *
from matplotlib import pyplot as plt
from tqdm import tqdm, trange

plt.style.use("tableau-colorblind10")


def parse_ii_inp_file_pairs(ii_filenames, inp_filenames, spacegroup=None, log=None):
    data = None
    sg_number = None
    cell = np.zeros(6)
    for file_number, (iiFN, inpFN) in enumerate(
        tqdm(list(zip(ii_filenames, inp_filenames)))
    ):
        df = rs.read_precognition(iiFN)
        with open(inpFN) as f:
            for line in f:
                if "Crystal" in line:
                    cell += np.array(list(map(float, line.split()[1:7]))) / len(
                        ii_filenames
                    )
                    if spacegroup is None:
                        sg_number = int(line.split()[7])
                if "Pixel" in line:
                    float(line.split()[1])
                    float(line.split()[2])
        #        del(df['Resolution'])

        # If log is given, use it to determine the image number
        if log is None:
            df["BATCH"] = file_number
        else:
            entry = log.loc[log.file == inpFN[:-4]]
            assert len(entry) == 1
            df["BATCH"] = entry.index.values[0]

        # Purge multiples from Precognition processing
        # These will be recomputed during merging later
        df = df.reset_index().groupby(["X", "Y"], as_index=False).first()
        df.cell = rs.dataset.gemmi.UnitCell(*cell)
        if sg_number is not None:
            df.spacegroup = rs.dataset.gemmi.find_spacegroup_by_number(sg_number)
        if file_number == 0:
            data = df
        else:
            data = rs.concat([data, df], check_isomorphous=False)

    del data["Multiplicity"]
    data["Precognition"] = [1] * len(data)
    data["Indexed"] = [True] * len(data)
    data["Intensity"] = [0] * len(data)
    data.set_index(["H", "K", "L"], inplace=True)
    data.infer_mtz_dtypes(inplace=True)
    return data


# expt_file = "aligned.expt"
# expt_file = "processing/monochromatic.expt"
# refl_file = "processing/monochromatic.refl"
expt_file = "processing/optimized.expt"
refl_file = "processing/optimized.refl"

# centroid distance cutoff in pixels
centroid_max_distance = 10.0

# ds = rs.read_precognition(ii_file).reset_index()
print("parsing precog files")
precog_df = parse_ii_inp_file_pairs(
    sorted(glob("precognition_files/intensities/e080_???.mccd.ii")),
    sorted(glob("precognition_files/input/e080_???.mccd.inp")),
).reset_index()

print("reading DIALS files")
from dials.array_family.flex import reflection_table
from dxtbx.model.experiment_list import ExperimentListFactory

elist = ExperimentListFactory.from_json_file(expt_file, False)
cryst = elist.crystals()[0]
unit_cell = cryst.get_unit_cell()
precog_df.spacegroup = gemmi.SpaceGroup(
    cryst.get_space_group().type().universal_hermann_mauguin_symbol()
)
precog_df.cell = gemmi.UnitCell(*unit_cell.parameters())
refls = reflection_table.from_file(refl_file)
# refls = refls.select(refls.get_flags(refls.flags.indexed))
# refls['wavelength'] = flex.double(np.zeros(len(refls)))
# refls = refls.select(refls.get_flags(refls.flags.used_in_refinement))

print("generating DIALS dataframe")
dials_df = rs.DataSet(
    {
        "X": refls["xyzobs.px.value"].parts()[0].as_numpy_array(),
        "Y": refls["xyzobs.px.value"].parts()[1].as_numpy_array(),
        "H": refls["miller_index"].as_vec3_double().parts()[0].as_numpy_array(),
        "K": refls["miller_index"].as_vec3_double().parts()[1].as_numpy_array(),
        "L": refls["miller_index"].as_vec3_double().parts()[2].as_numpy_array(),
        "Wavelength": refls["wavelength"].as_numpy_array(),
        "Resolution": np.linalg.norm(refls["rlp"].as_numpy_array(), axis=1) ** 2,
        "BATCH": refls["image_id"].as_numpy_array(),
        "Precognition": [0] * len(refls),
        "Indexed": [1]
        * len(refls),  # np.asarray(refls.get_flags(refls.flags.indexed), dtype=bool),
        "Intensity": refls["intensity.sum.value"].as_numpy_array(),
    },
    cell=precog_df.cell,
    spacegroup=precog_df.spacegroup,
).infer_mtz_dtypes()
dials_df["Indexed"] = dials_df["Indexed"].astype(bool)

print("initializing metrics")
percent_correct = np.zeros(len(elist))
percent_outliers = np.zeros(len(elist))
percent_misindexed = np.zeros(len(elist))
nspots = np.zeros(len(elist))
nmatch = np.zeros(len(elist))
h_diff = []
k_diff = []
l_diff = []
manhattan_hkl = []

both_df = None

# Iterate by frame and match HKLs, seeing what percentage are correct
for i in trange(len(elist)):
    # Get reflection indices from each batch
    im_pre = precog_df[precog_df["BATCH"] == i]
    im_dia = dials_df[dials_df["BATCH"] == i]
    nspots[i] = len(im_dia)

    # Removing the unindexed refls
    outliers = np.all(im_dia[["H", "K", "L"]] == 0.0, axis=1)
    percent_outliers[i] = 100.0 * outliers.sum() / len(outliers)
    xy_outliers = im_dia.loc[outliers, ["X", "Y"]].to_numpy(float)
    # im_dia = im_dia[~outliers]
    if len(im_dia) == 0:
        continue

    dmat = np.linalg.norm(
        im_dia[["X", "Y"]].to_numpy(float)[:, None, :]
        - im_pre[["X", "Y"]].to_numpy(float)[None, :, :],
        axis=-1,
    )

    # This prevents duplicated matches
    idx1, idx2 = np.where(
        (dmat == dmat.min(0))
        & (dmat == dmat.min(1)[:, None])
        & (dmat <= centroid_max_distance)
    )
    im_dia = im_dia.iloc[idx1]
    im_pre = im_pre.iloc[idx2]

    if len(im_dia) == 0:
        continue

    precog_hkl = im_pre[["H", "K", "L"]].to_numpy(float)
    dials_hkl = im_dia[["H", "K", "L"]].to_numpy(float)

    # Get XY positions for refls
    precog_xy = im_pre[["X", "Y"]].to_numpy(float)
    dials_xy = im_dia[["X", "Y"]].to_numpy(float)

    # Align precog to DIALS hkls
    aligned_hkls = align_hkls(
        dials_hkl, precog_hkl, precog_df.spacegroup
    )  # rs.utils.hkl_to_asu for mapping to same asu, need hkls + spacegroup
    #    dials_hkl = rs.utils.hkl_to_asu(dials_hkl, 19)[0]
    #    precog_hkl = rs.utils.hkl_to_asu(precog_hkl, 19)[0]
    #    aligned_hkls = precog_hkl

    # Check correctness of matching
    correct = is_ray_equivalent(aligned_hkls, dials_hkl)
    if len(correct) > 0:
        percent_misindexed[i] = (
            100.0 * sum(~correct[im_dia["Indexed"]]) / len(correct[im_dia["Indexed"]])
        )
        percent_correct[i] = (
            100.0 * sum(correct[im_dia["Indexed"]]) / len(correct[im_dia["Indexed"]])
        )

    nmatch[i] = len(im_dia)

    # Add this image to `both_df`
    im_pre.loc[:, ["H", "K", "L"]] = aligned_hkls
    im_pre.infer_mtz_dtypes(inplace=True)
    _both_df = im_dia.reset_index().join(
        im_pre.reset_index(), rsuffix="_pre", lsuffix="_dia"
    )
    # _both_df = im_dia.join(im_pre.set_index(['H', 'K', 'L']), on=['H', 'K', 'L'], lsuffix='_dia', rsuffix='_pre')
    _both_df["correct"] = np.all([correct, _both_df["Indexed_dia"]], axis=0)
    _both_df["incorrect"] = np.all([~correct, _both_df["Indexed_dia"]], axis=0)
    both_df = pd.concat((both_df, _both_df))

    # Get Manhattan distance of HKLs
    h_diff.append(both_df["H_dia"] - both_df["H_pre"])
    k_diff.append(both_df["K_dia"] - both_df["K_pre"])
    l_diff.append(both_df["L_dia"] - both_df["L_pre"])
    manhattan_hkl.append(np.abs(h_diff[-1]) + np.abs(k_diff[-1]) + np.abs(l_diff[-1]))

#    if i == 0:
#        x_diff = precog_xy[:,0] - dials_xy[:,0]
#        y_diff = precog_xy[:,1] - dials_xy[:,1]
#
#        plt.figure()
#        plt.scatter(dials_xy[:,0][correct], dials_xy[:,1][correct], c='b', alpha=1, label='Correct')
#        plt.scatter(dials_xy[:,0][incorrect], dials_xy[:,1][incorrect], c='r', alpha=1, label='Incorrect')
#        plt.xlabel('X (pixels)')
#        plt.ylabel('Y (pixels)')
#        plt.title('DIALS vs Precog Spot Centroids (Single Image)')
#        plt.legend()
#        plt.show()
#
#        plt.figure()
#        plt.plot(
#            _both_df.loc[_both_df.correct, 'Wavelength_pre'].to_numpy(),
#            _both_df.loc[_both_df.correct, 'Wavelength_dia'].to_numpy(),
#            'k.',
#            alpha=0.5,
#        )
#        plt.xlabel('$\lambda$ (Precognition)')
#        plt.ylabel('$\lambda$ (DIALS)')
#        plt.title('DIALS vs Precog Wavelengths (Single Image)')
#        plt.show()

both_df.to_csv("matched_hkls.csv", index=False)

plt.figure()
plt.plot(np.arange(len(percent_correct)), percent_correct)
plt.xlabel("Image Number")
plt.ylabel("Percent")
plt.title("Fraction Reflections Correctly Indexed by Image")
plt.savefig("correctly_indexed_frac_per_image.png")
plt.show()

frames_misindexed = pd.DataFrame(
    {"Frame": np.arange(len(percent_misindexed)), "% Misindexed": percent_misindexed}
)
plt.figure()
plt.plot(np.arange(len(percent_correct)), percent_misindexed)
plt.xlabel("Image Number")
plt.ylabel("Percent")
plt.title("Fraction Reflections Misindexed by Image")
plt.savefig("incorrectly_indexed_frac_per_image.png")
plt.show()

plt.figure()
plt.plot(np.arange(len(nspots)), nspots, label="Strong Spots")
plt.plot(np.arange(len(nmatch)), nmatch, label="Matched to Precog")
plt.xlabel("Image Number")
plt.ylabel("Count")
plt.title("Spots per Image")
plt.legend()
plt.savefig("matched_strong_spots_per_image.png")
plt.show()

plt.figure()
plt.plot(
    both_df.loc[both_df.correct, "Wavelength_pre"].to_numpy(),
    both_df.loc[both_df.correct, "Wavelength_dia"].to_numpy(),
    "k.",
    alpha=0.1,
)
plt.xlabel("$\lambda$ (Precognition)")
plt.ylabel("$\lambda$ (DIALS)")
plt.savefig("wavelength_correctly_indexed.png")
plt.show()

labels = ["Correct HKL", "Incorrect HKL", "Unindexed"]
plt.figure()
plt.scatter(
    both_df.loc[both_df.correct, "Resolution_dia"].to_numpy(),
    both_df.loc[both_df.correct, "Intensity_dia"].to_numpy(),
    alpha=0.2,
    label=labels[0],
)
plt.scatter(
    both_df.loc[both_df.incorrect, "Resolution_dia"].to_numpy(),
    both_df.loc[both_df.incorrect, "Intensity_dia"].to_numpy(),
    alpha=0.2,
    label=labels[1],
)
plt.scatter(
    both_df.loc[~both_df.Indexed_dia, "Resolution_dia"].to_numpy(),
    both_df.loc[~both_df.Indexed_dia, "Intensity_dia"].to_numpy(),
    alpha=0.2,
    label=labels[2],
)
ax = plt.gca()
ax.set_yscale("log")
plt.xlabel("DIALS rlp Resolution (1/d^2)")
plt.ylabel("DIALS Intensity")
plt.savefig("res_vs_intensity.png")
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())
plt.show()

plt.figure()
plt.scatter(
    both_df.loc[both_df.correct, "Resolution_dia"].to_numpy(),
    both_df.loc[both_df.correct, "Wavelength_pre"].to_numpy(),
    alpha=0.2,
)
plt.xlabel("DIALS rlp Resolution (1/d^2)")
plt.ylabel("Precognition Wavelength (A)")
plt.title("Correctly Indexed")
plt.savefig("res_vs_lambda_correct.png")
plt.show()

plt.figure()
plt.scatter(
    both_df.loc[both_df.incorrect, "Resolution_dia"].to_numpy(),
    both_df.loc[both_df.incorrect, "Wavelength_pre"].to_numpy(),
    alpha=0.2,
)
plt.xlabel("DIALS rlp Resolution (1/d^2)")
plt.ylabel("Precognition Wavelength (A)")
plt.title("Incorrectly Indexed")
plt.savefig("res_vs_lambda_incorrect.png")
plt.show()

plt.figure()
plt.scatter(
    both_df.loc[~both_df.Indexed_dia, "Resolution_dia"].to_numpy(),
    both_df.loc[~both_df.Indexed_dia, "Wavelength_pre"].to_numpy(),
    alpha=0.2,
)
plt.xlabel("DIALS rlp Resolution (1/d^2)")
plt.ylabel("Precognition Wavelength (A)")
plt.title("Unindexed")
plt.savefig("res_vs_lambda_unindexed.png")
plt.show()

plt.figure()
c1, c2, c3 = "#1b9e77", "#d95f02", "#7570b3"
alpha = 0.1
cor = both_df.correct
plt.plot(
    both_df.loc[cor, "H_pre"].to_numpy(),
    both_df.loc[cor, "H_dia"].to_numpy(),
    ".",
    color=c1,
    label="H (correct)",
    alpha=alpha,
)
plt.plot(
    both_df.loc[cor, "K_pre"].to_numpy(),
    both_df.loc[cor, "K_dia"].to_numpy(),
    ".",
    color=c2,
    label="K (correct)",
    alpha=alpha,
)
plt.plot(
    both_df.loc[cor, "L_pre"].to_numpy(),
    both_df.loc[cor, "L_dia"].to_numpy(),
    ".",
    color=c3,
    label="L (correct)",
    alpha=alpha,
)
plt.xlabel("Precognition")
plt.ylabel("DIALS")
plt.legend()
plt.savefig("correct_hkl_per_image.png")
plt.show()

lam_diffs = both_df.loc[cor, "Wavelength_dia"] - both_df.loc[cor, "Wavelength_pre"]
plt.figure()
plt.hist(lam_diffs, bins=100)
plt.xlabel("Wavelength Error (Angstroms)")
plt.ylabel("Number of Spots")
plt.savefig("wavelength_error.png")
plt.show()

from itertools import chain

manhattan_hkl_flattened = list(chain(*manhattan_hkl))
manhattan_hkl_flattened = [
    i for i in manhattan_hkl_flattened if i != 0
]  # Only misindexed spots
plt.figure()
plt.hist(manhattan_hkl_flattened, bins=len(np.unique(manhattan_hkl_flattened)))
plt.xlabel("HKL Manhattan Distance")
plt.ylabel("Number of Spots")
plt.savefig("manhattan_hkl.png")
plt.show()

manhattan_distances = np.zeros(len(manhattan_hkl))
for i in range(len(manhattan_distances)):
    misindexed = []
    misindexed = [
        j for j in manhattan_hkl[i] if j != 0
    ]  # Only misindexed spots per image
    manhattan_distances[i] = np.mean(misindexed)
plt.figure()
plt.plot(np.arange(len(manhattan_distances)), manhattan_distances)
plt.xlabel("Image Number")
plt.ylabel("Average HKL Manhattan Distance")
plt.title("Incorrectly Indexed Reflections")
plt.savefig("manhattan_hkl_avg_per_image.png")
plt.show()

manhattan_distances = np.zeros(len(manhattan_hkl))
for i in range(len(manhattan_distances)):
    all_refls = []
    all_refls = [j for j in manhattan_hkl[i]]
    manhattan_distances[i] = np.mean(all_refls)
plt.figure()
plt.plot(np.arange(len(manhattan_distances)), manhattan_distances)
plt.xlabel("Image Number")
plt.ylabel("Average HKL Manhattan Distance")
plt.title("All Indexed Reflections")
plt.savefig("manhattan_hkl_avg_per_image_all_refls.png")
plt.show()
