from scipy.spatial import KDTree
import pickle
from matplotlib.colors import colorConverter
import matplotlib as mpl
from matplotlib import pyplot as plt
import pandas as pd
from dials.array_family import flex
import reciprocalspaceship as rs
import numpy as np
from IPython import embed

# Define input files
refl_file = "predicted.refl"
precog_file = "precognition_files/intensities/e080_103.mccd.ii"
image_file = "data/e080_103.mccd"
mask_file = "pixels.mask"

# Various Settings
space_group = 19
match_radius = 5. # Reject dials-precog matches beyond this radius

if refl_file.endswith('predicted.refl'):
    xyzkey = 'xyzcal.px'
elif refl_file.endswith('integrated.refl'):
    xyzkey = 'xyzcal.px'
else:
    xyzkey = 'xyzobs.px.value'


# Do I/O
pixels = plt.imread(image_file)
refls = flex.reflection_table.from_file(refl_file)
df = pd.DataFrame(rs.read_precognition(precog_file).reset_index())
is_indexed = False
if 'miller_index' in refls:
    is_indexed = True

# Match dials spots and put them in df
#x,y = np.array(refls['xyzobs.px.value'])[:,:2].T
x,y = np.array(refls[xyzkey])[:,:2].T
df_dials = pd.DataFrame({
    'X' : x,
    'Y' : y,
})

if is_indexed:
    hkl = np.array(refls['miller_index'])
    df_dials['H'], df_dials['K'], df_dials['L'] = hkl.T

k = KDTree(df[['X', 'Y']])
df_dials['match_distance'],df_dials['precog_idx'] = k.query(df_dials[["X", "Y"]].to_numpy())
df = df.join(df_dials.set_index('precog_idx'), rsuffix='_DIALS')

# plot image with spots overlaid
cmap_name = 'Greys_r'
plotkw = {
    'marker':'o',
    'mfc':'none',
    'ls' : 'none',
    'alpha' : 0.5,
}
plt.matshow(pixels, vmin=1, vmax=50, cmap=cmap_name)
plt.plot(df.X, df.Y, 'b', label='precognition', **plotkw)
df['matched'] = df.match_distance <= match_radius
plt.plot(df[df.matched].X_DIALS, df[df.matched].Y_DIALS, 'y', label='dials (matched)', **plotkw)
plt.plot(df[~df.matched].X_DIALS, df[~df.matched].Y_DIALS, 'r', label='dials (unmatched)', **plotkw)
plt.legend()
plt.title(refl_file)

# overlay the mask
if mask_file is not None:
    mask = pickle.load(open(mask_file, "rb"))[0]
    mask = np.array(mask).reshape(pixels.shape)
    mask_color = 'red'
    transparent = colorConverter.to_rgba('white',alpha=0.0)
    opaque = colorConverter.to_rgba(mask_color, alpha=0.5)
    mask_cmap = mpl.colors.LinearSegmentedColormap.from_list('mask_cmap',[opaque, transparent], 2)
    plt.gca().matshow(mask, cmap=mask_cmap)


if is_indexed:
    matched = df[df.matched]
    H_precog = rs.utils.hkl_to_asu(matched[['H', 'K', 'L']].to_numpy(), space_group)[0]
    H_dials = rs.utils.hkl_to_asu(matched[['H_DIALS', 'K_DIALS', 'L_DIALS']].to_numpy(), space_group)[0]
    matched['assigned'] = (H_dials != 0).any(-1)
    matched['correct'] = (H_precog == H_dials).all(-1)
    matched['label'] = ''
    matched['label'][matched.assigned & matched.correct] = 'correct'
    matched['label'][matched.assigned & ~matched.correct] = 'incorrect'
    matched['label'][~matched.assigned] ='unassigned'
    labels = ['correct', 'incorrect', 'unassigned']
    colors = ['k', 'r', 'b']
    plt.matshow(pixels, vmin=1, vmax=50, cmap=cmap_name)
    for color,label in zip(colors, labels):
        idx = matched['label'] == label
        plt.plot(matched[idx].X_DIALS, matched[idx].Y_DIALS, color, label=label, **plotkw)

    plt.legend()
    plt.title(refl_file)

plt.show()
