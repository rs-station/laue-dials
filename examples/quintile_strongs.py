from diffgeolib import *
from tqdm import trange,tqdm
from matplotlib import pyplot as plt
from scitbx import matrix
from IPython import embed
from glob import glob
from dials.array_family import flex
import numpy as np
import reciprocalspaceship as rs
import gemmi
import pandas as pd

plt.style.use('tableau-colorblind10')

def parse_ii_inp_file_pairs(ii_filenames, inp_filenames, spacegroup=None, log=None):
    data = None
    sg_number = None
    p = []
    cell = np.zeros(6)
    for file_number, (iiFN,inpFN) in enumerate(tqdm(list(zip(ii_filenames, inp_filenames)))):
        df = rs.read_precognition(iiFN)
        with open(inpFN) as f:
            for line in f:
                if "Crystal" in line:
                    cell +=  np.array(list(map(float, line.split()[1:7])))/len(ii_filenames)
                    if spacegroup is None:
                        sg_number  = int(line.split()[7])
                if "Pixel" in line:
                    pixel_x = float(line.split()[1])
                    pixel_y = float(line.split()[2])

        # If log is given, use it to determine the image number
        if log is None:
            df['BATCH'] = file_number
        else:
            entry = log.loc[log.file == inpFN[:-4]]
            assert len(entry) == 1
            df['BATCH'] = entry.index.values[0]

        #Purge multiples from Precognition processing
        #These will be recomputed during merging later
        df = df.reset_index().groupby(['X', 'Y'], as_index=False).first()
        if file_number == 0:
            data = df
        else:
            data = rs.concat([data, df], check_isomorphous=False)

    del(data['H'])
    del(data['K'])
    del(data['L'])
    del(data['Wavelength'])
    del(data['I'])
    del(data['SigI'])
    data.infer_mtz_dtypes(inplace=True)
    return data

print('reading DIALS files')
#refl_file = "strong.refl"
refl_file = "monochromatic.refl"
from dials.array_family.flex import reflection_table
refls = reflection_table.from_file(refl_file)
refls = refls.select(refls.get_flags(refls.flags.indexed))

#centroid distance cutoff in pixels
centroid_max_distance = 10.

#ds = rs.read_precognition(ii_file).reset_index()
print('parsing precog files')
precog_df = parse_ii_inp_file_pairs(
    sorted(glob('precognition_files/intensities/e080_???.mccd.ii')),
    sorted(glob('precognition_files/input/e080_???.mccd.inp')),
).reset_index()

print('generating DIALS dataframe')
dials_df = rs.DataSet({
    'X' : refls['xyzobs.px.value'].parts()[0].as_numpy_array(),
    'Y' : refls['xyzobs.px.value'].parts()[1].as_numpy_array(),
    'BATCH' : refls['image_id'].as_numpy_array(),
}).infer_mtz_dtypes()
    #'BATCH' : refls['id'].as_numpy_array(),

print('initializing metrics')
quintiles_by_frame = []
frames = np.unique(dials_df['BATCH'])
matched_percent = np.zeros(len(frames))

# Iterate by frame and match HKLs, seeing what percentage are correct
for i in trange(len(frames)):
    # Get reflection indices from each batch
    im_pre = precog_df[precog_df['BATCH'] == i]
    im_dia = dials_df[dials_df['BATCH'] == i]

    dmat = np.linalg.norm(
        im_dia[['X', 'Y']].to_numpy(float)[:,None,:] - \
        im_pre[['X', 'Y']].to_numpy(float)[None,:,:],
        axis = -1
    )

    # This prevents duplicated matches
    total_refls = len(im_dia)
    idx1,idx2 = np.where((dmat == dmat.min(0)) & (dmat == dmat.min(1)[:,None]) & (dmat <= centroid_max_distance))
    im_dia = im_dia.iloc[idx1]
    im_pre = im_pre.iloc[idx2]
    matched_refls = len(im_dia)
    matched_percent[i] = 100*matched_refls/total_refls

    resolutions = np.asarray(im_pre['Resolution'])

    quintiles = np.quantile(resolutions, [0, 0.2, 0.4, 0.6, 0.8, 1.0])
    quintiles_by_frame.append(quintiles)

quintiles_by_frame = np.asarray(quintiles_by_frame)

plt.figure()
plt.hist(quintiles_by_frame[:,0])
plt.xlabel('Resolution')
plt.ylabel('# Reflections')
plt.title('0% Quintile Distribution')
plt.show()

plt.figure()
plt.hist(quintiles_by_frame[:,1])
plt.xlabel('Resolution')
plt.ylabel('# Reflections')
plt.title('20% Quintile Distribution')
plt.show()

plt.figure()
plt.plot(frames, quintiles_by_frame[:,0])
plt.plot(frames, quintiles_by_frame[:,1])
plt.plot(frames, quintiles_by_frame[:,2])
plt.plot(frames, quintiles_by_frame[:,3])
plt.plot(frames, quintiles_by_frame[:,4])
plt.plot(frames, quintiles_by_frame[:,5])
plt.xlabel('Frame')
plt.ylabel('Precognition Resolution (A)')
#plt.title('Strong Spot Resolution by Frame')
plt.title('Indexed Spot Resolution by Frame')
ax = plt.gca()
ax.axhline(1, linestyle='--')
ax.set_yscale('log')
plt.show()

plt.figure()
plt.plot(frames, matched_percent)
plt.xlabel('Frame')
plt.ylabel('% Reflections Matched to Precognition')
#plt.title('Matched DIALS Strong Spots by Frame')
plt.title('Matched DIALS Indexed Spots by Frame')
plt.show()
