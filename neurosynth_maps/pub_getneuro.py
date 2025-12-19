#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract peak MNI coordinates from NeuroSynth maps.

Steps:
1) Load and average NeuroSynth maps per region.
2) Extract peak MNI coordinates (global max; left/right for TPJ).
3) Save results to CSV.
4) Generate QC plot with markers.

Dependencies: nibabel, numpy, pandas, scipy, nilearn, seaborn, matplotlib
"""

import nibabel as nib
import numpy as np
import pandas as pd
import os, glob
from scipy.ndimage import maximum_filter
from nilearn.image import new_img_like
from nilearn import plotting
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# --- Functions ---
def load_nifti_files(file_paths):
    return [nib.load(fp) for fp in file_paths]

def average_nifti_images(images):
    data = np.mean([img.get_fdata() for img in images], axis=0)
    return new_img_like(images[0], data, images[0].affine)

def find_peak_mni_coords(image, per_hemisphere=False):
    data = image.get_fdata()
    filt = maximum_filter(data, size=3)
    coords = np.array(np.nonzero(data)).T
    mni = nib.affines.apply_affine(image.affine, coords)

    if per_hemisphere:
        left = coords[mni[:,0] < 0]
        right = coords[mni[:,0] >= 0]
        l_val = data[left[:,0], left[:,1], left[:,2]]
        r_val = data[right[:,0], right[:,1], right[:,2]]
        l_coord = nib.affines.apply_affine(image.affine, left[np.argmax(l_val)])
        r_coord = nib.affines.apply_affine(image.affine, right[np.argmax(r_val)])
        return l_coord, r_coord
    else:
        vals = data[filt == data]
        vox = np.array(np.nonzero(filt == data)).T[np.argmax(vals)]
        return nib.affines.apply_affine(image.affine, vox)

# --- Paths ---
path_main = '/Users/administrator/Documents/MATLAB/project_cb_tom/neurosynth_maps_dmn'
regions = {
    'vmpfc': glob.glob(os.path.join(path_main, 'vmpfc', '*.nii.gz')),
    'dmpfc': glob.glob(os.path.join(path_main, 'dmpfc', '*.nii.gz')),
    'precuneus': glob.glob(os.path.join(path_main, 'precuneus', '*.nii.gz')),
    'tpj': glob.glob(os.path.join(path_main, 'tpj', '*.nii.gz'))
}

# --- Process ---
results = []
for name, files in regions.items():
    imgs = load_nifti_files(files)
    avg = average_nifti_images(imgs)

    if name == 'tpj':
        l, r = find_peak_mni_coords(avg, per_hemisphere=True)
        results.append({'Region': 'tpj_left',  'X': l[0], 'Y': l[1], 'Z': l[2]})
        results.append({'Region': 'tpj_right', 'X': r[0], 'Y': r[1], 'Z': r[2]})
    else:
        peak = find_peak_mni_coords(avg)
        results.append({'Region': name, 'X': peak[0], 'Y': peak[1], 'Z': peak[2]})

df = pd.DataFrame(results)
df.to_csv(os.path.join(path_main, 'neurosynth_coord.csv'), index=False)
print(df)

# --- QC plot ---
palette = sns.color_palette("Set2", n_colors=len(df))
node_cmap = LinearSegmentedColormap.from_list("node_cmap", palette)
coords = df[['X','Y','Z']].values.tolist()
plotting.plot_markers(range(len(coords)), coords, node_cmap=node_cmap,
                      node_size=200, display_mode='ortho', colorbar=False)
print("Labels:", df['Region'].tolist())
