#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
from pathlib import Path
import tifffile as tiff
import numpy as np
import pandas as pd
import math
from skimage import filters, util
from skimage.measure import label as sk_label, regionprops_table
from skimage.morphology import remove_small_objects, disk, white_tophat
from skimage.restoration import denoise_wavelet
import torch
from tqdm import tqdm
import matplotlib.pyplot as plt
from cellpose import plot, models, core

# =============================================================
# Configuration
# =============================================================

# path for local testing on SSD
#input_dir = Path("/Users/felix/HMS Dropbox/Felix Kraus/Felix/Harvard/03_LSD-PD/Microscopy/20250613_HeLa-Ctrl-ASAH1_EndoLysoStains/cpsamEL/split_channels_renamed")
#output_dir = Path("/Users/felix/HMS Dropbox/Felix Kraus/Felix/Harvard/03_LSD-PD/Microscopy/20250613_HeLa-Ctrl-ASAH1_EndoLysoStains/cpsamEL/out_masks")

# remote HDD paths:
#input_dir = Path("/Volumes/Expansion/Harvard/NIC/LucilleTwo/Project_neuroLSD/20250613_HeLa-Ctrl-ASAH1_EndoLysoStains/split_channels_renamed/VPS35-488_EEA1-568_HA-647")
#input_dir = Path("/Volumes/Expansion/Harvard/NIC/LucilleTwo/Project_neuroLSD/20250613_HeLa-Ctrl-ASAH1_EndoLysoStains/split_channels_renamed/VPS35-488_LAMP1-568_HA-647")
input_dir = Path("/Volumes/Expansion/Harvard/NIC/LucilleTwo/Project_neuroLSD/20250613_HeLa-Ctrl-ASAH1_EndoLysoStains/split_channels_renamed/VPS35-488_RAB7-568_HA-647")
output_dir = Path("/Volumes/Expansion/Harvard/NIC/LucilleTwo/Project_neuroLSD/20250613_HeLa-Ctrl-ASAH1_EndoLysoStains/out_masks")
metrics_dir = output_dir / "metrics"
qc_dir = output_dir / "qc_overlays"
filtered_tiff_dir = output_dir / "filtered_organelle_tiffs"

for d in [output_dir, metrics_dir, qc_dir, filtered_tiff_dir]:
    d.mkdir(parents=True, exist_ok=True)

device = torch.device("mps")
use_gpu = core.use_gpu()
nuclei_model = models.CellposeModel(gpu=use_gpu, model_type='cpsam')

# Channel configuration and area filters (diameter in px)
CHANNELS = {
    "DNA":    {"min_diam": 100, "max_diam": 1000, "model": nuclei_model},
    "HA":     {"min_diam": 6,   "max_diam": 100},
    "EEA1":   {"min_diam": 3,   "max_diam": 25},
    "VPS35":  {"min_diam": 3,   "max_diam": 25},
    "RAB7":  {"min_diam": 3,   "max_diam": 25},
    "LAMP1":  {"min_diam": 6,   "max_diam": 100},
}
for ch in CHANNELS.values():
    ch["min_area"] = math.pi * (ch["min_diam"]/2)**2
    ch["max_area"] = math.pi * (ch["max_diam"]/2)**2



# =============================================================
# Helper functions
# Image file name parsing, flot conversion of matplotlib-compaible plotting, Image preprocessing
# adjust disk-radius for modulation of pre-filtering. Larger is harsher C
# =============================================================
regex_pat = re.compile(r"^(?P<prefix>.+_MIP)_(?P<channel>[A-Za-z0-9+-]+)-(?P<wavelength>\d{3,4})\.tif$")

def parse_filename(fname: str):
    m = regex_pat.match(fname)
    if m:
        d = m.groupdict()
        d['field'] = 1
        return d
    return None

def prepare_img_for_plot(img):
    if img is None:
        return img
    img = np.asarray(img)
    img_min, img_max = img.min(), img.max()
    if img_max > img_min:
        img_norm = (img - img_min) / (img_max - img_min)
    else:
        img_norm = np.zeros_like(img, dtype=np.float32)
    img_uint8 = (img_norm * 255).astype(np.uint8)
    return img_uint8

def preprocess_image(img, channel_name):
    img = np.asarray(img)
    img_8bit = util.img_as_ubyte(img / img.max() if img.max() > 0 else img)
    if channel_name == "DNA":
        background = filters.rank.mean(img_8bit, disk(10))
        img_bgsub = img_8bit.astype(np.int32) - background.astype(np.int32)
        img_bgsub = np.clip(img_bgsub, 0, 255).astype(np.uint8)
        return img_bgsub
    elif channel_name == "EEA1":
        background = filters.rank.mean(img_8bit, disk(25))
        img_bgsub = img_8bit.astype(np.int32) - background.astype(np.int32)
        img_bgsub = np.clip(img_bgsub, 0, 255).astype(np.uint8)
        img_sharp = white_tophat(img_bgsub, footprint=disk(3))
        return img_sharp
    elif channel_name == "VPS35":
        # Denoise before filtering for VPS35
        # Remove deprecated multichannel/convert2ycbcr, use channel_axis=None for grayscale
        img_denoised = denoise_wavelet(img_8bit, channel_axis=None, rescale_sigma=True)
        img_denoised = (img_denoised * 255).astype(np.uint8)
        background = filters.rank.mean(img_denoised, disk(25))
        img_bgsub = img_denoised.astype(np.int32) - background.astype(np.int32)
        img_bgsub = np.clip(img_bgsub, 0, 255).astype(np.uint8)
        img_sharp = white_tophat(img_bgsub, footprint=disk(3))
        return img_sharp
    else:
        background = filters.rank.mean(img_8bit, disk(15))
        img_bgsub = img_8bit.astype(np.int32) - background.astype(np.int32)
        img_bgsub = np.clip(img_bgsub, 0, 255).astype(np.uint8)
        img_sharp = white_tophat(img_bgsub, footprint=disk(3))
        return img_sharp



# =============================================================
# Main Processing
# =============================================================

all_tifs = [p for p in input_dir.glob("*.tif") if not p.name.startswith("._")]
if len(all_tifs) == 0:
    print(f"No .tif files found in {input_dir}. Did you split your composite TIFFs?")
else:
    print(f"Found {len(all_tifs)} .tif files in {input_dir}")
    for p in all_tifs[:5]:
        print(f"Example file: {p.name}")

parsed = [ (p, parse_filename(p.name)) for p in all_tifs ]
parsed = [ (p, meta) for p, meta in parsed if meta is not None ]
if len(parsed) == 0:
    print("No files matched the expected filename pattern. Check your regex and filenames.")
else:
    print(f"Parsed {len(parsed)} files with expected pattern.")

from collections import defaultdict
groups = defaultdict(dict)
all_channels_found = set()
for p, meta in parsed:
    if meta is None:
        continue
    prefix = meta['prefix']
    field = int(meta['field'])
    raw_channel = meta['channel']
    groups[(prefix, field)][raw_channel] = p
    all_channels_found.add(raw_channel)

present_stains = [ch for ch in CHANNELS if ch in all_channels_found]
print("Detected stains in dataset:", present_stains)

summary_rows = []
overlap_rows = []

image_items = sorted(groups.items(), key=lambda x: (x[0][0], x[0][1]))

for (prefix, field), chdict in tqdm(image_items, desc="Processing fields"):
    results = {
        'prefix': Path(prefix).name,
        'field': field,
    }
    nuclei_count = np.nan
    masks = {}
    avg_sizes = {}

    # --- Nucleus segmentation with Cellpose ---
    if 'DNA' in chdict:
        img_nuclei = tiff.imread(chdict['DNA'])
        img_nuclei_proc = preprocess_image(img_nuclei, "DNA")
        nmasks, nflows, _ = nuclei_model.eval(img_nuclei_proc, flow_threshold=0.4, cellprob_threshold=0.0, channels=[0,0])
        nmasks = sk_label(nmasks > 0)
        # Area filtering by diameter
        props = regionprops_table(nmasks, properties=("label", "area"))
        valid_labels = [lbl for lbl, area in zip(props["label"], props["area"]) if CHANNELS["DNA"]["min_area"] < area < CHANNELS["DNA"]["max_area"]]
        filtered = np.zeros_like(nmasks)
        for lbl in valid_labels:
            filtered[nmasks == lbl] = lbl
        nmasks = sk_label(filtered > 0)
        nuclei_count = int(nmasks.max())
        results["n_nuclei"] = nuclei_count

        # Nucleus QC overlay
        img_plot = prepare_img_for_plot(img_nuclei_proc)
        if (
            img_plot is not None
            and img_plot.size > 0
            and np.isfinite(img_plot).all()
            and img_plot.ndim == 2
            and img_plot.shape[0] > 1
            and img_plot.shape[1] > 1
            and img_plot.dtype in (np.uint8, np.uint16, np.float32, np.float64)
        ):
            fig, ax = plt.subplots(figsize=(8,4))
            ax.imshow(img_plot, cmap='gray')
            ax.contour(nmasks, colors='r', linewidths=0.5)
            ax.set_axis_off()
            qc_path = qc_dir / f"{Path(prefix).name}_F{field:03d}_DNA_overlay.png"
            plt.tight_layout()
            plt.savefig(qc_path, dpi=300)
            plt.close(fig)
        else:
            print(f"Skipping overlay for {prefix} DNA: invalid image for plotting.")

    # --- Organelle segmentation (simple thresholding) ---
    for stain in present_stains:
        if stain == "DNA" or stain not in chdict:
            continue
        img = tiff.imread(chdict[stain])
        img_proc = preprocess_image(img, stain)
        # Otsu threshold
        thresh = filters.threshold_otsu(img_proc)
        bw = img_proc > thresh
        bw = remove_small_objects(bw, min_size=int(CHANNELS[stain]["min_area"]))
        seg = sk_label(bw)
        props = regionprops_table(seg, properties=("label", "area", "centroid"))
        valid_labels = [lbl for lbl, area in zip(props["label"], props["area"]) if CHANNELS[stain]["min_area"] < area < CHANNELS[stain]["max_area"]]
        filtered = np.zeros_like(seg)
        valid_areas = []
        for lbl, area in zip(props["label"], props["area"]):
            if lbl in valid_labels:
                filtered[seg == lbl] = lbl
                valid_areas.append(area)
        seg = sk_label(filtered > 0)
        n_objs = int(seg.max())
        total_area = float(np.sum(seg > 0))
        results[f"n_{stain}"] = n_objs
        results[f"area_{stain}"] = total_area
        results[f"{stain}_per_nucleus"] = n_objs / nuclei_count if nuclei_count and nuclei_count > 0 else np.nan
        masks[stain] = seg
        # Add average size for this stain
        results[f"avg_area_{stain}"] = np.mean(valid_areas) if valid_areas else np.nan

        # Overlay
        img_plot = prepare_img_for_plot(img_proc)
        if (
            img_plot is not None
            and img_plot.size > 0
            and np.isfinite(img_plot).all()
            and img_plot.ndim == 2
            and img_plot.shape[0] > 1
            and img_plot.shape[1] > 1
            and img_plot.dtype in (np.uint8, np.uint16, np.float32, np.float64)
        ):
            fig, ax = plt.subplots(figsize=(8,4))
            ax.imshow(img_plot, cmap='gray')
            ax.contour(seg, colors='r', linewidths=0.5)
            ax.set_axis_off()
            qc_path = qc_dir / f"{Path(prefix).name}_F{field:03d}_{stain}_overlay.png"
            plt.tight_layout()
            plt.savefig(qc_path, dpi=300)
            plt.close(fig)
        else:
            print(f"Skipping overlay for {prefix} {stain}: invalid image for plotting.")

# --- # Save filtered TIFF for organelle channels only        
        # if stain != "DNA":
        #     filtered_tiff_path = filtered_tiff_dir / f"{Path(prefix).name}_F{field:03d}_{stain}_filtered.tif"
        #     tiff.imwrite(filtered_tiff_path, img_proc.astype(np.uint8))

    # --- Colocalization for all pairs of non-DNA stains ---
    stains_for_coloc = [s for s in present_stains if s != "DNA" and s in masks]
    for i, stainA in enumerate(stains_for_coloc):
        for j, stainB in enumerate(stains_for_coloc):
            if i >= j:
                continue
            maskA = masks[stainA]
            maskB = masks[stainB]
            propsA = regionprops_table(maskA, properties=("label", "centroid"))
            centA = np.array(list(zip(propsA["centroid-0"], propsA["centroid-1"]))) if "centroid-0" in propsA else np.zeros((0,2))
            propsB = regionprops_table(maskB, properties=("label", "centroid"))
            centB = np.array(list(zip(propsB["centroid-0"], propsB["centroid-1"]))) if "centroid-0" in propsB else np.zeros((0,2))

            COLOC_RADIUS_PX = 500/160
            coloc_pairs = 0
            for ca in centA:
                dists = np.sqrt(((ca - centB)**2).sum(axis=1))
                if np.any(dists <= COLOC_RADIUS_PX):
                    coloc_pairs += 1
            coloc_count = coloc_pairs
            n_objs_A = int(maskA.max())
            n_objs_B = int(maskB.max())
            overlap_area = np.sum((maskA > 0) & (maskB > 0))
            total_area_A = np.sum(maskA > 0)
            total_area_B = np.sum(maskB > 0)

            results[f"n_coloc_{stainA}_{stainB}"] = coloc_count

            # Overlay for colocalization pairs
            imgA = tiff.imread(chdict[stainA])
            imgB = tiff.imread(chdict[stainB])
            imgA_proc = prepare_img_for_plot(preprocess_image(imgA, stainA))
            imgB_proc = prepare_img_for_plot(preprocess_image(imgB, stainB))
            rgb = np.zeros((*imgA_proc.shape, 3), dtype=np.uint8)
            rgb[..., 0] = imgA_proc
            rgb[..., 1] = imgB_proc
            fig, ax = plt.subplots(figsize=(8,4))
            ax.imshow(rgb)
            if len(centA) > 0:
                ax.scatter(centA[:,1], centA[:,0], c='r', s=1, label=f'{stainA} centroids')
            if len(centB) > 0:
                ax.scatter(centB[:,1], centB[:,0], c='g', s=1, label=f'{stainB} centroids')
            if len(centA) > 0 and len(centB) > 0:
                for idxA, ca in enumerate(centA):
                    dists = np.sqrt(((ca - centB)**2).sum(axis=1))
                    close_idx = np.where(dists <= COLOC_RADIUS_PX)[0]
                    for idxB in close_idx:
                        ax.plot([ca[1], centB[idxB][1]], [ca[0], centB[idxB][0]], 'y-', linewidth=0.5)
            ax.set_axis_off()
            ax.legend(loc='upper right', fontsize=6)
            qc_path = qc_dir / f"{Path(prefix).name}_F{field:03d}_{stainA}_vs_{stainB}_coloc_overlay.png"
            plt.tight_layout()
            plt.savefig(qc_path, dpi=300)
            plt.close(fig)

            overlap_rows.append({
                "prefix": Path(prefix).name,
                "field": field,
                "stainA": stainA,
                "stainB": stainB,
                "nuclei_count": nuclei_count,
                "n_objs_A": n_objs_A,
                "n_objs_B": n_objs_B,
                "coloc_count": coloc_count,
                "overlap_area": overlap_area,
                "total_area_A": total_area_A,
                "total_area_B": total_area_B,
                "frac_coloc_A": coloc_count / n_objs_A if n_objs_A > 0 else 0,
                "frac_coloc_B": coloc_count / n_objs_B if n_objs_B > 0 else 0,
                "coloc_per_nucleus": coloc_count / nuclei_count if nuclei_count > 0 else 0,
            })

    summary_rows.append(results)

# =============================================================
# Write summary CSV
# =============================================================
if summary_rows:
    df = pd.DataFrame(summary_rows)
    df.to_csv(metrics_dir / "per_image_summary.csv", index=False)
if overlap_rows:
    df_overlap = pd.DataFrame(overlap_rows)
    df_overlap.to_csv(metrics_dir / "pairwise_colocalization.csv", index=False)

print("Done. Results in:", metrics_dir)