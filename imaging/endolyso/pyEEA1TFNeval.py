#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
from pathlib import Path
import tifffile as tiff
import numpy as np
import pandas as pd
from skimage import filters, util
from skimage.measure import label as sk_label, regionprops_table
from skimage.morphology import remove_small_objects, disk, white_tophat
import matplotlib.pyplot as plt
from tqdm import tqdm

# =============================================================
# Configuration
# =============================================================
# local path
#root_dir = Path("/Users/felix/Desktop/channels")

root_dir = Path("/Volumes/Expansion/Harvard/NIC/LucilleTwo/Project_neuroLSD/20250121_HeLa-Ctrl-ASAH1_TFN-488_uptakeAssay/channels")
rep_dirs = [root_dir / f"rep{i}" for i in range(1, 4)]
output_dir = Path("/Volumes/Expansion/Harvard/NIC/LucilleTwo/Project_neuroLSD/20250121_HeLa-Ctrl-ASAH1_TFN-488_uptakeAssay/pyeval")
qc_dir = output_dir / "qc_overlays"
metrics_dir = output_dir / "metrics"
for d in [output_dir, qc_dir, metrics_dir]:
    d.mkdir(parents=True, exist_ok=True)

# =============================================================
# Helper functions
# =============================================================

fname_pat = re.compile(
    r"^(?P<sample>.+_MIP\.tif)(?P<channel>EEA1|TFN)-(?P<idx>\d+)\.tif$"
)

def parse_pair_filename(fname):
    m = fname_pat.match(fname)
    if m:
        return m.groupdict()
    return None

def preprocess(img, channel_name=None):
    img = np.asarray(img)
    img_8bit = util.img_as_ubyte(img / img.max() if img.max() > 0 else img)
    # Use larger disk for background subtraction for EEA1/TFN, else default
    if channel_name == "EEA1":
        background = filters.rank.mean(img_8bit, disk(15))
        img_bgsub = img_8bit.astype(np.int32) - background.astype(np.int32)
        img_bgsub = np.clip(img_bgsub, 0, 255).astype(np.uint8)
        img_sharp = white_tophat(img_bgsub, footprint=disk(3))
        return img_sharp
    elif channel_name == "TFN":
        background = filters.rank.mean(img_8bit, disk(10))
        img_bgsub = img_8bit.astype(np.int32) - background.astype(np.int32)
        img_bgsub = np.clip(img_bgsub, 0, 255).astype(np.uint8)
        img_sharp = white_tophat(img_bgsub, footprint=disk(3))
        return img_sharp
    else:
        background = filters.rank.mean(img_8bit, disk(15))
        img_bgsub = img_8bit.astype(np.int32) - background.astype(np.int32)
        img_bgsub = np.clip(img_bgsub, 0, 255).astype(np.uint8)
        img_sharp = white_tophat(img_bgsub, footprint=disk(3))
        return img_sharp

def segment(img_proc, min_area=30, max_area=2000):
    thresh = filters.threshold_otsu(img_proc)
    bw = img_proc > thresh
    bw = remove_small_objects(bw, min_size=min_area)
    seg = sk_label(bw)
    props = regionprops_table(seg, properties=("label", "area", "centroid"))
    valid_labels = [lbl for lbl, area in zip(props["label"], props["area"]) if min_area < area < max_area]
    filtered = np.zeros_like(seg)
    for lbl in valid_labels:
        filtered[seg == lbl] = lbl
    seg = sk_label(filtered > 0)
    return seg, props

# =============================================================
# Main Processing: Only EEA1/TFN
# =============================================================

summary_rows = []
overlap_rows = []

for rep_dir in rep_dirs:
    rep_summary_rows = []
    rep_overlap_rows = []
    all_tifs = [p for p in rep_dir.glob("*.tif") if not p.name.startswith("._")]
    parsed = [ (p, parse_pair_filename(p.name)) for p in all_tifs ]
    parsed = [ (p, meta) for p, meta in parsed if meta is not None ]
    from collections import defaultdict
    groups = defaultdict(dict)
    for p, meta in parsed:
        key = (meta["sample"], meta["idx"])
        groups[key][meta["channel"]] = p

    image_items = sorted(groups.items(), key=lambda x: (x[0][0], int(x[0][1])))

    for (sample, idx), chdict in tqdm(image_items, desc=f"Processing {rep_dir.name}"):
        if "EEA1" not in chdict or "TFN" not in chdict:
            continue  # Only process pairs

        img_eea1 = tiff.imread(chdict["EEA1"])
        img_tfn = tiff.imread(chdict["TFN"])
        img_eea1_proc = preprocess(img_eea1, channel_name="EEA1")
        img_tfn_proc = preprocess(img_tfn, channel_name="TFN")

        seg_eea1, props_eea1 = segment(img_eea1_proc)
        seg_tfn, props_tfn = segment(img_tfn_proc)

        n_eea1 = int(seg_eea1.max())
        n_tfn = int(seg_tfn.max())
        area_eea1 = float(np.sum(seg_eea1 > 0))
        area_tfn = float(np.sum(seg_tfn > 0))
        img_area = img_eea1.shape[0] * img_eea1.shape[1]

        cent_eea1 = np.array(list(zip(props_eea1["centroid-0"], props_eea1["centroid-1"]))) if "centroid-0" in props_eea1 else np.zeros((0,2))
        cent_tfn = np.array(list(zip(props_tfn["centroid-0"], props_tfn["centroid-1"]))) if "centroid-0" in props_tfn else np.zeros((0,2))
        COLOC_RADIUS_PX = 500/160  # ~3.125 px

        coloc_pairs = 0
        if len(cent_tfn) > 0:
            for ca in cent_eea1:
                dists = np.sqrt(((ca - cent_tfn)**2).sum(axis=1))
                if np.any(dists <= COLOC_RADIUS_PX):
                    coloc_pairs += 1
        # else: coloc_pairs remains 0

        # Overlay PNG
        rgb = np.zeros((*img_eea1_proc.shape, 3), dtype=np.uint8)
        rgb[..., 0] = img_eea1_proc
        rgb[..., 1] = img_tfn_proc
        fig, ax = plt.subplots(figsize=(8,4))
        ax.imshow(rgb)
        if len(cent_eea1) > 0:
            ax.scatter(cent_eea1[:,1], cent_eea1[:,0], c='r', s=2, label='EEA1')
        if len(cent_tfn) > 0:
            ax.scatter(cent_tfn[:,1], cent_tfn[:,0], c='g', s=2, label='TFN')
        if len(cent_eea1) > 0 and len(cent_tfn) > 0:
            for ca in cent_eea1:
                dists = np.sqrt(((ca - cent_tfn)**2).sum(axis=1))
                close_idx = np.where(dists <= COLOC_RADIUS_PX)[0]
                for idxB in close_idx:
                    ax.plot([ca[1], cent_tfn[idxB][1]], [ca[0], cent_tfn[idxB][0]], 'y-', linewidth=0.5)
        ax.set_axis_off()
        ax.legend(loc='upper right', fontsize=6)
        qc_path = qc_dir / f"{rep_dir.name}_{sample}_series{idx}_overlay.png"
        plt.tight_layout()
        plt.savefig(qc_path, dpi=200)
        plt.close(fig)

        # Save summary
        row = {
            "replicate": rep_dir.name,
            "sample": sample,
            "series": idx,
            "n_eea1": n_eea1,
            "n_tfn": n_tfn,
            "area_eea1": area_eea1,
            "area_tfn": area_tfn,
            "img_area": img_area,
            "n_eea1_per_img": n_eea1 / img_area,
            "n_tfn_per_img": n_tfn / img_area,
            "area_eea1_per_img": area_eea1 / img_area,
            "area_tfn_per_img": area_tfn / img_area,
            "coloc_pairs": coloc_pairs,
            "coloc_pairs_per_img": coloc_pairs / img_area,
        }
        rep_summary_rows.append(row)
        summary_rows.append(row)
        overlap_row = {
            "replicate": rep_dir.name,
            "sample": sample,
            "series": idx,
            "n_eea1": n_eea1,
            "n_tfn": n_tfn,
            "coloc_pairs": coloc_pairs,
            "frac_coloc_eea1": coloc_pairs / n_eea1 if n_eea1 > 0 else 0,
            "frac_coloc_tfn": coloc_pairs / n_tfn if n_tfn > 0 else 0,
        }
        rep_overlap_rows.append(overlap_row)
        overlap_rows.append(overlap_row)

    # Write per-repeat CSVs
    if rep_summary_rows:
        df_rep = pd.DataFrame(rep_summary_rows)
        df_rep.to_csv(metrics_dir / f"per_image_summary_{rep_dir.name}.csv", index=False)
    if rep_overlap_rows:
        df_overlap_rep = pd.DataFrame(rep_overlap_rows)
        df_overlap_rep.to_csv(metrics_dir / f"pairwise_colocalization_{rep_dir.name}.csv", index=False)

# Write CSVs
if summary_rows:
    df = pd.DataFrame(summary_rows)
    df.to_csv(metrics_dir / "per_image_summary.csv", index=False)
if overlap_rows:
    df_overlap = pd.DataFrame(overlap_rows)
    df_overlap.to_csv(metrics_dir / "pairwise_colocalization.csv", index=False)

print("Done. Results in:", metrics_dir)