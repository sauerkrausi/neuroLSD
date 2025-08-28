#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Segmentation and evaluation for SynMarker MIP composites.

Inputs:
- Root folders with MIP TIFFs named like:
  WellA1_Seq0000_1_iN_d38_Ctrl_Basson-488_SYP-568_tubulin-647.nd2 (series 01)_MIP.tif

Outputs per root:
- ./_out/metrics/per_object.csv
- ./_out/metrics/per_image_summary.csv
- ./_out/qc_overlays/*_overlay.png
- ./_out/cleaned/*_cleaned.tif and *_cleaned_bgsub.tif

Defaults locked to our spec.
"""

from __future__ import annotations
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
from tifffile import imread, imwrite
from skimage import filters, exposure
from skimage.filters import threshold_otsu, gaussian
from skimage.filters import threshold_yen, threshold_triangle, threshold_li
from skimage.morphology import remove_small_objects, remove_small_holes, binary_dilation, disk, white_tophat
from skimage.measure import label as sk_label, regionprops_table
from skimage.segmentation import find_boundaries
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

# -----------------------------
# User paths
# -----------------------------
ROOT_DIRS = [
    #Path("/Users/felix/Documents/PostDoc/Harvard/03_LSD-PD/Microscopy/20241225_diff129_iN-iDA_d23-d30-d38/Data/test"),
    # replace with your paths
    # Path("/Users/felix/Documents/PostDoc/Harvard/03_LSD-PD/Microscopy/20241225_diff129_iN-iDA_d23-d30-d38/Data/20241130_diff129_iNiDA_d30_ASAH1_plate1_SynMarker_MIP"),
     Path("/Users/felix/Documents/PostDoc/Harvard/03_LSD-PD/Microscopy/20241225_diff129_iN-iDA_d23-d30-d38/Data/20241130_diff129_iNiDA_d30_Ctrl_plate1_SynMarker_MIP"),
    # Path("/Users/felix/Documents/PostDoc/Harvard/03_LSD-PD/Microscopy/20241225_diff129_iN-iDA_d23-d30-d38/Data/20241224_diff129_iNiDA_d38_ASAH1_plate2_SynMarkers_MIP"),
    # Path("/Users/felix/Documents/PostDoc/Harvard/03_LSD-PD/Microscopy/20241225_diff129_iN-iDA_d23-d30-d38/Data/20241224_diff129_iNiDA_d38_Ctrl_plate2_SynMarkers_MIP"),
]

# -----------------------------
# Fixed channel policy and markers
# -----------------------------
CHANNEL_ORDER_WL = [568, 488, 405, 647]  # 405 is optional
CYTO_KEYS = {"tubulin", "tub", "nef", "nefh", "microtubule"}
SPOT_KEYS = {"eea1", "vps35"}
DIVERSE_KEYS = {"tom20", "syp","syn1", "bassoon"}

# -----------------------------
# Parameters
# -----------------------------
mask_gauss_sigma_px: float = 1.5  # Gaussian blur sigma for mask smoothing (in pixels)
dilation_radius_px: int = int(round(2.0 * mask_gauss_sigma_px))  # Dilation radius for mask (in pixels)
bg_tophat_radius_px: int = 12  # Radius for white tophat background subtraction (in pixels)
# Optional per-stain cap for spot areas (in pixels); None disables the cap for that stain
SPOT_MAX_AREA_PX_PER_STAIN = {
    "eea1": None,     # allow larger EEA1 foci (do not cap at the global 300 px)
    "vps35": 300,     # keep tight for VPS35
}
spot_min_area_px: int = 3  # Minimum area for spot segmentation (in pixels)
spot_max_area_px: int = 300  # Maximum area for spot segmentation (in pixels)
diverse_min_area_px: int = 3  # Minimum area for diverse marker segmentation (in pixels)
log_scales_px: List[float] = [1.0, 1.6, 2.4, 3.2, 4.5]  # expanded with a larger scale to capture bigger foci
# Set pixel size: 6.1538 pixels = 1 micron, so 1 px = 1/6.1538 micron
pixel_size_um: Optional[float] = 1.0 / 6.1538  # pixel size in microns

# Neuronal mask tuning
mask_threshold_method: str = "li"  # 'li' tends to yield more contiguous FG than Otsu; options: 'li','triangle','yen','otsu'
mask_closing_radius_px: int = 3    # morphological closing radius to connect gaps
mask_hole_area_px: int = 256       # fill holes up to this area
mask_min_object_area_px: int = 64  # drop tiny specks below this area

# Bassoon-specific parameters
BASSOON_BG_TOPHAT_RADIUS_PX = 24  # Stronger background subtraction for Bassoon
BASSOON_TOPHAT_RADIUS_PX = 6      # Additional tophat filter for Bassoon
BASSOON_MAX_AREA_UM2 = 100        # Reject Bassoon objects larger than 100 µm²

# Channel-specific size filters (in µm²)
SIZE_FILTERS_UM2 = {
    "bassoon":   {"min": 0.5, "max": 10},
    "syp":       {"min": 2,    "max": 450},
    "syn1":      {"min": 2,    "max": 450},
    "eea1":      {"min": 3,    "max": None},  
    "tubulin":   {"min": 5,    "max": None},
    "nefh":      {"min": 5,    "max": None},
}

def apply_size_filter(df: pd.DataFrame, stain_name: str) -> pd.DataFrame:
    # Apply min/max area_um2 filter for given stain_name if defined
    key = stain_name.lower()
    if pixel_size_um is None or "area_um2" not in df.columns:
        return df
    filt = SIZE_FILTERS_UM2.get(key)
    if not filt:
        return df
    if filt["min"] is not None:
        df = df[df["area_um2"] >= filt["min"]]
    if filt["max"] is not None:
        df = df[df["area_um2"] <= filt["max"]]
    return df

def segment_diverse(img_float: np.ndarray, stain_name: str = "") -> np.ndarray:
    # Segment 'diverse' markers using Gaussian smoothing and Otsu threshold.
    # For Bassoon, apply additional tophat filtering.
    img_proc = img_float
    if stain_name.lower() == "bassoon":
        img_proc = white_tophat(img_proc, footprint=disk(BASSOON_TOPHAT_RADIUS_PX))
    sm = gaussian(img_proc, sigma=1.0, preserve_range=True)
    thr = threshold_otsu(sm)
    bw = sm > thr
    # For SYN1, reject small objects more aggressively (handled by size filter, so use diverse_min_area_px here)
    bw = remove_small_objects(bw, min_size=diverse_min_area_px)
    bw = remove_small_holes(bw, area_threshold=16)
    lbl = sk_label(bw)
    return lbl

# -----------------------------
# Filename parsing
# -----------------------------
# Example:
# WellA1_Seq0000_1_iN_d38_Ctrl_Basson-488_SYP-568_tubulin-647.nd2 (series 01)_MIP.tif
# This regex extracts metadata from the filename.
FNAME_RE = re.compile(
    r"^Well(?P<well_row>[A-Z])(?P<well_col>\d+)_Seq(?P<seq>\d+)(?:_\d+)?_(?P<neuron_type>iN|iDA)_(?P<day>d\d+)_(?P<genotype>Ctrl|ASAH1)_(?P<stains>.+?)\.nd2\s+\(series\s+(?P<series>\d+)\)_MIP\.tif$",
    re.IGNORECASE,
)

WL_IN_TOKEN_RE = re.compile(r"(?P<name>[A-Za-z0-9]+)-(?P<wl>\d{3})", re.IGNORECASE)

def parse_metadata_from_name(fname: str) -> Optional[Dict]:
    # Parse metadata from the filename using regex, with a permissive fallback.
    name = Path(fname).name

    # 1) Try strict pattern first
    m = FNAME_RE.match(name)
    if m:
        gd = m.groupdict()
        stain_block = gd["stains"]
        stain_tokens = [t for t in stain_block.split("_") if "-" in t]
        stains: List[Tuple[str, int]] = []
        for tok in stain_tokens:
            mm = WL_IN_TOKEN_RE.search(tok)
            if mm:
                nm = mm.group("name")
                wl = int(mm.group("wl"))
                stains.append((nm, wl))
        panel_name = "_".join([f"{n}-{w}" for n, w in stains])
        return {
            "well_id": f"{gd['well_row']}{gd['well_col']}",
            "seq_id": gd["seq"],
            "series_id": gd["series"],
            "neuron_type": gd["neuron_type"],
            "genotype": gd["genotype"],
            "day": gd["day"],
            "stains": stains,
            "panel_name": panel_name,
        }

    # 2) Fallback: permissive parse for variants like 'Seq00103_iDA_d30_...'
    try:
        # Well
        mw = re.search(r"Well([A-Z])(\d+)", name, flags=re.IGNORECASE)
        # Seq (optionally followed by _digits)
        ms = re.search(r"Seq(\d+)(?:_(\d+))?", name, flags=re.IGNORECASE)
        # Neuron type, day, genotype
        mt = re.search(r"\b(iN|iDA)\b", name, flags=re.IGNORECASE)
        md = re.search(r"\b(d\d+)\b", name, flags=re.IGNORECASE)
        mg = re.search(r"\b(Ctrl|ASAH1)\b", name, flags=re.IGNORECASE)
        # Series
        mr = re.search(r"\(series\s+(\d+)\)", name, flags=re.IGNORECASE)

        if not (mw and ms and mt and md and mg and mr):
            return None

        well_id = f"{mw.group(1).upper()}{mw.group(2)}"
        seq_id = ms.group(1)
        series_id = mr.group(1)
        neuron_type = mt.group(1)
        day = md.group(1)
        genotype = mg.group(1)

        # Stain block: text after genotype_ and before '.nd2'
        # Find the last occurrence of genotype token index, then capture until '.nd2'
        gidx = name.lower().rfind(genotype.lower())
        after_g = name[gidx + len(genotype):]
        msb = re.search(r"_(.+?)\.nd2", after_g, flags=re.IGNORECASE)
        stain_block = msb.group(1) if msb else ""
        stain_tokens = [t for t in stain_block.split("_") if "-" in t]
        stains: List[Tuple[str, int]] = []
        for tok in stain_tokens:
            mm = WL_IN_TOKEN_RE.search(tok)
            if mm:
                nm = mm.group("name")
                wl = int(mm.group("wl"))
                stains.append((nm, wl))
        panel_name = "_".join([f"{n}-{w}" for n, w in stains])

        return {
            "well_id": well_id,
            "seq_id": seq_id,
            "series_id": series_id,
            "neuron_type": neuron_type,
            "genotype": genotype,
            "day": day,
            "stains": stains,
            "panel_name": panel_name,
        }
    except Exception:
        return None

def stains_to_wavelengths(stains: List[Tuple[str,int]]) -> List[int]:
    # Given a list of (name, wavelength), return wavelengths in fixed order.
    present = [wl for _, wl in stains]
    ordered = [wl for wl in CHANNEL_ORDER_WL if wl in present]
    return ordered

def detect_channel_axis(img: np.ndarray) -> Optional[int]:
    # Detect which axis is the channel axis (if any).
    if img.ndim == 2:
        return None
    # Consider small channel counts (2, 3, 4) as possible channel axes
    candidate_axes = [ax for ax, n in enumerate(img.shape) if n in (2, 3, 4)]
    # Prefer last axis if included
    for ax in candidate_axes:
        if ax == img.ndim - 1:
            return ax
    return candidate_axes[0] if candidate_axes else None

def extract_channel(img: np.ndarray, ch_axis: int, index_in_axis: int) -> np.ndarray:
    # Extract a single channel from a multi-channel image.
    slicer = [slice(None)] * img.ndim
    slicer[ch_axis] = index_in_axis
    return img[tuple(slicer)]

# -----------------------------
# Utilities
# -----------------------------
def to_float01(img: np.ndarray) -> np.ndarray:
    # Normalize image to float [0,1] using 99.9th percentile.
    img = img.astype(np.float32)
    vmax = np.percentile(img, 99.9) if img.max() > 0 else 1.0
    if vmax <= 0:
        vmax = 1.0
    return np.clip(img / vmax, 0.0, 1.0)

def normalize_and_sum(images: List[np.ndarray]) -> np.ndarray:
    # Normalize and sum a list of images, clipping to [0,1].
    if not images:
        raise ValueError("No images to combine")
    acc = np.zeros_like(to_float01(images[0]), dtype=np.float32)
    for im in images:
        acc += to_float01(im)
    acc = np.clip(acc, 0.0, 1.0)
    return acc

def classify_marker(name_raw: str) -> str:
    # Classify marker as 'spot', 'diverse', or 'other' based on name.
    nm = name_raw.lower()
    if any(k in nm for k in SPOT_KEYS):
        return "spot"
    if any(k in nm for k in DIVERSE_KEYS):
        return "diverse"
    return "other"

def is_cyto_marker(name_raw: str) -> bool:
    # Return True if marker is a cytoskeletal marker.
    nm = name_raw.lower()
    return any(k in nm for k in CYTO_KEYS)

# -----------------------------
# Robust thresholding helper
# -----------------------------
def robust_thresh(img: np.ndarray, method: str = "yen", mad_k: float = 3.0) -> float:
    """
    Return a conservative threshold for bright-foreground images.
    Guards against degenerate inputs (all zeros/constant/NaN) to avoid numerical warnings.
    Supported methods: 'yen', 'triangle', 'li', 'otsu', 'mad'.
    """
    a = np.asarray(img)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return 0.0
    vmin = float(a.min())
    vmax = float(a.max())
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin + 1e-12:
        # constant or invalid image: pick a mid/high threshold so BW is empty unless there is signal
        return vmax * 0.5 if vmax > 0 else 0.0

    m = method.lower()
    thr = None
    with np.errstate(divide="ignore", invalid="ignore"):
        if m == "yen":
            thr = float(threshold_yen(a))
        elif m == "triangle":
            thr = float(threshold_triangle(a))
        elif m == "li":
            thr = float(threshold_li(a))
        elif m == "mad":
            med = float(np.median(a))
            mad = float(np.median(np.abs(a - med))) * 1.4826
            thr = med + mad_k * mad
        else:
            thr = float(threshold_otsu(a))

    if not np.isfinite(thr):
        # robust fallback
        mu = float(a.mean())
        sd = float(a.std())
        thr = mu + 3.0 * sd
        if not np.isfinite(thr):
            thr = vmax * 0.5
    return float(thr)

# -----------------------------
# Masking and preprocessing
# -----------------------------
def build_neuronal_mask(stains: List[Tuple[str,int]],
                        channel_wls_in_order: List[int],
                        img_float_by_wl: Dict[int, np.ndarray]) -> np.ndarray:
    # Build a binary mask for neuronal regions using cytoskeletal channels.
    cyto_imgs = []
    for stain_name, wl in stains:
        if is_cyto_marker(stain_name):
            if wl in img_float_by_wl:
                cyto_imgs.append(img_float_by_wl[wl])
    if not cyto_imgs:
        # fallback to any tubulin-like wavelength present
        for stain_name, wl in stains:
            if "tub" in stain_name.lower() and wl in img_float_by_wl:
                cyto_imgs.append(img_float_by_wl[wl])
    if not cyto_imgs:
        # final fallback: use last channel if nothing matched
        any_wl = channel_wls_in_order[-1]
        cyto_imgs = [img_float_by_wl[any_wl]]
    mask_img = normalize_and_sum(cyto_imgs)
    blurred = gaussian(mask_img, sigma=mask_gauss_sigma_px, preserve_range=True)
    # More contiguous threshold than Otsu
    thr = robust_thresh(blurred, method=mask_threshold_method)
    bw = blurred > thr
    # Connect small gaps
    if mask_closing_radius_px and mask_closing_radius_px > 0:
        bw = binary_dilation(bw, footprint=disk(mask_closing_radius_px))
        bw = binary_dilation(bw, footprint=disk(mask_closing_radius_px))  # two passes approximate closing
    # Fill holes more aggressively
    bw = remove_small_holes(bw, area_threshold=mask_hole_area_px)
    # Remove tiny specks
    bw = remove_small_objects(bw, min_size=mask_min_object_area_px)
    # Final gentle dilation for safety margin
    if dilation_radius_px and dilation_radius_px > 0:
        bw = binary_dilation(bw, footprint=disk(dilation_radius_px))
    return bw.astype(bool)

def cleaned_and_bgsub(img_float: np.ndarray, mask_bool: np.ndarray, bg_tophat_radius: int = None) -> Tuple[np.ndarray, np.ndarray]:
    # Apply mask and perform background subtraction using white tophat.
    if bg_tophat_radius is None:
        bg_tophat_radius = bg_tophat_radius_px
    cleaned = np.where(mask_bool, img_float, 0.0).astype(np.float32)
    # background subtraction on cleaned
    cleaned_8 = (np.clip(cleaned, 0.0, 1.0) * 255.0).astype(np.uint8)
    bg_removed = white_tophat(cleaned_8, footprint=disk(bg_tophat_radius)).astype(np.float32) / 255.0
    return cleaned, bg_removed

# -----------------------------
# Segmentation
# -----------------------------
def segment_spot(img_float: np.ndarray, stain_name: str = "") -> np.ndarray:
    """
    Segment 'spot' markers using multiscale white-tophat + Yen threshold.
    - VPS35: clip low intensities, build a foreground mask with Otsu, then apply multiscale tophat inside FG.
    - Others (EEA1, SYP, etc.): multiscale tophat over the full image.
    """
    nm = stain_name.lower().strip()
    if nm == "vps35":
        img_adj = np.copy(img_float)
        img_adj[img_adj < (100.0 / 255.0)] = 0.0
        otsu_thr = threshold_otsu(img_adj)
        fg = img_adj > otsu_thr
        responses = []
        for s in log_scales_px:
            selem = disk(int(round(s)))
            responses.append(white_tophat(img_adj * fg, footprint=selem))
        resp_stack = np.max(np.stack(responses, axis=0), axis=0)
        resp_stack = exposure.rescale_intensity(resp_stack, in_range="image", out_range=(0.0, 1.0))
        thr = robust_thresh(resp_stack, method="yen")
        bw = (resp_stack > thr) & fg
    else:
        responses = []
        for s in log_scales_px:
            selem = disk(int(round(s)))
            responses.append(white_tophat(img_float, footprint=selem))
        resp_stack = np.max(np.stack(responses, axis=0), axis=0)
        resp_stack = exposure.rescale_intensity(resp_stack, in_range="image", out_range=(0.0, 1.0))
        thr = robust_thresh(resp_stack, method="yen")
        bw = resp_stack > thr

    bw = remove_small_objects(bw, min_size=spot_min_area_px)
    lbl = sk_label(bw)
    # Apply a stain-specific maximum area cap if configured
    eff_max = SPOT_MAX_AREA_PX_PER_STAIN.get(nm, spot_max_area_px)
    if eff_max is not None and eff_max > 0:
        rp = regionprops_table(lbl, properties=("label", "area"))
        bad = {l for l, a in zip(rp["label"], rp["area"]) if a > eff_max}
        if bad:
            bw2 = np.copy(lbl)
            for b in bad:
                bw2[bw2 == b] = 0
            lbl = sk_label(bw2 > 0)
    return lbl

def segment_diverse(img_float: np.ndarray, stain_name: str = "") -> np.ndarray:
    # Segment 'diverse' markers using Gaussian smoothing and Otsu threshold.
    # For Bassoon, apply additional tophat filtering.
    img_proc = img_float
    if stain_name.lower() == "bassoon":
        img_proc = white_tophat(img_proc, footprint=disk(BASSOON_TOPHAT_RADIUS_PX))
    sm = gaussian(img_proc, sigma=1.0, preserve_range=True)
    thr = threshold_otsu(sm)
    bw = sm > thr
    # For SYN1, reject small objects more aggressively (handled by size filter, so use diverse_min_area_px here)
    bw = remove_small_objects(bw, min_size=diverse_min_area_px)
    bw = remove_small_holes(bw, area_threshold=16)
    lbl = sk_label(bw)
    return lbl

# -----------------------------
# Measurement
# -----------------------------
def compute_coloc_stats(syp_df, bassoon_df, radius_px=5):
    # Returns: coloc_count, not_coloc_count, percent_coloc
    if syp_df.empty or bassoon_df.empty:
        return 0, len(syp_df), 0.0
    syp_centroids = syp_df[["centroid-0", "centroid-1"]].values
    bassoon_centroids = bassoon_df[["centroid-0", "centroid-1"]].values
    tree = cKDTree(bassoon_centroids)
    dists, idxs = tree.query(syp_centroids, distance_upper_bound=radius_px)
    coloc_mask = dists <= radius_px
    coloc_count = np.count_nonzero(coloc_mask)
    not_coloc_count = len(syp_df) - coloc_count
    percent_coloc = 100.0 * coloc_count / len(syp_df) if len(syp_df) > 0 else 0.0
    return coloc_count, not_coloc_count, percent_coloc

def measure_regions(lbl: np.ndarray, inten_raw: np.ndarray, inten_cleaned: np.ndarray) -> pd.DataFrame:
    # Measure region properties for each labeled object.
    # Returns a DataFrame with area, intensity, and shape metrics.
    if lbl.max() == 0:
        return pd.DataFrame(columns=[
            "label", "area_px", "perimeter_px", "eccentricity", "solidity",
            "feret_diameter_max_px", "mean_intensity_raw", "integrated_intensity_raw",
            "mean_intensity_cleaned", "integrated_intensity_cleaned"
        ])
    props_common = ("label", "area", "perimeter", "eccentricity", "solidity", "centroid")
    # try feret if available
    try:
        df_raw = pd.DataFrame(regionprops_table(lbl, intensity_image=inten_raw,
                                                properties=props_common + ("feret_diameter_max", "intensity_mean")))
        feret_key = "feret_diameter_max"
    except Exception:
        df_raw = pd.DataFrame(regionprops_table(lbl, intensity_image=inten_raw,
                                                properties=props_common + ("major_axis_length", "intensity_mean")))
        df_raw.rename(columns={"major_axis_length": "feret_diameter_max"}, inplace=True)
        feret_key = "feret_diameter_max"
    df_raw.rename(columns={
        "area": "area_px",
        "perimeter": "perimeter_px",
        "intensity_mean": "mean_intensity_raw",
        feret_key: "feret_diameter_max_px",
    }, inplace=True)
    df_clean = pd.DataFrame(regionprops_table(lbl, intensity_image=inten_cleaned,
                                              properties=("label", "intensity_mean")))
    df_clean.rename(columns={"intensity_mean": "mean_intensity_cleaned"}, inplace=True)
    df = df_raw.merge(df_clean, on="label", how="left")
    df["integrated_intensity_raw"] = df["mean_intensity_raw"] * df["area_px"]
    df["integrated_intensity_cleaned"] = df["mean_intensity_cleaned"] * df["area_px"]
    if pixel_size_um is not None:
        df["area_um2"] = (pixel_size_um ** 2) * df["area_px"]
        df["feret_diameter_max_um"] = pixel_size_um * df["feret_diameter_max_px"]
    return df

# -----------------------------
# QC overlay
# -----------------------------
def save_qc_overlay_png(out_png: Path,
                        img_raw: np.ndarray,
                        img_clean: np.ndarray,
                        img_bgsub: np.ndarray,
                        lbl: np.ndarray,
                        outline_thickness_px: int = 2,
                        keep_labels: Optional[np.ndarray] = None) -> None:
    # Save a 2x2 panel PNG showing raw, cleaned, bgsub, and overlay with labels.
    # If keep_labels is provided, only show those labels in the overlay.
    if keep_labels is not None:
        # Mask out all labels not in keep_labels
        mask = np.isin(lbl, keep_labels)
        lbl_overlay = lbl * mask
    else:
        lbl_overlay = lbl
    boundaries = find_boundaries(lbl_overlay, mode="outer")
    thick_boundaries = boundaries
    if outline_thickness_px and outline_thickness_px > 1:
        thick_boundaries = binary_dilation(boundaries, footprint=disk(outline_thickness_px))
    fig, axs = plt.subplots(2, 2, figsize=(8, 8), dpi=600)
    axs[0,0].imshow(img_raw, cmap="gray")
    axs[0,0].set_title("raw")
    axs[0,1].imshow(img_clean, cmap="gray")
    axs[0,1].set_title("cleaned")
    axs[1,0].imshow(img_bgsub, cmap="gray")
    axs[1,0].set_title("cleaned_bgsub")
    axs[1,1].imshow(img_raw, cmap="gray")
    overlay = np.zeros((*thick_boundaries.shape, 4), dtype=float)
    overlay[thick_boundaries] = [1, 0, 0, 0.4]  # semi-transparent red (RGBA)
    axs[1,1].imshow(overlay)
    axs[1,1].set_title("overlay labels")
    for ax in axs.ravel():
        ax.set_axis_off()
    plt.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)

def save_syp_bassoon_overlay_png(out_png: Path, syp_img: np.ndarray, bassoon_img: np.ndarray, dpi: int = 300):
    # Save a color overlay PNG: SYP (green), Bassoon (magenta)
    # Inputs are float images [0,1]
    rgb = np.zeros((*syp_img.shape, 3), dtype=np.float32)
    rgb[..., 1] = np.clip(syp_img, 0, 1)  # green
    rgb[..., 0] = np.clip(bassoon_img, 0, 1)  # red
    rgb[..., 2] = np.clip(bassoon_img, 0, 1)  # blue (magenta = R+B)
    plt.figure(figsize=(syp_img.shape[1]/dpi, syp_img.shape[0]/dpi), dpi=dpi)
    plt.imshow(rgb)
    plt.axis('off')
    plt.tight_layout(pad=0)
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight', pad_inches=0)
    plt.close()

# -----------------------------
# Main per file
# -----------------------------
def process_file(tif_path: Path, out_root: Path,
                 per_object_rows: List[Dict],
                 per_image_rows: List[Dict]) -> None:
    # Process a single TIFF file: parse metadata, extract channels, mask, segment, measure, and save outputs.
    meta = parse_metadata_from_name(tif_path.name)
    if meta is None:
        print(f"[skip] name did not match: {tif_path.name}")
        return

    # out_clean_dir = out_root / "cleaned"
    out_qc_dir = out_root / "qc_overlays"
    # out_clean_dir.mkdir(parents=True, exist_ok=True)
    out_qc_dir.mkdir(parents=True, exist_ok=True)

    img = imread(str(tif_path))
    ch_axis = detect_channel_axis(img)
    if ch_axis is None:
        print(f"[warn] single channel image: {tif_path.name}")
        img = img[..., np.newaxis]
        ch_axis = img.ndim - 1

    # build wavelength ordered list present in this image based on stains
    wl_present = stains_to_wavelengths(meta["stains"])
    # If there are multiple axes and one matches the expected number of wavelengths, prefer it
    if img.ndim >= 3:
        expected_n = len(wl_present)
        matching_axes = [ax for ax, n in enumerate(img.shape) if n == expected_n and expected_n in (2, 3, 4)]
        if matching_axes:
            # prefer last axis among matches
            ch_axis = matching_axes[-1]

    # Validate against actual channel count; handle single-channel exports gracefully
    n_channels = img.shape[ch_axis] if ch_axis is not None else 1
    if n_channels != len(wl_present):
        print(f"[warn] {tif_path.name}: channel count mismatch (image has {n_channels}, name suggests {len(wl_present)}). Using first {n_channels} wavelengths from fixed order.")
        wl_present = wl_present[:n_channels]
        # Filter stains to those wavelengths we kept, preserving order
        meta["stains"] = [(n, w) for (n, w) in meta["stains"] if w in set(wl_present)]

    # slice channels into float [0,1]
    img_float_by_wl: Dict[int, np.ndarray] = {}
    if ch_axis is None:
        # true single-channel image; map the sole channel to the first wavelength in wl_present
        single_wl = wl_present[0] if wl_present else CHANNEL_ORDER_WL[-1]
        img_float_by_wl[single_wl] = to_float01(img[..., 0] if img.ndim == 3 else img)
    else:
        for idx, wl in enumerate(wl_present):
            ch_img = extract_channel(img, ch_axis, idx)
            img_float_by_wl[wl] = to_float01(ch_img)

    # neuronal mask
    neuron_mask = build_neuronal_mask(meta["stains"], wl_present, img_float_by_wl)

    # Store per-stain object DataFrames for coloc
    stain_obj_dfs = {}
    # Store per-stain images for overlay
    stain_img_floats = {}

    # iterate stains that are present and in scope
    for stain_name, wl in meta["stains"]:
        if wl not in img_float_by_wl:
            continue
        # Skip gal3
        if stain_name.lower() == "gal3":
            continue
        marker_class = classify_marker(stain_name)
        # Tubulin and NEFH: only background subtract, no mask
        if stain_name.lower() in {"tubulin", "nefh"}:
            raw_f = img_float_by_wl[wl]
            cleaned_f = raw_f
            cleaned_bg = white_tophat((np.clip(raw_f, 0.0, 1.0) * 255.0).astype(np.uint8),
                                      footprint=disk(bg_tophat_radius_px)).astype(np.float32) / 255.0
        # Bassoon-specific background subtraction
        elif stain_name.lower() == "bassoon":
            raw_f = img_float_by_wl[wl]
            cleaned_f, cleaned_bg = cleaned_and_bgsub(raw_f, neuron_mask, bg_tophat_radius=BASSOON_BG_TOPHAT_RADIUS_PX)
        else:
            raw_f = img_float_by_wl[wl]
            cleaned_f, cleaned_bg = cleaned_and_bgsub(raw_f, neuron_mask)

        # Always segment EEA1 as spot (not diverse)
        if marker_class == "spot" or stain_name.lower() == "eea1":
            lbl = segment_spot(cleaned_bg, stain_name=stain_name)
        else:
            lbl = segment_diverse(cleaned_bg, stain_name=stain_name)
        # measurements
        df_objs = measure_regions(lbl, inten_raw=raw_f, inten_cleaned=cleaned_bg)
        df_objs = apply_size_filter(df_objs, stain_name)
        # Store for coloc if SYP or Bassoon
        if stain_name.lower() in {"syp", "bassoon"}:
            stain_obj_dfs[stain_name.lower()] = df_objs.copy()
        # Store for overlay if SYP or Bassoon
        if stain_name.lower() in {"syp", "bassoon"}:
            stain_img_floats[stain_name.lower()] = img_float_by_wl[wl]
        # Only keep labels that passed the filter for overlay
        keep_labels = df_objs["label"].values if not df_objs.empty else None
        # add metadata columns
        if not df_objs.empty:
            df_objs.insert(0, "file_name", tif_path.name)
            df_objs.insert(1, "well_id", meta["well_id"])
            df_objs.insert(2, "neuron_type", meta["neuron_type"])
            df_objs.insert(3, "genotype", meta["genotype"])
            df_objs.insert(4, "day", meta["day"])
            df_objs.insert(5, "panel_name", meta["panel_name"])
            df_objs.insert(6, "stain_name", stain_name)
            df_objs.insert(7, "wavelength", wl)
            df_objs.insert(8, "series_id", meta["series_id"])
            df_objs.insert(9, "dataset_dir", str(tif_path.parent))
            per_object_rows.extend(df_objs.to_dict(orient="records"))

        per_image_rows.append({
            "file_name": tif_path.name,
            "well_id": meta["well_id"],
            "neuron_type": meta["neuron_type"],
            "genotype": meta["genotype"],
            "day": meta["day"],
            "panel_name": meta["panel_name"],
            "stain_name": stain_name,
            "wavelength": wl,
            "series_id": meta["series_id"],
            "dataset_dir": str(tif_path.parent),
            "object_count": int(lbl.max()),
        })

        # QC overlay: only show labels that passed size filter
        base = tif_path.stem.replace(" ", "_")
        save_qc_overlay_png(
            out_qc_dir / f"{base}_{stain_name}-{wl}_overlay.png",
            img_raw=raw_f, img_clean=cleaned_f, img_bgsub=cleaned_bg, lbl=lbl,
            keep_labels=keep_labels
        )

    # After all stains: compute SYP/Bassoon coloc if both present
    if "syp" in stain_obj_dfs and "bassoon" in stain_obj_dfs:
        syp_df = stain_obj_dfs["syp"]
        bassoon_df = stain_obj_dfs["bassoon"]
        coloc_count, not_coloc_count, percent_coloc = compute_coloc_stats(syp_df, bassoon_df, radius_px=5)
        # Add coloc stats to per_image_rows for SYP and Bassoon
        for row in per_image_rows:
            if row["stain_name"].lower() == "syp":
                row["syp_bassoon_coloc_count"] = coloc_count
                row["syp_bassoon_not_coloc_count"] = not_coloc_count
                row["syp_bassoon_coloc_percent"] = percent_coloc
            if row["stain_name"].lower() == "bassoon":
                row["syp_bassoon_coloc_count"] = coloc_count
                row["syp_bassoon_not_coloc_count"] = not_coloc_count
                row["syp_bassoon_coloc_percent"] = percent_coloc
        # Save color overlay PNG of SYP (green) and Bassoon (magenta)
        if "syp" in stain_img_floats and "bassoon" in stain_img_floats:
            base = tif_path.stem.replace(" ", "_")
            out_png = out_root / "qc_overlays" / f"{base}_SYP-Bassoon_color_overlay.png"
            save_syp_bassoon_overlay_png(
                out_png,
                syp_img=stain_img_floats["syp"],
                bassoon_img=stain_img_floats["bassoon"],
                dpi=300
            )

# -----------------------------
# Runner
# -----------------------------
def main():
    # Main entry point: process all root directories, aggregate results, and write CSVs.
    if not ROOT_DIRS:
        print("Set ROOT_DIRS to your four dataset folders.")
        return
    all_object_rows: List[Dict] = []
    all_image_rows: List[Dict] = []
    for root in ROOT_DIRS:
        if not root.exists():
            print(f"[skip] missing root: {root}")
            continue
        out_root = root / "_out"
        out_root.mkdir(parents=True, exist_ok=True)

        per_object_rows: List[Dict] = []
        per_image_rows: List[Dict] = []

        tifs = sorted([p for p in root.glob("*.tif") if p.name.endswith("_MIP.tif") and not p.name.startswith("._")])
        if not tifs:
            print(f"[warn] no MIP tifs in {root}")
        for p in tifs:
            try:
                process_file(p, out_root, per_object_rows, per_image_rows)
            except Exception as e:
                print(f"[error] {p.name}: {e}")

        # write per-root CSVs
        metrics_dir = out_root / "metrics"
        metrics_dir.mkdir(parents=True, exist_ok=True)
        if per_object_rows:
            df_obj = pd.DataFrame(per_object_rows)
            df_obj.to_csv(metrics_dir / "per_object.csv", index=False)
        else:
            df_obj = pd.DataFrame()
        if per_image_rows:
            df_img = pd.DataFrame(per_image_rows)
            # Add per-image averages of numeric columns from per_object
            if not df_obj.empty:
                group_keys = ["file_name", "stain_name", "wavelength"]
                # Exclude group keys from numeric columns before aggregation
                numeric_cols = [c for c in df_obj.select_dtypes(include=[np.number]).columns if c not in group_keys]
                if numeric_cols:
                    agg_df = df_obj.groupby(group_keys, as_index=False)[numeric_cols].mean()
                    # Prefix only the aggregated value columns
                    agg_df = agg_df.rename(columns={col: f"avg_{col}" for col in numeric_cols})
                    # Merge averages into per_image_summary
                    df_img = df_img.merge(agg_df, on=group_keys, how="left")
            df_img.to_csv(metrics_dir / "per_image_summary.csv", index=False)

        all_object_rows.extend(per_object_rows)
        all_image_rows.extend(per_image_rows)

    # combined CSV across roots (next to first root if present)
    if ROOT_DIRS:
        combined_base = ROOT_DIRS[0].parent
    else:
        combined_base = Path(".")
    if all_object_rows:
        pd.DataFrame(all_object_rows).to_csv(combined_base / "combined_per_object.csv", index=False)
    if all_image_rows:
        pd.DataFrame(all_image_rows).to_csv(combined_base / "combined_per_image_summary.csv", index=False)
    print("Done.")

if __name__ == "__main__":
    main()