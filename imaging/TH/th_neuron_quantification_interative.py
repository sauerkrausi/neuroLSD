# %% Imports and setup
import os
from pathlib import Path
import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt
from cellpose import models, io
from skimage.measure import regionprops, label
from skimage.color import label2rgb
from skimage.morphology import remove_small_objects
from skimage.segmentation import clear_border
from skimage import exposure
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import disk, white_tophat, black_tophat
import torch
import pandas as pd
import seaborn as sns

# %% Define paths and device
# --- Set paths
input_root = Path("tiff_channels")
#input_root = Path("test_channels")
output_dir = Path("/Users/felix/HMS Dropbox/Felix Kraus/Felix/Harvard/03_LSD-PD/Microscopy/20240417_diff118_iN-iDA_test_d23_Ctrl-SMPD1-ASAH1/THeval/results")
output_dir.mkdir(parents=True, exist_ok=True)

summary_data = []

# Use MPS on macOS if available
if torch.backends.mps.is_available():
    device = torch.device("mps")
    print("Using MPS backend")
else:
    device = torch.device("cpu")
    print("MPS not available, using CPU")

# %% Load Cellpose models
# Load custom model for MAP2/TH (neuron bodies)
model_path_MAP2 = "/Users/felix/.cellpose/models/FK_MAP2_20250713_v2"
model_path_TH = "/Users/felix/.cellpose/models/cpsam_20250713_TH_v1"
model_path_nuclei = "/Users/felix/.cellpose/models/FK_nuclei_20250713_v1"


model_map2 = models.CellposeModel(pretrained_model=model_path_MAP2, gpu=torch.backends.mps.is_available())
model_th = models.CellposeModel(pretrained_model=model_path_TH, gpu=torch.backends.mps.is_available())
model_nuclei = models.CellposeModel(pretrained_model=model_path_nuclei, gpu=torch.backends.mps.is_available())

# %% Preprocessing function
# --- Image preprocessing function
def preprocess_image(img, sigma=1.0, background_method='white_tophat', disk_size=50):
    """
    Preprocess image with background subtraction and Gaussian smoothing
    
    Parameters:
    - img: input image
    - sigma: standard deviation for Gaussian smoothing
    - background_method: 'white_tophat', 'black_tophat', or 'rolling_ball'
    - disk_size: size of structuring element for morphological operations
    """
    # Ensure image is float for processing
    img_float = img.astype(np.float32)
    
    # Background subtraction
    if background_method == 'white_tophat':
        # White top-hat to remove large background structures
        selem = disk(disk_size)
        background = white_tophat(img_float, selem)
        img_bg_sub = img_float - background
    elif background_method == 'black_tophat':
        # Black top-hat for dark background
        selem = disk(disk_size)
        background = black_tophat(img_float, selem)
        img_bg_sub = img_float + background
    elif background_method == 'rolling_ball':
        # Simple rolling ball background subtraction approximation
        from scipy.ndimage import uniform_filter
        background = uniform_filter(img_float, size=disk_size)
        img_bg_sub = img_float - background
    else:
        img_bg_sub = img_float
    
    # Clip negative values
    img_bg_sub = np.clip(img_bg_sub, 0, img_bg_sub.max())
    
    # Gaussian smoothing
    if sigma > 0:
        img_smooth = gaussian(img_bg_sub, sigma=sigma, preserve_range=True)
    else:
        img_smooth = img_bg_sub
    
    # Convert back to original dtype range
    if img.dtype == np.uint16:
        img_processed = np.clip(img_smooth, 0, 65535).astype(np.uint16)
    else:
        img_processed = np.clip(img_smooth, 0, 255).astype(np.uint8)
    
    return img_processed

# %% Collect image sets
# --- Process each image set
# Find all TH files and group them by condition and base pattern
th_files = list(input_root.glob("*MIP.tifTH-*.tif"))

# Group files by their base pattern
file_groups = {}
for th_file in th_files:
    # Extract the base pattern and index
    th_match = th_file.name
    # Find corresponding MAP2 and nuclei files
    th_index = th_match.split("TH-")[-1].split(".tif")[0]
    
    # Construct patterns for MAP2 and nuclei files
    base_pattern = th_match.replace(f"TH-{th_index}.tif", "")
    map2_file = input_root / f"{base_pattern}MAP2-{th_index}.tif"
    nuclei_file = input_root / f"{base_pattern}nuclei-{th_index}.tif"
    
    if map2_file.exists() and nuclei_file.exists():
        # Extract condition from filename
        name_parts = th_file.name.split("_")
        condition_match = "_".join(name_parts[0:3])  # e.g., "diff118_WellC1_iDA_Ctrl"
        
        file_groups[th_file.stem] = {
            'condition': condition_match,
            'th_file': th_file,
            'map2_file': map2_file,
            'nuclei_file': nuclei_file,
            'base_pattern': base_pattern,
            'index': th_index
        }

print(f"Found {len(file_groups)} complete image sets to process")

# %% Process image sets
for group_name, file_info in file_groups.items():
    print(f"Processing {group_name}...")
    
    # Load individual channel images
    th_img = tiff.imread(file_info['th_file'])
    map2_img = tiff.imread(file_info['map2_file'])
    nuclei_img = tiff.imread(file_info['nuclei_file'])

    
    # %% Preprocess images
    print(f"  Preprocessing images...")
    # Different preprocessing parameters for different channels
    th_img_processed = preprocess_image(th_img, sigma=1.0, background_method='rolling_ball', disk_size=20)
    map2_img_processed = map2_img
    nuclei_img_processed = preprocess_image(nuclei_img, sigma=0.8, background_method='white_tophat', disk_size=20)
    
    condition = file_info['condition']

    # %% Segment MAP2 and TH directly with Cellpose SAM
    masks_map2, _, _ = model_map2.eval(map2_img_processed, diameter=None)
    masks_th, _, _ = model_th.eval(th_img_processed, diameter=None)
    # Filter out small TH signals (<50 px)
    masks_th = remove_small_objects(masks_th.astype(bool), min_size=50).astype(int)

    # Segment nuclei
    masks_nuclei, _, _ = model_nuclei.eval(nuclei_img_processed, diameter=None)

    # %% Create combined cell mask and filter nuclei
    combined_cell_mask = (masks_map2 > 0) | (masks_th > 0)

    # Keep only nuclei whose centroid falls inside the combined cell mask and area >= min_area
    labeled_nuclei = label(masks_nuclei > 0)
    regions = regionprops(labeled_nuclei, intensity_image=nuclei_img)
    valid_nuclei = []
    min_area = 50
    for r in regions:
        cy, cx = [int(c) for c in r.centroid]
        if combined_cell_mask[cy, cx] and r.area >= min_area:
            valid_nuclei.append(r)
    
    # Merge nuclei closer than 100 pixels
    merged_labels = np.zeros_like(labeled_nuclei)
    used = set()
    merged_nuclei = []
    distance_thresh = 300

    for i, r1 in enumerate(valid_nuclei):
        if r1.label in used:
            continue
        cy1, cx1 = r1.centroid
        group = [r1]
        used.add(r1.label)
        for j, r2 in enumerate(valid_nuclei):
            if r2.label in used:
                continue
            cy2, cx2 = r2.centroid
            dist = np.sqrt((cy1 - cy2)**2 + (cx1 - cx2)**2)
            if dist < distance_thresh:
                group.append(r2)
                used.add(r2.label)
        # Create a merged region mask
        merged_mask = np.isin(labeled_nuclei, [r.label for r in group])
        merged_label = len(merged_nuclei) + 1
        merged_labels[merged_mask] = merged_label
        merged_nuclei.append((merged_label, merged_mask))

    total_nuclei = len(merged_nuclei)

    # Recompute TH positive using merged nuclei
    th_positive = 0
    for label_id, nucleus_mask in merged_nuclei:
        overlap_ratio = np.sum((masks_th > 0) & nucleus_mask) / np.sum(nucleus_mask)
        if overlap_ratio > 0.05:
            th_positive += 1

    percent_th_positive = 100 * th_positive / total_nuclei if total_nuclei > 0 else 0

    # %% Visualization overlays
    base_image_norm = exposure.rescale_intensity(map2_img, out_range=(0, 1))
    overlay_map2 = label2rgb(masks_map2, image=base_image_norm, bg_label=0, alpha=0.3)
    overlay_th = label2rgb(masks_th, image=base_image_norm, bg_label=0, alpha=0.3)
    overlay_cellmask = label2rgb(combined_cell_mask.astype(int), image=base_image_norm, bg_label=0, alpha=0.3)

    # Build nuclei overlay only for valid nuclei
    filtered_nuclei_mask = np.zeros_like(labeled_nuclei)
    for r in valid_nuclei:
        filtered_nuclei_mask[labeled_nuclei == r.label] = r.label
    overlay_nuclei = label2rgb(filtered_nuclei_mask, image=base_image_norm, bg_label=0, alpha=0.3)

    # Plot and save
    fig, axs = plt.subplots(1, 4, figsize=(20, 5))
    axs[0].imshow(overlay_map2)
    axs[0].set_title("MAP2 mask")
    axs[1].imshow(overlay_th)
    axs[1].set_title("TH mask")
    axs[2].imshow(overlay_cellmask)
    axs[2].set_title("Cell mask")
    axs[3].imshow(overlay_nuclei)
    axs[3].set_title("Filtered nuclei")
    for ax in axs:
        ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_dir / f"{group_name}_MAP2_TH_CellMask_Nuclei.png", dpi=200)
    plt.close()

    summary_data.append({
        "file": file_info['th_file'].name,
        "condition": condition,
        "total_nuclei": total_nuclei,
        "th_positive": th_positive,
        "percent_th_positive": percent_th_positive
    })

    # Save summary for this file
    with open(output_dir / f"{group_name}_{condition}_summary.txt", "w") as f:
        f.write(f"Condition: {condition}\n")
        f.write(f"File: {file_info['th_file'].name}\n")
        f.write(f"Total nuclei: {total_nuclei}\n")
        f.write(f"TH+ nuclei: {th_positive}\n")
        f.write(f"Percent TH+: {percent_th_positive:.2f}\n")

# %% Save summary CSV
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(output_dir / "TH_quantification_summary.csv", index=False)


# %% Plot results
# Plotting
# Add columns for genotype and neuron type
summary_df[["genotype", "neurontype"]] = summary_df["condition"].str.extract(r"(\w+)_([a-zA-Z]+)")

# Set colors
palette = {
    ("Ctrl", "iN"): "lightgrey",
    ("Ctrl", "iDA"): "dimgrey",
    ("ASAH1", "iN"): "lightblue",
    ("ASAH1", "iDA"): "blue",
    ("SMPD1", "iN"): "lightskyblue",
    ("SMPD1", "iDA"): "steelblue",
}

# Plot
plt.figure(figsize=(8, 6))
sns.barplot(data=summary_df, x="genotype", y="percent_th", hue="neurontype",
            palette=[palette[(g, n)] for g, n in summary_df[["genotype", "neurontype"]].drop_duplicates().values],
            ci="sd", capsize=0.1)
sns.stripplot(data=summary_df, x="genotype", y="percent_th", hue="neurontype",
              dodge=True, palette=[palette[(g, n)] for g, n in summary_df[["genotype", "neurontype"]].drop_duplicates().values],
              alpha=0.6, jitter=True, linewidth=0.5, edgecolor="black")

plt.ylabel("% TH+ cells")
plt.xlabel("Genotype")
plt.legend(title="Neuron Type", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(output_dir / "TH_quantification_plot.pdf")
plt.close()