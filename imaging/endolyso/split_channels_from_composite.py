import re
from pathlib import Path
import tifffile as tiff
import numpy as np

# Directory containing composite TIFFs
input_dir = Path("/Volumes/Expansion/Harvard/NIC/LucilleTwo/Project_neuroLSD/20250613_HeLa-Ctrl-ASAH1_EndoLysoStains/20250613_HeLa-Ctrl-ASAH1_EndoLysoStains_MIPs")
output_dir = Path("/Volumes/Expansion/Harvard/NIC/LucilleTwo/Project_neuroLSD/20250613_HeLa-Ctrl-ASAH1_EndoLysoStains/split_channels_renamed")
output_dir.mkdir(parents=True, exist_ok=True)

# Regex to extract channel-stain mapping from filename
channel_pat = re.compile(r'([A-Za-z0-9]+)-(\d{3,4})')

# Define the order of channels in the TIFF (1-based index)
channel_order = [568, 488, 405, 647]

# DEBUG: Print found files
tiff_files = list(input_dir.glob("*.tif"))
print(f"Found {len(tiff_files)} TIFF files in {input_dir}")
if not tiff_files:
    print("No TIFF files found. Check your input directory path and file extensions.")

for tiff_path in tiff_files:
    print(f"Processing: {tiff_path.name}")
    matches = channel_pat.findall(tiff_path.name)
    mapping = {}
    for stain, wl in matches:
        mapping[int(wl)] = stain
    mapping[405] = "DNA"  # Always assign 405nm to DNA

    try:
        img = tiff.imread(tiff_path)
    except Exception as e:
        print(f"Could not read {tiff_path}: {e}")
        continue

    # Print shape for debugging
    print(f"Image shape: {img.shape}")

    if img.ndim != 3 or img.shape[0] != 4:
        print(f"Skipping {tiff_path.name}: expected 4-channel image, got shape {img.shape}")
        continue

    for idx, wl in enumerate(channel_order):
        stain = mapping.get(wl, f"ch{idx+1}")
        out_name = f"{tiff_path.stem}_{stain}-{wl}.tif"
        tiff.imwrite(output_dir / out_name, img[idx])
        print(f"Saved: {output_dir / out_name}")
