from cellpose import models, core, io, plot
import nd2
import matplotlib.pyplot as plt
import tifffile
import numpy as np
from pathlib import Path
from tqdm import tqdm
import pandas as pd
from skimage.measure import regionprops_table, label as sk_label

# This script processes a folder of .nd2 files by applying CellposeSAM segmentation,
# saves the resulting label masks, and computes average intensity over time per label.

io.logger_setup() # run this to get printing of progress


# Set up input and output directories for ND2 files and label outputs
#folder = Path("nd2")
folder = Path("/Volumes/Expansion/Harvard/NIC/Station8/Project_LSD-PD/20241202_diff129_iDA-d33_Ctrl-ASAH1e2_Calcium/nd2")
output = Path("nd2_labels")
if not folder.exists():
    raise FileNotFoundError(f"Input folder {folder} does not exist.")
if not output.exists():
    output.mkdir(parents=True)
# Optionally, check if there are any .nd2 files to process
if not any(folder.glob('*.nd2')):
    print(f"Warning: No .nd2 files found in {folder}")


# models are stored under
# /Users/felix/.cellpose/models/... 

use_gpu = core.use_gpu()
#model = models.CellposeModel(gpu=use_gpu)
model = models.CellposeModel(gpu=True,
                             pretrained_model="/Users/felix/.cellpose/models/cpsam_20250530_neuronSoma_v2")

nd2_files = sorted([f for f in folder.glob("*.nd2") if not f.name.startswith("._")])
for nd2_file in tqdm(nd2_files, desc="Running Cellpose"):
    print(nd2_file.name)
    t_stack = nd2.imread(nd2_file)

    image = np.max(t_stack, axis=0)

    flow_threshold = 0.4
    cellprob_threshold = 0.0
    tile_norm_blocksize = 0
    batch_size=8

    masks, flows, styles = model.eval(image, batch_size=batch_size, flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold,
                                      normalize={"tile_norm_blocksize": tile_norm_blocksize})

    fig = plt.figure(figsize=(12,5))
    plot.show_segmentation(fig, image, masks, flows[0])
    plt.tight_layout()
    plt.savefig(output / f"{nd2_file.name}_labels.png", dpi=600)
    # plt.show()
    plt.close()

    tifffile.imwrite(output / f"{nd2_file.name}_labels.tif", masks, imagej=True)

    print(masks.max())

    # Define thresholds for valid areas
    MIN_AREA = 100
    MAX_AREA = 10000
    
    # Compute region properties and filter by area
    props = regionprops_table(sk_label(masks), properties=("label", "area"))
    valid_labels = [label for label, area in zip(props["label"], props["area"]) if MIN_AREA < area < MAX_AREA]
    valid_labels = [label for label, area in zip(props["label"], props["area"]) if 100 < area < 10000]

    # Compute centroids for each valid label and save to CSV
    centroids = regionprops_table(sk_label(masks), properties=("label", "centroid"))
    centroid_data = [
        {
            "file": nd2_file.stem,
            "label": label,
            "centroid_x": x,
            "centroid_y": y
        }
        for label, y, x in zip(centroids["label"], centroids["centroid-0"], centroids["centroid-1"])
        if label in valid_labels
    ]
    centroid_df = pd.DataFrame(centroid_data)
    centroid_df.to_csv(output / f"{nd2_file.stem}_centroids.csv", index=False)

    intensity_data = []

    for label in valid_labels:
        label_mask = masks == label
        indices = np.where(label_mask)  # Get indices of the labeled region
        mean_intensity = t_stack[:, indices[0], indices[1]].mean(axis=1)  # Compute mean intensity over time
        for t, intensity in enumerate(mean_intensity):
            intensity_data.append({
                "file": nd2_file.stem,
                "label": label,
                "time": t,
                "intensity": intensity
            })
    # save centroids only of valid labels, meaning that within the range of area
    valid_centroid_df = pd.DataFrame(centroid_data)
    valid_centroid_df.to_csv(output / f"{nd2_file.stem}_centroids_valid.csv", index=False)

    df = pd.DataFrame(intensity_data)
    # Optionally include filename as a column in the pivot table
    include_filename = True  # Set to False to exclude filename column
    # Pivot the DataFrame to have time as rows and label as columns, with intensity as values
    pivot_df = df.pivot(index="time", columns="label", values="intensity")
    if include_filename:
        pivot_df.insert(0, "file", nd2_file.stem)  # optional: include filename as a column
    pivot_df.to_csv(output / f"{nd2_file.stem}_intensity.csv")