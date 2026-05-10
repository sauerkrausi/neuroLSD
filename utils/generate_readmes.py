import csv
import os
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(SCRIPT_DIR, "raw_readme.csv")

metadata = {}
panels = []

with open(CSV_FILE, newline="", encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    rows = list(reader)

# Metadata rows until blank line
i = 0
while i < len(rows):
    row = rows[i]
    if all(cell.strip() == "" for cell in row):
        i += 1
        break
    if row[0].strip():
        metadata[row[0].strip()] = row[1].strip()
    i += 1

# Header row
header = [col.strip() for col in rows[i]]
i += 1

# Panel rows
for row in rows[i:]:
    if all(cell.strip() == "" for cell in row):
        continue
    entry = {header[j]: row[j].strip() if j < len(row) else "" for j in range(len(header))}
    panels.append(entry)

# Group by methodology_collection
groups = defaultdict(list)
for panel in panels:
    key = panel.get("methodology_collection", "").strip()
    if key:
        groups[key].append(panel)

def val(p, key):
    v = p.get(key, "").strip()
    return v if v and v != "N/A" else ""

def make_readme(method_name, method_panels, metadata, out_dir):
    filename = os.path.join(out_dir, f"_README_{method_name.replace(' ', '_')}.md")
    lines = []

    lines += [
        "## GENERAL INFORMATION",
        "",
        f"**Title:** {metadata.get('title', '')}",
        "",
        f"**Citation:** {metadata.get('bioRvix', '') or metadata.get('citation', '')}",
        "",
        f"**Principal Investigator:** {metadata.get('principal_investigator', '')}",
        "",
        f"**Authors:** {metadata.get('authors', '')}",
        "",
        f"**License:** {metadata.get('license', '')}",
        "",
        "---",
        "",
        "## FIGURE PANEL OVERVIEW",
        "",
    ]

    for p in method_panels:
        line = f"- **{val(p, 'figure_panel')}** | {val(p, 'method_name')}"
        if val(p, 'instrument_info'):
            line += f" | {val(p, 'instrument_info')}"
        lines.append(line)

    lines += ["", "---", "", "## DATA SPECIFIC INFORMATION", ""]

    for p in method_panels:
        lines.append(f"### Panel {val(p, 'figure_panel')}")
        lines.append("")
        lines.append(f"**Method:** {val(p, 'method_name')}")
        if val(p, 'instrument_info'):
            lines.append(f"**Instrument:** {val(p, 'instrument_info')}")
        if val(p, 'code_file_1'):
            lines.append(f"**Code File:** {val(p, 'code_file_1')}")
        if val(p, 'code_file_2'):
            lines.append(f"**Code File:** {val(p, 'code_file_2')}")
        if val(p, 'code_module'):
            lines.append(f"**Code Module:** {val(p, 'code_module')}")
        if val(p, 'code_info'):
            lines.append(f"**Code Info:** {val(p, 'code_info')}")
        if val(p, 'code_link'):
            lines.append(f"**Code Repository:** {val(p, 'code_link')}")
        if val(p, 'dependencies'):
            lines.append(f"**Dependencies:** {val(p, 'dependencies')}")
        if val(p, 'dataset'):
            lines.append(f"**Dataset ID:** {val(p, 'dataset')}")
        if val(p, 'dataset_info'):
            lines.append(f"**Dataset Description:** {val(p, 'dataset_info')}")
        if val(p, 'data_link'):
            lines.append(f"**Data Link:** {val(p, 'data_link')}")
        lines.append("")

    with open(filename, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"Written: {filename}")

def make_master_readme(groups, metadata, out_dir):
    filename = os.path.join(out_dir, "_README_index.md")
    lines = []

    biorxiv = metadata.get('bioRvix', '')
    zenodo  = metadata.get('zenodo', '')

    lines += [
        "## GENERAL INFORMATION",
        "",
        f"**Title:** {metadata.get('title', '')}",
        "",
        f"**Citation:** {biorxiv if biorxiv else metadata.get('citation', '')}",
        "",
        f"**Principal Investigator:** {metadata.get('principal_investigator', '')}",
        "",
        f"**Authors:** {metadata.get('authors', '')}",
        "",
        f"**License:** {metadata.get('license', '')}",
        "",
        "---",
        "",
        "## DATA AVAILABILITY",
        "",
        "- Proteomic data (.RAW files) for nDIA of iN/iDA neurons: MASSIVE MSV000099237 (https://massive.ucsd.edu)",
        "- HeLa LysoIP TMT data: ProteomeXchange PXD067219",
        "- Lipidomics data: Metabolomics Workbench ST004217 (http://dx.doi.org/10.21228/M8Z556)",
        "- CryoET tomograms: EMDB EMD-55210 (Control), EMD-55211 (ASAH1-/-)",
        f"- Key Resource Table, source data, segmentation models: Zenodo {zenodo}" if zenodo else "- Key Resource Table, source data, segmentation models: Zenodo 10.5281/zenodo.16733440",
        "- Analysis scripts: GitHub https://github.com/sauerkrausi/neuroLSD",
        f"- Preprint: {biorxiv}" if biorxiv else "",
        "",
        "---",
        "",
        "## README INDEX",
        "",
    ]

    for method_name in sorted(groups.keys()):
        fname = f"_README_{method_name.replace(' ', '_')}.md"
        panels = [p.get("figure_panel", "") for p in groups[method_name]]
        lines.append(f"### [{method_name}]({fname})")
        lines.append(f"Panels: {', '.join(panels)}")
        lines.append("")

    with open(filename, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"Written: {filename}")

out_dir = SCRIPT_DIR
for method_name, method_panels in groups.items():
    make_readme(method_name, method_panels, metadata, out_dir)

make_master_readme(groups, metadata, out_dir)
print("Done.")
