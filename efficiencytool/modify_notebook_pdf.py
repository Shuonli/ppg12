#!/usr/bin/env python3
import json

# Read the notebook
with open('eff_plot.ipynb', 'r') as f:
    nb = json.load(f)

# Find the Truth Purity Plots cell (it's the 4th cell, index 3)
cell = nb['cells'][3]

# Add output directory creation after model_name definition
source = cell['source']

# Find the line with root_file_path and add output directory creation after it
for i, line in enumerate(source):
    if 'root_file_path = f' in line:
        # Insert output directory creation after this line
        source.insert(i+2, '    "# Create output directory with model name prefix\\n",')
        source.insert(i+3, '    "output_dir = f\\"plot_eff_{model_name}\\"\\n",')
        source.insert(i+4, '    "if not os.path.exists(output_dir):\\n",')
        source.insert(i+5, '    "    os.makedirs(output_dir)\\n",')
        source.insert(i+6, '    "    print(f\\"Created output directory: {output_dir}\\")\\n",')
        source.insert(i+7, '    "\\n",')
        break

# Update PNG save paths to PDF and add output directory
for i, line in enumerate(source):
    if 'png_name = f"truth_purity_{i}.png"' in line:
        source[i] = '    "            png_name = f\\"{output_dir}/truth_purity_{i}.pdf\\"\\n",'
    elif 'png_name_total = f"total_tight_iso_clusters_{i}.png"' in line:
        source[i] = '    "            png_name_total = f\\"{output_dir}/total_tight_iso_clusters_{i}.pdf\\"\\n",'
    elif 'png_name_rel_err = f"total_tight_iso_clusters_rel_err_{i}.png"' in line:
        source[i] = '    "            png_name_rel_err = f\\"{output_dir}/total_tight_iso_clusters_rel_err_{i}.pdf\\"\\n",'
    elif 'comb_png = "truth_purity_all.png"' in line:
        source[i] = '    "            comb_png = f\\"{output_dir}/truth_purity_all.pdf\\"\\n",'
    elif 'comb_total_png = "total_tight_iso_clusters_all.png"' in line:
        source[i] = '    "            comb_total_png = f\\"{output_dir}/total_tight_iso_clusters_all.pdf\\"\\n",'
    elif 'comb_rel_err_png = "total_tight_iso_clusters_rel_err_all.png"' in line:
        source[i] = '    "            comb_rel_err_png = f\\"{output_dir}/total_tight_iso_clusters_rel_err_all.pdf\\"\\n",'
    elif 'output_root_path = "vertex_truth_purity_output.root"' in line:
        source[i] = '    "                output_root_path = f\\"{output_dir}/vertex_truth_purity_output.root\\"\\n",'

# Write the notebook back
with open('eff_plot.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

print('Successfully updated the Truth Purity Plots cell with output directory and PDF functionality')


