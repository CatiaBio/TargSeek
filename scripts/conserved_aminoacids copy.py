import logomaker
import pandas as pd
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
import os
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, leaves_list

# ---------------- CONFIG ----------------

alignment_file = "data/sequence_alg.fasta"
#alignment_file = "results/gram_negative/bamA_alg.fasta"

min_region_length = 10  # Minimum continuous region length to plot
max_allowed_gaps = 0    # Maximum number of allowed gaps per column
min_stack_height = 0.5  # Minimum information (bits) per column to plot
logo_type = "information"  # Options: "probability", "information", "log_odds"
background_model = "BLOSUM62"  # Options: "uniform", "BLOSUM62"
start_position_offset = 1  # Starting position (e.g., 1 for mature proteins)

# Clustering (optional)
clustering_method = None  # Options: None, "hierarchical"

# ---------------- BACKGROUND FREQUENCIES ----------------

blosum62_background = {
    'A': 0.078, 'R': 0.051, 'N': 0.042, 'D': 0.053, 'C': 0.019,
    'Q': 0.040, 'E': 0.063, 'G': 0.074, 'H': 0.026, 'I': 0.068,
    'L': 0.099, 'K': 0.058, 'M': 0.025, 'F': 0.047, 'P': 0.039,
    'S': 0.068, 'T': 0.058, 'W': 0.013, 'Y': 0.032, 'V': 0.073
}

# ---------------- LOAD MSA ----------------

base_folder = os.path.dirname(alignment_file)
filename = os.path.basename(alignment_file)
gene_name = filename.replace("_aligned.fasta", "").replace("_alg.fasta", "")

output_folder = os.path.join(base_folder, f"{gene_name}_sequence_logos_flexible")
os.makedirs(output_folder, exist_ok=True)

output_table_file = os.path.join(output_folder, f"{gene_name}_aa_frequencies.tsv")
output_regions_file = os.path.join(output_folder, f"{gene_name}_detected_regions.txt")

alignment = AlignIO.read(alignment_file, "fasta")
aligned_seqs = [str(record.seq) for record in alignment]
n_sequences = len(aligned_seqs)

real_lengths = [len(seq.replace("-", "")) for seq in aligned_seqs]
avg_real_length = sum(real_lengths) / len(real_lengths)
print(f"Average ungapped sequence length: {avg_real_length:.2f} residues")

# ---------------- BUILD COUNT MATRIX ----------------

counts = {}
for seq in aligned_seqs:
    for i, aa in enumerate(seq):
        if i not in counts:
            counts[i] = {}
        counts[i][aa] = counts[i].get(aa, 0) + 1

df_counts = pd.DataFrame.from_dict(counts, orient="index").fillna(0)

# Keep only standard amino acids
standard_aas = list("ACDEFGHIKLMNPQRSTVWY")
df_counts = df_counts[[aa for aa in standard_aas if aa in df_counts.columns]]


# ---------------- FLEXIBLE GAP FILTER ----------------
gaps_per_position = n_sequences - df_counts.sum(axis=1)
positions_to_keep = gaps_per_position[gaps_per_position <= max_allowed_gaps].index

df_counts_filtered = df_counts.loc[positions_to_keep]
df_counts_filtered.index = range(1, len(df_counts_filtered) + 1)  # <-- 1-based numbering

# ---------------- CALCULATE FREQUENCIES ----------------
df_freq = df_counts_filtered.div(n_sequences)
df_freq = df_freq.round(4)

# Save frequency table
df_freq.to_csv(output_table_file, sep="\t")

# ---------------- FIND REGIONS ----------------

positions = sorted(df_freq.index)
regions = []

start = None
for i in range(len(positions)):
    if start is None:
        start = positions[i]

    if i == len(positions) - 1 or positions[i+1] != positions[i] + 1:
        end = positions[i] + 1  # exclusive
        if end - start >= min_region_length:
            regions.append((start, end))
        start = None

# Save regions summary
with open(output_regions_file, "w") as f:
    if regions:
        for idx, (region_start, region_end) in enumerate(regions, 1):
            f.write(f"Region {idx}: Positions {region_start}–{region_end - 1}\n")
    else:
        f.write("No regions found with the given thresholds.\n")

# ---------------- PLOT EACH REGION ----------------

for idx, (region_start, region_end) in enumerate(regions):
    df_freq_zoom = df_freq.loc[(df_freq.index >= region_start) & (df_freq.index < region_end)]

    # Apply clustering if selected
    if clustering_method == "hierarchical":
        dists = pdist(df_freq_zoom.fillna(0).T, metric='cosine')
        linkage_matrix = linkage(dists, method='average')
        order = leaves_list(linkage_matrix)
        df_freq_zoom = df_freq_zoom.iloc[:, order]
    
    # Build background frequencies matching df_freq_zoom shape
    if background_model == "BLOSUM62":
        background_array = np.array([blosum62_background[aa] for aa in df_freq_zoom.columns])
        background_df = pd.DataFrame(
            np.tile(background_array, (df_freq_zoom.shape[0], 1)), 
            index=df_freq_zoom.index, 
            columns=df_freq_zoom.columns
        )
    else:
        background_df = None

    # Transform matrix
    if logo_type == "information":
        df_logo = logomaker.transform_matrix(
            df_freq_zoom,
            from_type='probability',
            to_type='information',
            background=background_df
        )
    elif logo_type == "log_odds":
        df_logo = logomaker.transform_matrix(
            df_freq_zoom,
            from_type='probability',
            to_type='log_odds',
            background=background_df
        )
    else:
        df_logo = df_freq_zoom.copy()

    
    # Filter by minimum stack height
    stack_heights = df_logo.sum(axis=1)
    positions_to_plot = stack_heights[stack_heights >= min_stack_height].index
    df_logo_filtered = df_logo.loc[positions_to_plot]

    # Skip empty plots
    if df_logo_filtered.empty:
        continue

    # Plot
    fig, ax = plt.subplots(figsize=(max(10, (region_end - region_start) * 0.5), 4))

    logo = logomaker.Logo(
        df_logo_filtered,
        ax=ax,
        font_name='Arial',
        color_scheme='chemistry',
        show_spines=False
    )

    logo.style_spines(visible=False)
    logo.ax.set_ylabel("Information (bits)", fontsize=12)
    logo.ax.set_xlabel(f"Position {region_start + start_position_offset}-{region_end - 1 + start_position_offset}", fontsize=12)
    logo.ax.set_xticks(
        range(region_start + start_position_offset, region_end + start_position_offset),
        minor=False
    )
    logo.ax.tick_params(labelsize=10)
    logo.ax.set_ylim(0, 5)

    plt.tight_layout()
    output_logo_file = os.path.join(output_folder, f"{gene_name}_sequence_logo_region{idx+1}.png")
    plt.savefig(output_logo_file, dpi=300)
    plt.close()

    print(f"Saved logo for {gene_name} region {region_start}-{region_end-1}")

if not regions:
    print(f"No good regions found for {gene_name}")

# ---------------- SUMMARY ----------------
print(f"\nSUMMARY for {gene_name}:")
print(f"Total sequences: {n_sequences}")
print(f"Average ungapped sequence length: {avg_real_length:.2f}")
print(f"Total positions kept after gap filtering: {len(positions_to_keep)}")
print(f"Regions detected: {len(regions)}\n")