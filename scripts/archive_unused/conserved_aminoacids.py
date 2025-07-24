import logomaker
import pandas as pd
from Bio import AlignIO
import matplotlib.pyplot as plt
import os

# ---------------- CONFIG ----------------

alignment_files = snakemake.input  

min_region_length = 10  # Minimum continuous region length to plot
max_allowed_gaps = 2    # Maximum number of allowed gaps per column

# ---------------- PROCESS EACH ALIGNMENT FILE ----------------
for alignment_file in alignment_files:
    print(f"\nProcessing alignment: {alignment_file}")

    base_folder = os.path.dirname(alignment_file)
    filename = os.path.basename(alignment_file)
    gene_name = filename.replace("_aligned.fasta", "")

    output_folder = os.path.join(base_folder, f"{gene_name}_sequence_logos_flexible")
    os.makedirs(output_folder, exist_ok=True)

    output_table_file = os.path.join(output_folder, f"{gene_name}_aa_frequencies.tsv")
    output_regions_file = os.path.join(output_folder, f"{gene_name}_detected_regions.txt")

    # ---------------- LOAD MSA ----------------
    alignment = AlignIO.read(alignment_file, "fasta")
    aligned_seqs = [str(record.seq) for record in alignment]
    n_sequences = len(aligned_seqs)

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

        df_info_zoom = logomaker.transform_matrix(df_freq_zoom, from_type='probability', to_type='information')

        fig, ax = plt.subplots(figsize=(max(10, (region_end - region_start) * 0.5), 4))

        logo = logomaker.Logo(df_info_zoom,
                              ax=ax,
                              font_name='Arial',
                              color_scheme='chemistry',
                              show_spines=False)

        logo.style_spines(visible=False)
        logo.ax.set_ylabel("Information (bits)", fontsize=12)
        logo.ax.set_xlabel(f"Position {region_start}-{region_end-1}", fontsize=12)
        logo.ax.set_xticks(range(region_start, region_end, max(1, (region_end - region_start) // 10)))
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
    print(f"Total positions kept after gap filtering: {len(positions_to_keep)}")
    print(f"Regions detected: {len(regions)}\n")
