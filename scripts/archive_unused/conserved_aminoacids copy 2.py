import logomaker
import pandas as pd
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
import os

# ---------------- CONFIG ----------------

alignment_file = "results/gram_negative/bamA_alg_mod.fasta"
output_folder = "results/gram_negative/bamA_sequence_logos_complete"
os.makedirs(output_folder, exist_ok=True)

chunk_size = 30  # number of amino acids per plot
max_allowed_gaps = 5
min_ic_per_position = 0.3  # bits threshold
min_aa_prob_to_plot = 0.1

# ---------------- LOAD MSA ----------------
alignment = AlignIO.read(alignment_file, "fasta")
aligned_seqs = [str(record.seq) for record in alignment]
n_sequences = len(aligned_seqs)
sequence_length = alignment.get_alignment_length()

# ---------------- BUILD FREQUENCY MATRIX ----------------
counts = {}
for seq in aligned_seqs:
    for i, aa in enumerate(seq):
        if i not in counts:
            counts[i] = {}
        counts[i][aa] = counts[i].get(aa, 0) + 1

df_counts = pd.DataFrame.from_dict(counts, orient="index").fillna(0)
standard_aas = list("ACDEFGHIKLMNPQRSTVWY")
df_counts = df_counts[[aa for aa in standard_aas if aa in df_counts.columns]]

# ---------------- GAP FILTER ----------------
gaps_per_position = n_sequences - df_counts.sum(axis=1)
positions_to_keep = gaps_per_position[gaps_per_position <= max_allowed_gaps].index
df_counts_filtered = df_counts.loc[positions_to_keep]

# ---------------- FREQUENCIES ----------------
df_freq = df_counts_filtered.div(n_sequences)

# ---------------- INFORMATION CONTENT ----------------
def compute_information_content(freq_df):
    ic_values = []
    for _, row in freq_df.iterrows():
        probs = row[row > 0]
        entropy = -np.sum(probs * np.log2(probs))
        ic = np.log2(20) - entropy
        ic_values.append(ic)
    return pd.Series(ic_values, index=freq_df.index)

df_ic = compute_information_content(df_freq)

# ---------------- PLOT ENTIRE SEQUENCE IN CHUNKS ----------------
for chunk_start in range(0, sequence_length, chunk_size):
    chunk_end = min(chunk_start + chunk_size, sequence_length)

    df_freq_chunk = df_freq.loc[(df_freq.index >= chunk_start) & (df_freq.index < chunk_end)]
    if df_freq_chunk.empty:
        continue

    df_info_chunk = logomaker.transform_matrix(df_freq_chunk, from_type='probability', to_type='information')

    ic_per_position = df_info_chunk.sum(axis=1)
    df_info_chunk_filtered = df_info_chunk.loc[ic_per_position >= min_ic_per_position]

    if df_info_chunk_filtered.empty:
        print(f"Skipping low-information chunk: {chunk_start}-{chunk_end-1}")
        continue

    df_info_chunk_filtered = df_info_chunk_filtered.where(df_freq_chunk >= min_aa_prob_to_plot, other=0)

    fig_width = max(10, len(df_info_chunk_filtered) * 0.5)
    fig, ax = plt.subplots(figsize=(fig_width, 4))

    logo = logomaker.Logo(df_info_chunk_filtered,
                          ax=ax,
                          font_name='Arial',
                          color_scheme='chemistry')

    ax.set_ylabel("Information (bits)", fontsize=14)
    ax.set_xlabel(f"Positions {df_info_chunk_filtered.index[0]}–{df_info_chunk_filtered.index[-1]}", fontsize=14)
    ax.tick_params(labelsize=12)
    ax.set_ylim(0, 5)

    tick_interval = max(1, len(df_info_chunk_filtered) // 10)
    ax.set_xticks(df_info_chunk_filtered.index[::tick_interval])

    plt.tight_layout()

    output_logo_file = os.path.join(output_folder, f"sequence_logo_{df_info_chunk_filtered.index[0]}-{df_info_chunk_filtered.index[-1]}.png")
    plt.savefig(output_logo_file, dpi=300)
    plt.close()

    print(f"Saved complete-sequence logo: {df_info_chunk_filtered.index[0]}–{df_info_chunk_filtered.index[-1]}")