#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO


# Input files


CSV_IN = "/Users/hugo/Desktop/vp1_parachovirus/dates_treetime_v11.csv"
FASTA_IN = "/Users/hugo/Desktop/vp1_parachovirus/vp1_headers_updated_remove_v11.fasta"


# Output files


CSV_OUT = "/Users/hugo/Desktop/vp1_parachovirus/dates_treetime_v11_subsampled.csv"
FASTA_OUT = "/Users/hugo/Desktop/vp1_parachovirus/vp1_headers_updated_remove_v11_subsampled.fasta"


# Parameters


MAX_SEQUENCES_PER_YEAR = 10
RANDOM_SEED = 42


# Load metadata

df = pd.read_csv(CSV_IN)

# Convert decimal dates to integer years
df["year"] = df["date"].astype(float).astype(int)


# Subsample sequences by year

print("Performing temporal subsampling")

subsampled_groups = []

for year, group in df.groupby("year"):

    n_sequences = len(group)

    if n_sequences <= MAX_SEQUENCES_PER_YEAR:

        # Keep all sequences for sparsely sampled years
        selected = group.copy()

    else:

        # Randomly select a fixed number of sequences
        selected = group.sample(
            n=MAX_SEQUENCES_PER_YEAR,
            random_state=RANDOM_SEED
        )

    subsampled_groups.append(selected)

subsampled_df = pd.concat(subsampled_groups)

# Sort chronologically
subsampled_df = subsampled_df.sort_values("year")


# Save subsampled metadata


subsampled_df.to_csv(CSV_OUT, index=False)

print(f"\nSubsampled metadata saved to:")
print(CSV_OUT)


# Extract selected sequence names


selected_names = set(subsampled_df["name"])

print(f"\nSelected sequences: {len(selected_names)}")


# Filter FASTA file

print("\nFiltering FASTA file...")

selected_records = []

for record in SeqIO.parse(FASTA_IN, "fasta"):

    if record.id in selected_names:
        selected_records.append(record)

SeqIO.write(selected_records, FASTA_OUT, "fasta")


# Summary statistics


print("\nSubsampling summary:")
print("-" * 40)

summary = (
    subsampled_df["year"]
    .value_counts()
    .sort_index()
)

print(summary)

print("\nTotal metadata records retained:")
print(len(subsampled_df))

print("\nTotal FASTA sequences retained:")
print(len(selected_records))

print("\nOutput FASTA:")
print(FASTA_OUT)

print("\nDone.")
