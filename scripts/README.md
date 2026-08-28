 Scripts

This directory contains the scripts used for sequence processing, metadata preparation, phylogenetic analysis, genotype assignment, and figure generation for the Human Parechovirus (HPeV) analyses.

## Genotype assignment

`relabel_tips.py` assigns genotypes based on phylogenetic proximity. For each sequence with an unknown genotype, the script identifies the *k* nearest sequences with known genotypes within a defined distance cutoff and assigns the genotype supported by the majority of these neighbors.
