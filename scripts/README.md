 # Scripts

This directory contains the scripts used for sequence processing, metadata preparation, phylogenetic analysis, genotype assignment, and figure generation for the Human Parechovirus (HPeV) analyses.

### Description of scripts 

`relabel_tips.py` assigns genotypes based on phylogenetic proximity. For each sequence with an unknown genotype, the script identifies the *k* nearest sequences with known genotypes within a defined distance cutoff and assigns the genotype supported by the majority of these neighbors.

`figure_2_genetics_analysis.R` The script analyzes the genetic diversity of HPeV using VP1 sequences. It compares how different the genotypes are using genetic distances, PCoA, and PERMANOVA; shows which genotypes are most similar using a dendrogram; and calculates the diversity within each genotype using haplotypes, nucleotide diversity (π), and Tajima’s D.



`figure_2_genetics_analysis.R` downsamples a FASTA file of HPeV3 sequences by selecting a maximum of 10 sequences per geographic region and year combination.
