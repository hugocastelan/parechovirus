 # Scripts

This directory contains the scripts used for sequence processing, metadata preparation, phylogenetic analysis, genotype assignment, and figure generation for the Human Parechovirus (HPeV) analyses.

### Description of scripts 

`relabel_tips.py` assigns genotypes based on phylogenetic proximity. For each sequence with an unknown genotype, the script identifies the *k* nearest sequences with known genotypes within a defined distance cutoff and assigns the genotype supported by the majority of these neighbors.

`Figure_1_map_and_bubble_plot.R` This script generates a global distribution map and a temporal bubble plot for HPeV genotypes using metadata. It includes data cleaning, visualization functions, and exports the final figure as a PDF.

`figure2_Tree_panelA_panelB_pevA_VP1_SRA.R`  Loads the PeV-A VP1 phylogenetic tree and its metadata, applies midpoint rooting, and matches each tip to its genotype, country, year, and sequence type (SRA vs non-SRA). Builds Panel A with the complete tree, black branches, and all tips colored by genotype, using a legend that includes the total number of sequences per genotype. Builds Panel B with the same tree but showing only SRA tips as colored points. 

`figure_3_genetics_analysis.R` The script analyzes the genetic diversity of HPeV using VP1 sequences. It compares how different the genotypes are using genetic distances, PCoA, and PERMANOVA; shows which genotypes are most similar using a dendrogram; and calculates the diversity within each genotype using haplotypes, nucleotide diversity (π), and Tajima’s D.

`subsampling_genotype_3.R` downsamples a FASTA file of HPeV3 sequences by selecting a maximum of 10 sequences per geographic region and year combination.


