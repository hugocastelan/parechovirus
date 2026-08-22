###############################################################
# HPeV VP1 PHYLOGENETIC TREE
###############################################################
rm(list = ls())

library(ape)
library(phytools)
library(ggtree)
library(treeio)
library(tidyverse)

setwd("/Users/hugo/Desktop/vp1_parachovirus")

tree <- read.tree("vp1_headers_updated_v2.nwk")

cat("Number of tips:", Ntip(tree), "\n")
cat("Number of internal nodes:", tree$Nnode, "\n")

tree_rooted <- midpoint.root(tree)
tree_rooted <- ladderize(tree_rooted, right = FALSE)

meta <- tibble(label = tree_rooted$tip.label) %>%
  separate(label, into = c("Accession", "Virus", "Genotype", "Country", "Year"), sep = "\\|", remove = FALSE, fill = "right", extra = "merge") %>%
  mutate(Country = str_replace_all(Country, "_", " "), Year = suppressWarnings(as.numeric(Year)))

ordered_genotypes <- paste0("HPeV", 1:17)
meta <- meta %>% mutate(Genotype = factor(Genotype, levels = ordered_genotypes))

cols <- list(
  "HPeV1" = "#4E79A7",
  "HPeV2" = "#F28E2B",
  "HPeV3" = "#59A14F",
  "HPeV4" = "#E15759",
  "HPeV5" = "#76B7B2",
  "HPeV6" = "#EDC948",
  "HPeV7" = "#B07AA1",
  "HPeV8" = "#FF9DA7",
  "HPeV9" = "#9C755F",
  "HPeV10" = "#BAB0AC",
  "HPeV11" = "#8CD17D",
  "HPeV12" = "#499894",
  "HPeV13" = "#D37295",
  "HPeV14" = "#86BCB6",
  "HPeV15" = "#A0CBE8",
  "HPeV16" = "#FFBE7D",
  "HPeV17" = "#8A89A6"
)

genotype_counts <- meta %>% count(Genotype, name = "VP1_sequences") %>% filter(!is.na(Genotype))
print(genotype_counts)

genotypes_present <- ordered_genotypes[ordered_genotypes %in% as.character(unique(na.omit(meta$Genotype)))]

p_tree <- ggtree(tree_rooted, linewidth = 0.20) %<+% meta +
  geom_tippoint(aes(color = Genotype), size = 0.75, alpha = 0.95) +
  scale_color_manual(values = cols, breaks = genotypes_present, limits = genotypes_present, drop = TRUE, name = "Genotype") +
  geom_treescale(x = 0, y = 0, width = 0.1, fontsize = 3, linesize = 0.5) +
  theme_tree() +
  theme(legend.position = "right", legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 8), legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"), plot.margin = margin(10, 20, 10, 10)) +
  guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3, alpha = 1)))

p_tree