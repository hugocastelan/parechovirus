# PeV-A VP1 PHYLOGENY - FIGURE 2

library(ape)
library(phytools)
library(ggtree)
library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)

# Input files 
tree_file <- paste0(
  "/Users/hugo/Desktop/vp1_parachovirus/",
  "manual_v3_coverage50_ambiguity10_noduplicates_relabel_v5.1.tree"
)
metadata_file <- paste0(
  "/Users/hugo/Desktop/vp1_parachovirus/",
  "manual_v3_coverage50_ambiguity10_noduplicates_relabel_v5_metadata_A18_A19.csv"
)

# Read data
tree <- read.tree(tree_file)
meta <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

cat("Tree tips:", Ntip(tree), "\n")
cat("Internal nodes:", tree$Nnode, "\n")
cat("Metadata rows:", nrow(meta), "\n")
cat("Node support labels:", length(tree$node.label), "\n")

# Midpoint root 
tree <- midpoint.root(tree)

# Prepare metadata
meta$accession_clean <- sub("\\.[0-9]+$", "", trimws(meta$accession))
meta$genotype_plot <- trimws(meta$genotype)
meta$genotype_plot <- gsub("^HPeV[-_ ]?([0-9]+)$", "PeV-A\\1", meta$genotype_plot)
meta$genotype_plot <- gsub("^PeV-A[-_ ]?([0-9]+)$", "PeV-A\\1", meta$genotype_plot)

genotype_order <- paste0("PeV-A", 1:19)
meta$genotype_plot <- factor(meta$genotype_plot, levels = genotype_order)

# Identify SRA 
meta$is_sra <- grepl("^(SRR|ERR|DRR)", meta$accession_clean)
cat("Total SRA sequences:", sum(meta$is_sra), "\n")

# Genotype colors 
genotype_colors <- c(
  "PeV-A1"="#4C78A8","PeV-A2"="#F58518","PeV-A3"="#2CA02C","PeV-A4"="#E45756",
  "PeV-A5"="#72B7B2","PeV-A6"="#ECA400","PeV-A7"="#B279A2","PeV-A8"="#FF9DA6",
  "PeV-A9"="#9C755F","PeV-A10"="#8C564B","PeV-A11"="#59A14F","PeV-A12"="#76B7B2",
  "PeV-A13"="#EDC948","PeV-A14"="#4E9F95","PeV-A15"="#E15759","PeV-A16"="#F28E2B",
  "PeV-A17"="#6B6ECF","PeV-A18"="#B07AA1","PeV-A19"="#17BECF"
)

#Tree data 
p_base <- ggtree(tree, color = "grey65", linewidth = 0.28)
tree_data <- p_base$data

tree_data$accession_clean <- NA_character_
tree_data$accession_clean[tree_data$isTip] <- sub(
  "\\|.*$", "", tree_data$label[tree_data$isTip]
)
tree_data$accession_clean[tree_data$isTip] <- sub(
  "\\.[0-9]+$", "", tree_data$accession_clean[tree_data$isTip]
)

tree_data <- tree_data %>% left_join(meta, by = "accession_clean")

cat("Tips with metadata:",
    sum(tree_data$isTip & !is.na(tree_data$genotype_plot)), "\n")
cat("Tips without metadata:",
    sum(tree_data$isTip & is.na(tree_data$genotype_plot)), "\n")

# Node support
tree_data$support <- NA_real_
tree_data$support[!tree_data$isTip] <- suppressWarnings(
  as.numeric(tree_data$label[!tree_data$isTip])
)
tree_data$support_percent <- NA_real_
tree_data$support_percent[!tree_data$isTip] <- round(
  tree_data$support[!tree_data$isTip] * 100
)

support_threshold <- 95
support_data <- tree_data %>%
  filter(!isTip, !is.na(support_percent), support_percent >= support_threshold)

cat("Nodes with support >=", support_threshold, ":", nrow(support_data), "\n")

# Counts for Panel A 
counts_all <- tree_data %>%
  filter(isTip, !is.na(genotype_plot)) %>%
  count(genotype_plot, .drop = FALSE)

legend_labels_A <- setNames(
  paste0(counts_all$genotype_plot, " (n = ", counts_all$n, ")"),
  as.character(counts_all$genotype_plot)
)

#  Panel A 
pA <- ggtree(tree, color = "grey65", linewidth = 0.28)
pA$data <- tree_data

pA <- pA +
  geom_point(
    data = tree_data %>% filter(isTip, !is.na(genotype_plot)),
    aes(x = x, y = y, color = genotype_plot),
    inherit.aes = FALSE, size = 0.85, alpha = 0.95
  ) +
  geom_text(
    data = support_data,
    aes(x = x, y = y, label = support_percent),
    inherit.aes = FALSE, color = "black", size = 1.65,
    fontface = "plain", hjust = -0.20, vjust = -0.15, check_overlap = TRUE
  ) +
  scale_color_manual(
    values = genotype_colors,
    breaks = genotype_order,
    labels = legend_labels_A[genotype_order],
    drop = FALSE,
    name = "Parechovirus A genotype"
  ) +
  geom_treescale(width = 0.1, fontsize = 3) +
  theme_tree() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.height = unit(0.30, "cm"),
    plot.margin = margin(5, 15, 5, 5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

#SRA data for Panel B
sra_tip_data <- tree_data %>%
  filter(isTip, !is.na(is_sra), is_sra, !is.na(genotype_plot))

sra_counts <- sra_tip_data %>% count(genotype_plot, .drop = FALSE) %>% filter(n > 0)

sra_genotypes <- genotype_order[
  genotype_order %in% as.character(sra_counts$genotype_plot)
]

legend_labels_B <- setNames(
  paste0(sra_counts$genotype_plot, " (n = ", sra_counts$n, ")"),
  as.character(sra_counts$genotype_plot)
)

# ---------- Bars for Panel B ----------
minimum_sra_for_bar <- 4
bar_genotypes <- sra_counts %>%
  filter(n >= minimum_sra_for_bar) %>%
  pull(genotype_plot) %>% as.character()

bar_positions <- sra_tip_data %>%
  filter(as.character(genotype_plot) %in% bar_genotypes) %>%
  group_by(genotype_plot) %>%
  summarise(
    ymin = min(y), ymax = max(y),
    ymid = (min(y) + max(y)) / 2,
    n = n(), .groups = "drop"
  )

xmin_tree <- min(tree_data$x, na.rm = TRUE)
xmax_tree <- max(tree_data$x, na.rm = TRUE)
tree_width <- xmax_tree - xmin_tree

if (nrow(bar_positions) > 0) {
  bar_positions <- bar_positions %>%
    arrange(ymid) %>%
    mutate(
      bar_number = row_number(),
      bar_x = xmax_tree + tree_width * (0.035 + (bar_number - 1) * 0.025),
      text_x = bar_x + tree_width * 0.012
    )
}

if (nrow(bar_positions) > 0) {
  plot_xmax <- max(bar_positions$text_x) + tree_width * 0.04
} else {
  plot_xmax <- xmax_tree + tree_width * 0.05
}

#Panel B 
pB <- ggtree(tree, color = "grey65", linewidth = 0.28)
pB$data <- tree_data

pB <- pB +
  geom_point(
    data = sra_tip_data,
    aes(x = x, y = y, color = genotype_plot),
    inherit.aes = FALSE, size = 2.2, alpha = 1
  ) +
  geom_segment(
    data = bar_positions,
    aes(x = bar_x, xend = bar_x, y = ymin, yend = ymax, color = genotype_plot),
    inherit.aes = FALSE, linewidth = 0.8, show.legend = FALSE
  ) +
  geom_text(
    data = bar_positions,
    aes(x = text_x, y = ymid, label = genotype_plot, color = genotype_plot),
    inherit.aes = FALSE, angle = 90, size = 3,
    hjust = 0.5, vjust = 0.5, show.legend = FALSE
  ) +
  scale_color_manual(
    values = genotype_colors,
    breaks = sra_genotypes,
    labels = legend_labels_B[sra_genotypes],
    drop = TRUE,
    name = "Parechovirus A genotype (SRA)"
  ) +
  geom_treescale(width = 0.1, fontsize = 3) +
  coord_cartesian(xlim = c(xmin_tree, plot_xmax), clip = "off") +
  theme_tree() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.35, "cm"),
    plot.margin = margin(5, 25, 5, 5)
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 4)))

# Combine 
fig <- pA + pB +
  plot_layout(widths = c(1, 1.12)) +
  plot_annotation(tag_levels = "A")

#Save 
output_pdf <- paste0(
  "/Users/hugo/Desktop/vp1_parachovirus/",
  "Figure_PeV-A_VP1_FastTree_support.pdf"
)

ggsave(
  filename = output_pdf, plot = fig,
  width = 15, height = 8, units = "in", device = cairo_pdf
)

