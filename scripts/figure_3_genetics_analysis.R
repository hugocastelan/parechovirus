
# PeV-A VP1 GENETIC DIVERSI
# 1. LIBRARIES

library(ape)
library(pegas)
library(ggplot2)
library(stringr)
library(vegan)
library(factoextra)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggdendro)


# 2. SETTINGS

alignment_file <- paste0(
  "/Users/hugo/Desktop/vp1_parachovirus/",
  "manual_v3_coverage50_ambiguity10_noduplicates_relabel_v5.1.fasta"
)

ordered_genotypes <- paste0("PeV-A", 1:19)
original_genotypes <- paste0("HPeV", 1:19)

colors <- c(
  "PeV-A1"  = "#4E79A7",
  "PeV-A2"  = "#F28E2B",
  "PeV-A3"  = "#59A14F",
  "PeV-A4"  = "#E15759",
  "PeV-A5"  = "#76B7B2",
  "PeV-A6"  = "#EDC948",
  "PeV-A7"  = "#B07AA1",
  "PeV-A8"  = "#FF9DA7",
  "PeV-A9"  = "#9C755F",
  "PeV-A10" = "#BAB0AC",
  "PeV-A11" = "#8CD17D",
  "PeV-A12" = "#499894",
  "PeV-A13" = "#D37295",
  "PeV-A14" = "#86BCB6",
  "PeV-A15" = "#A0CBE8",
  "PeV-A16" = "#FFBE7D",
  "PeV-A17" = "#8A89A6",
  "PeV-A18" = "#B07AA1",
  "PeV-A19" = "#17BECF"
)


# 3. READ ALIGNMENT

alignment <- read.dna(alignment_file, format = "fasta")
sequence_names <- rownames(alignment)

cat("\n========================================\n")
cat("ALIGNMENT\n")
cat("========================================\n")
cat("Total sequences:", length(sequence_names), "\n")


# 4. EXTRACT METADATA FROM FASTA HEADERS

split_names <- str_split(sequence_names, "\\|")

metadata <- data.frame(
  name = sequence_names,
  Genotype_original = sapply(split_names, function(x) {
    if(length(x) >= 3) x[3] else NA
  }),
  stringsAsFactors = FALSE
)


# 5. CLEAN GENOTYPE NAMES

metadata$Genotype_original <- trimws(metadata$Genotype_original)

metadata <- metadata %>%
  filter(
    !is.na(Genotype_original),
    Genotype_original != "NA",
    Genotype_original != "",
    Genotype_original != "Unknown",
    Genotype_original %in% original_genotypes
  )

metadata <- metadata %>%
  mutate(Genotype = sub("^HPeV", "PeV-A", Genotype_original))

metadata$Genotype <- factor(metadata$Genotype, levels = ordered_genotypes)


# 6. FILTER ALIGNMENT

alignment <- alignment[metadata$name, ]
metadata <- droplevels(metadata)
present_genotypes <- levels(metadata$Genotype)


# 7. CHECK GENOTYPES

cat("GENOTYPES\n")

print(table(metadata$Genotype))

cat("\nGenotypes present:\n")
print(present_genotypes)

cat("\nPeV-A18 / PeV-A19 sequences:\n")
print(metadata %>% filter(Genotype %in% c("PeV-A18", "PeV-A19")))


# 8. PANEL A: PCoA + PERMANOVA

distance_matrix <- as.matrix(
  dist.dna(alignment, model = "raw", pairwise.deletion = TRUE)
)

metadata <- metadata[match(rownames(distance_matrix), metadata$name), ]

adonis_results <- adonis2(
  as.dist(distance_matrix) ~ Genotype,
  data = metadata,
  permutations = 999
)

r2_value <- round(adonis_results$R2[1], 3)
p_value <- adonis_results$`Pr(>F)`[1]

p_label <- ifelse(
  p_value <= 0.001,
  "p ≤ 0.001",
  paste0("p = ", round(p_value, 4))
)

adonis_annotation <- paste0(
  "PERMANOVA (ADONIS)\n",
  "R² = ", r2_value, "\n",
  p_label
)

pcoa_results <- cmdscale(as.dist(distance_matrix), k = 2, eig = TRUE)

pcoa_data <- data.frame(
  Axis1 = pcoa_results$points[, 1],
  Axis2 = pcoa_results$points[, 2],
  Genotype = metadata$Genotype
)

positive_eigen <- pcoa_results$eig[pcoa_results$eig > 0]
variance_explained <- round(100 * pcoa_results$eig / sum(positive_eigen), 2)

pcoa_plot <- ggplot(
  pcoa_data,
  aes(x = Axis1, y = Axis2, color = Genotype)
) +
  geom_point(size = 2.6, alpha = 0.75) +
  scale_color_manual(
    values = colors,
    breaks = ordered_genotypes,
    limits = ordered_genotypes,
    drop = TRUE
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "A. Sequence diversity",
    x = paste0("PCoA 1 (", variance_explained[1], "%)"),
    y = paste0("PCoA 2 (", variance_explained[2], "%)"),
    color = "Genotype"
  ) +
  annotate(
    "text",
    x = min(pcoa_data$Axis1) * 0.85,
    y = max(pcoa_data$Axis2) * 0.85,
    label = adonis_annotation,
    hjust = 0, vjust = 1, size = 3.3,
    fontface = "italic", color = "black"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 6),
    plot.title = element_text(face = "plain")
  )


# 9. PANEL B: GENOTYPE DENDROGRAM

genotype_levels <- levels(droplevels(metadata$Genotype))

genotype_distances <- matrix(
  0,
  nrow = length(genotype_levels),
  ncol = length(genotype_levels),
  dimnames = list(genotype_levels, genotype_levels)
)

for(i in genotype_levels){
  for(j in genotype_levels){
    ids_i <- metadata$name[metadata$Genotype == i]
    ids_j <- metadata$name[metadata$Genotype == j]
    genotype_distances[i, j] <- mean(
      distance_matrix[ids_i, ids_j],
      na.rm = TRUE
    )
  }
}

hierarchical_clustering <- hclust(
  as.dist(genotype_distances),
  method = "average"
)

dendrogram <- as.dendrogram(hierarchical_clustering)
dendrogram_data <- dendro_data(dendrogram)

label_data <- dendrogram_data$labels
label_data$Genotype <- factor(label_data$label, levels = ordered_genotypes)
label_data$color <- colors[as.character(label_data$Genotype)]

max_y <- max(dendrogram_data$segments$y)
offset_value <- max_y * 0.02

dendrogram_plot <- ggplot() +
  geom_segment(
    data = dendrogram_data$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.6, color = "black"
  ) +
  geom_text(
    data = label_data,
    aes(x = x, y = 0 - offset_value, label = label, color = Genotype),
    size = 4.2, angle = 90, hjust = 1, vjust = 0.5, fontface = "bold"
  ) +
  scale_color_manual(values = colors) +
  labs(
    title = "B. Genotype clustering",
    y = "Mean genetic distance",
    x = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.20, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 50, 10),
    plot.title = element_text(face = "plain", size = 12)
  )


# 10. POPULATION GENETICS

results_table <- data.frame()

for(genotype in present_genotypes){
  ids <- metadata$name[metadata$Genotype == genotype]
  genotype_alignment <- alignment[ids, ]
  sequence_count <- length(ids)
  if(sequence_count < 3) next
  haplotypes <- haplotype(genotype_alignment)
  haplotype_count <- nrow(haplotypes)
  nucleotide_diversity <- nuc.div(genotype_alignment)
  tajima_test <- tryCatch(
    tajima.test(genotype_alignment),
    error = function(e) NA
  )
  tajima_d_value <- if(is.list(tajima_test)) tajima_test$D else NA
  results_table <- rbind(
    results_table,
    data.frame(
      Genotype = genotype,
      Haplotypes_Count = haplotype_count,
      Nucleotide_Diversity = nucleotide_diversity,
      Tajima_D = tajima_d_value
    )
  )
}

results_table$Genotype <- factor(
  results_table$Genotype,
  levels = ordered_genotypes
)

cat("\n========================================\n")
cat("POPULATION GENETICS RESULTS\n")
cat("========================================\n")
print(results_table)


# 11FORMAT

results_long <- results_table %>%
  pivot_longer(
    cols = c(Haplotypes_Count, Nucleotide_Diversity, Tajima_D),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = case_when(
      Metric == "Haplotypes_Count" ~ "C. Haplotype Richness (H)",
      Metric == "Nucleotide_Diversity" ~ "D. Nucleotide Diversity (π)",
      Metric == "Tajima_D" ~ "E. Tajima's D Index"
    ),
    Metric = factor(
      Metric,
      levels = c(
        "C. Haplotype Richness (H)",
        "D. Nucleotide Diversity (π)",
        "E. Tajima's D Index"
      )
    )
  )


# 12. PANELS C, D, E

lollipop_plots <- ggplot(
  results_long,
  aes(x = Genotype, y = Value, color = Genotype)
) +
  geom_segment(
    aes(xend = Genotype, y = 0, yend = Value),
    linewidth = 0.6
  ) +
  geom_point(size = 3) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = colors, drop = TRUE) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )


# 13. FINAL FIGURE

final_layout <- "
AABB
LLLL
"

master_figure <- pcoa_plot + dendrogram_plot + lollipop_plots +
  plot_layout(
    design = final_layout,
    heights = c(1.2, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot_title = element_text(face = "plain", size = 13),
      plot_subtitle = element_text(color = "dimgray", size = 9.5)
    )
  )

master_figure

# 14. SAVE

output_file <- paste0(
  "/Users/hugo/Desktop/vp1_parachovirus/",
  "PeV-A_VP1_diversity_analysis_A1_A19.pdf"
)
