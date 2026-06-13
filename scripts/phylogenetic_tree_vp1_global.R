############################################################
# PARECHOVIRUS VP1 TREE
# PANEL A = COUNTRY
# PANEL B = GENOTYPE
# SAME COLORS USED IN ALL FIGURES
############################################################

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)


# READ TREE


tree <- read.nexus(
  "/Users/hugo/Desktop/vp1_parachovirus/vp1_headers_updated.tree"
)

if(inherits(tree, "multiPhylo")){
  tree <- tree[[1]]
}


# EXTRACT METADATA


tips <- gsub(
  "'",
  "",
  tree$tip.label
)

fields <- str_split_fixed(
  tips,
  "\\|",
  5
)

meta <- data.frame(
  label     = tree$tip.label,
  accession = fields[,1],
  species   = fields[,2],
  genotype  = fields[,3],
  country   = fields[,4],
  year      = fields[,5],
  stringsAsFactors = FALSE
)


# CLEAN


meta$country <- gsub(
  "_",
  " ",
  meta$country
)

meta$country <- trimws(meta$country)
meta$genotype <- trimws(meta$genotype)


# COUNTRY ORDER


meta$country <- factor(
  meta$country,
  levels = names(
    sort(
      table(meta$country),
      decreasing = TRUE
    )
  )
)


# GENOTYPE ORDER


genotype_levels <- c(
  paste0("HPeV", 1:17),
  "Unknown"
)

genotype_levels <- genotype_levels[
  genotype_levels %in% unique(meta$genotype)
]

meta$genotype <- factor(
  meta$genotype,
  levels = genotype_levels
)


# COUNTRY PALETTE

country_levels <- levels(meta$country)

country_cols <- c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728",
  "#9467bd","#8c564b","#e377c2","#7f7f7f",
  "#bcbd22","#17becf","#393b79","#637939",
  "#8c6d31","#843c39","#7b4173","#3182bd",
  "#31a354","#756bb1","#636363","#e6550d",
  "#969696","#6baed6","#74c476","#fd8d3c",
  "#9e9ac8","#bdbdbd","#c7e9c0","#fdae6b",
  "#c994c7","#5254a3","#6b6ecf","#9c9ede",
  "#637939","#8ca252","#b5cf6b"
)

country_cols <- rep(
  country_cols,
  length.out = length(country_levels)
)

names(country_cols) <- country_levels


# GENOTYPE PALETTE


genotype_cols <- c(
  "HPeV1"  = "#1f77b4",
  "HPeV2"  = "#ff7f0e",
  "HPeV3"  = "#e41a1c",
  "HPeV4"  = "#4daf4a",
  "HPeV5"  = "#984ea3",
  "HPeV6"  = "#a65628",
  "HPeV7"  = "#f781bf",
  "HPeV8"  = "#00bfc4",
  "HPeV9"  = "#66c2a5",
  "HPeV10" = "#ffd92f",
  "HPeV11" = "#8da0cb",
  "HPeV12" = "#e78ac3",
  "HPeV13" = "#a6d854",
  "HPeV14" = "#fc8d62",
  "HPeV15" = "#e5c494",
  "HPeV16" = "#1b9e77",
  "HPeV17" = "#d95f02",
  "Unknown" = "#7F7F7F"
)

genotype_cols <- genotype_cols[
  names(genotype_cols) %in% genotype_levels]


# FIGURE A
# COUNTRY


figA <- ggtree(tree) %<+% meta +
  geom_tippoint(
    aes(color = country),
    size = 2
  ) +
  scale_color_manual(
    values = country_cols
  ) +
  theme_tree2() +
  labs(color = "Country") +
  ggtitle("A") +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    )
  )


# FIGURE B
# GENOTYPE

figB <- ggtree(tree) %<+% meta +
  geom_tippoint(
    aes(color = genotype),
    size = 2
  ) +
  scale_color_manual(
    values = genotype_cols,
    drop = FALSE
  ) +
  theme_tree2() +
  labs(color = "Genotype") +
  ggtitle("B") +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    )
  )


# COMBINED FIGURE


final_plot <- figA | figB


# SHOW


final_plot


# SAVE

ggsave(
  "Parechovirus_Tree_Country_Genotype.pdf",
  final_plot,
  width = 18,
  height = 10
)

ggsave(
  "Parechovirus_Tree_Country_Genotype.png",
  final_plot,
  width = 18,
  height = 10,
  dpi = 600
)

