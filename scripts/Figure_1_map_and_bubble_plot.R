# HPeV GLOBAL DISTRIBUTION
# PANEL A = WORLD MAP WITH DONUTS + SEQUENCE COUNTS
# PANEL B = TEMPORAL BUBBLE PLOT

rm(list = ls())


library(dplyr)
library(tidyr)
library(maps)
library(ggfree)

# 2. LOAD METADATA
metadata <- read.csv("/Users/hugo/Desktop/vp1_parachovirus/manual_v3_coverage50_ambiguity10_noduplicates_relabel_v5_metadata.csv",
  stringsAsFactors = FALSE)


# 3. CLEAN METADATA

metadata2 <- metadata %>%
  filter(
    !is.na(country), !is.na(genotype),
    country != "", genotype != "",
    !tolower(country) %in% c("unknown", "na"),
    !tolower(genotype) %in% c("unknown", "hpev_unknown", "hpevunknown", "na")
  ) %>%
  mutate(
    country = recode(
      country,
      "Bangladeshi" = "Bangladesh",
      "Czech" = "Czech Republic",
      "Ivoire Coast" = "Cote d'Ivoire",
      "United States" = "USA",
      "United States of America" = "USA",
      "UK" = "United Kingdom"
    ),
    year = suppressWarnings(floor(as.numeric(year)))
  )


# 4. GENOTYPES

genotypes <- paste0("HPeV", 1:17)
metadata2$genotype <- factor(metadata2$genotype, levels = genotypes)


# 5. GENOTYPE COLORS

genotype_colors <- c(
  HPeV1  = "#4477AA", HPeV2  = "#EE6677", HPeV3  = "#CC3311",
  HPeV4  = "#228833", HPeV5  = "#009988", HPeV6  = "#E69F00",
  HPeV7  = "#AA4499", HPeV8  = "#CC79A7", HPeV9  = "#8C6D31",
  HPeV10 = "#999999", HPeV11 = "#0072B2", HPeV12 = "#D55E00",
  HPeV13 = "#117733", HPeV14 = "#882255", HPeV15 = "#661100",
  HPeV16 = "#56B4E9", HPeV17 = "#D4C34F", HPeV18 = "#332288",
  HPeV19 = "#009E73"
)


# 6. TOTAL SEQUENCES PER GENOTYPE 

genotype_totals <- metadata2 %>%
  filter(!is.na(genotype)) %>%
  count(genotype, name = "n")

genotype_labels <- sapply(
  genotypes,
  function(g) {
    n <- genotype_totals$n[genotype_totals$genotype == g]
    if (length(n) == 0) n <- 0
    sprintf("%s (%d)", g, n)
  }
)


# 7. COUNTRY COORDINATES
# long / lat = approximate true geographic position

country_coordinates <- data.frame(
  country = c(
    "Australia", "Bangladesh", "Belarus", "Belgium", "Bolivia", "Brazil",
    "Bulgaria", "China", "Cote d'Ivoire", "Czech Republic", "Denmark",
    "Ethiopia", "Finland", "France", "Germany", "Ghana", "Greece",
    "Hong Kong", "Hungary", "India", "Ireland", "Japan", "Malawi",
    "Mexico", "Netherlands", "Norway", "Pakistan", "South Korea",
    "South Sudan", "Sri Lanka", "Taiwan", "Thailand", "United Kingdom",
    "USA", "Venezuela", "Vietnam"
  ),
  long = c(
    134.0, 90.3, 28.0, 4.7, -64.7, -51.9, 25.5, 104.2, -5.5, 15.5,
    9.5, 40.5, 26.0, 2.2, 10.4, -1.0, 22.0, 114.2, 19.5, 78.9,
    -8.0, 138.0, 34.3, -102.0, 5.5, 8.5, 69.3, 128.0, 30.2, 80.7,
    121.0, 101.0, -3.0, -100.0, -66.6, 108.0
  ),
  lat = c(
    -25.0, 23.7, 53.7, 50.8, -16.3, -14.2, 42.7, 35.9, 7.5, 49.8,
    56.0, 9.1, 64.0, 46.2, 51.2, 7.9, 39.0, 22.3, 47.2, 21.0,
    53.0, 37.0, -13.3, 23.5, 52.2, 61.0, 30.4, 36.0, 7.0, 7.9,
    23.7, 15.5, 54.0, 38.0, 6.4, 16.0
  ),
  stringsAsFactors = FALSE
)


# 8. INITIAL DISPLAY COORDINATES

country_coordinates <- country_coordinates %>%
  mutate(plot_long = long, plot_lat = lat, external = FALSE)


# 9. EXTERNAL DONUT POSITIONS FOR EUROPE

europe_positions <- data.frame(
  country = c(
    "Ireland", "United Kingdom", "Belgium", "Netherlands", "France",
    "Germany", "Denmark", "Norway", "Finland", "Czech Republic",
    "Hungary", "Bulgaria", "Greece", "Belarus"
  ),
  plot_long_new = c(
    -30, -20, -18, -8, -23, -6, 5, -9, 18, 14, 26, 37, 31, 38
  ),
  plot_lat_new = c(
    61, 69, 46, 73, 37, 41, 73, 76, 76, 68, 67, 60, 36, 72
  ),
  stringsAsFactors = FALSE
)


# 10. ADD EXTERNAL POSITIONS

country_coordinates <- country_coordinates %>%
  left_join(europe_positions, by = "country") %>%
  mutate(
    external = !is.na(plot_long_new),
    plot_long = ifelse(external, plot_long_new, long),
    plot_lat  = ifelse(external, plot_lat_new, lat)
  ) %>%
  select(country, long, lat, plot_long, plot_lat, external)


# 11. COUNTRY per GENOTYPE COUNTS

country_counts <- metadata2 %>%
  filter(!is.na(genotype)) %>%
  count(country, genotype, name = "n")


# 12. CHECK COUNTRIES

cat("\nCountries present in metadata:\n")
print(sort(unique(metadata2$country)))

cat("\nCountries without coordinates:\n")
print(setdiff(unique(country_counts$country), country_coordinates$country))


# 13. COMPLETE COUNTRY per GENOTYPE MATRIX

country_counts_complete <- country_counts %>%
  complete(
    country,
    genotype = factor(genotypes, levels = genotypes),
    fill = list(n = 0)
  ) %>%
  left_join(country_coordinates, by = "country") %>%
  filter(!is.na(plot_long), !is.na(plot_lat))


# 14. TOTAL SEQUENCES PER COUNTRY

country_totals <- country_counts_complete %>%
  group_by(country, long, lat, plot_long, plot_lat, external) %>%
  summarise(total = sum(n), .groups = "drop")

cat("\nSequences by country:\n")
print(country_totals %>% arrange(desc(total)))


# 15. TEMPORAL DATA

temporal_data <- metadata2 %>%
  filter(!is.na(year), year >= 1950, year <= 2026, !is.na(genotype)) %>%
  count(year, genotype, name = "n") %>%
  mutate(
    genotype = factor(genotype, levels = genotypes),
    genotype_position = as.numeric(genotype)
  )

cat("\nTemporal range:\n")
print(range(temporal_data$year, na.rm = TRUE))

cat("\nMaximum sequences in one genotype-year:\n")
print(max(temporal_data$n, na.rm = TRUE))


# 16. PANEL A — WORLD MAP plus DONUTS

draw_panel_A <- function() {
  par(mar = c(0, 0, 0, 0), xpd = NA, pty = "m")
  
  # WORLD MAP (reduced vertical range, no Antarctica)
  map(
    "world", fill = TRUE, col = "grey94", border = "white", lwd = 0.30,
    xlim = c(-175, 180), ylim = c(-38, 78), asp = NA
  )
  
  # Country borders
  map("world", add = TRUE, fill = FALSE, col = "grey82", lwd = 0.20)
  
  max_total <- max(country_totals$total, na.rm = TRUE)
  external_data <- country_totals %>% filter(external)
  
  # LEADER LINES
  if (nrow(external_data) > 0) {
    for (i in seq_len(nrow(external_data))) {
      segments(
        x0 = external_data$long[i], y0 = external_data$lat[i],
        x1 = external_data$plot_long[i], y1 = external_data$plot_lat[i],
        col = "grey45", lwd = 0.75
      )
      points(
        x = external_data$long[i], y = external_data$lat[i],
        pch = 21, bg = "white", col = "grey30", lwd = 0.70, cex = 0.65
      )
    }
  }
  
  # DRAW DONUTS
  for (i in seq_len(nrow(country_totals))) {
    current_country <- country_totals$country[i]
    
    values <- country_counts_complete %>%
      filter(country == current_country) %>%
      arrange(match(as.character(genotype), genotypes))
    
    donut_values <- values$n
    names(donut_values) <- as.character(values$genotype)
    donut_values <- donut_values[donut_values > 0]
    
    if (length(donut_values) == 0) next
    
    donut_colors <- genotype_colors[names(donut_values)]
    total_seq <- country_totals$total[i]
    
    # Donut size: sqrt scaling
    r1 <- 4.0 + 3.5 * sqrt(total_seq / max_total)
    r0 <- r1 * 0.56
    
    ggfree::ringplot(
      donut_values,
      x = country_totals$plot_long[i],
      y = country_totals$plot_lat[i],
      r0 = r0, r1 = r1,
      col = donut_colors, border = "white", lwd = 0.45,
      use.names = FALSE
    )
    
    text(
      x = country_totals$plot_long[i],
      y = country_totals$plot_lat[i],
      labels = total_seq, cex = 0.66, font = 2, col = "black"
    )
  }
  
  mtext("A", side = 3, adj = 0, line = -1.1, font = 2, cex = 1.50)
}


# 17. GENOTYPE LEGEND

draw_genotype_legend <- function() {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = genotype_labels,
    fill = genotype_colors[genotypes],
    border = NA,
    ncol = 10,       # 19 genotypes in two rows
    cex = 0.90,
    bty = "n",
    x.intersp = 0.65,
    y.intersp = 1.15
  )
}

# 18. PANEL B. TEMPORAL BUBBLE PLOT

draw_panel_B <- function() {
  par(mar = c(4.2, 4.8, 0.5, 4.8), xpd = NA)
  
  year_min <- min(temporal_data$year, na.rm = TRUE)
  year_max <- max(temporal_data$year, na.rm = TRUE)
  
  plot(
    NA,
    xlim = c(year_min - 1, year_max + 2),
    ylim = c(0.5, 19.5),
    xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n"
  )
  
  abline(h = 1:19, col = "grey90", lwd = 0.60)
  
  year_breaks <- seq(
    floor(year_min / 5) * 5,
    ceiling(year_max / 5) * 5,
    by = 5
  )
  
  abline(v = year_breaks, col = "grey94", lwd = 0.55)
  
  temporal_data$bubble_cex <- pmin(
    0.55 + 0.35 * sqrt(temporal_data$n),
    4.8
  )
  
  points(
    x = temporal_data$year,
    y = temporal_data$genotype_position,
    pch = 21,
    bg = genotype_colors[as.character(temporal_data$genotype)],
    col = adjustcolor("black", alpha.f = 0.50),
    lwd = 0.40,
    cex = temporal_data$bubble_cex
  )
  
  axis(1, at = year_breaks, labels = year_breaks, las = 2, cex.axis = 0.74, tck = -0.012)
  mtext("Year", side = 1, line = 3.3, cex = 1)
  
  axis(2, at = 1:19, labels = genotypes, las = 1, cex.axis = 0.76, tick = FALSE)
  
  mtext("B", side = 3, adj = 0, line = -0.35, font = 2, cex = 1.50)
  
  # BUBBLE SIZE LEGEND
  possible_breaks <- c(1, 5, 10, 25, 50, 100)
  size_breaks <- possible_breaks[
    possible_breaks <= max(temporal_data$n, na.rm = TRUE)
  ]
  
  if (length(size_breaks) > 0) {
    legend_cex <- pmin(0.55 + 0.35 * sqrt(size_breaks), 4.8)
    legend(
      x = year_max + 4,
      y = 18.5,
      legend = size_breaks,
      pch = 21,
      pt.bg = "grey80",
      col = "grey40",
      pt.cex = legend_cex,
      title = "Sequences",
      cex = 0.68,
      bty = "n",
      xpd = NA,
      y.intersp = 1.15
    )
  }
}

#
# 19. DRAW OF FIGURE
draw_complete_figure <- function() {
  layout(
    matrix(c(1, 2, 3), nrow = 3),
    heights = c(4.40, 0.95, 3.40)
  )
  draw_panel_A()
  draw_genotype_legend()
  draw_panel_B()
}
draw_complete_figure()


# 21. EXPORT PDF

pdf("HPeV_global_distribution_FINAL.pdf",  width = 15, height = 8,
  useDingbats = FALSE, family = "Helvetica")
draw_complete_figure()
dev.off()
