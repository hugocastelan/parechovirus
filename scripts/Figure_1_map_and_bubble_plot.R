
# A = Geographic distribution
# B = Temporal distribution by genotype


library(dplyr)
library(tidyr)
library(maps)
library(ggfree)

#read metadata 
metadata <- read.csv(
  "/Users/hugo/Desktop/vp1_parachovirus/manual_v3_coverage50_ambiguity10_noduplicates_relabel_v5_metadata_A18_A19.csv",
  stringsAsFactors = FALSE
)

# Standardize 
metadata <- metadata %>%
  mutate(
    accession = as.character(accession),
    virus     = as.character(virus),
    genotype  = as.character(genotype),
    country   = as.character(country),
    year      = as.character(year),
    accession = sub("\\.[0-9]+$", "", accession)
  )

# Clean 
metadata_clean <- metadata %>%
  filter(
    !is.na(genotype), genotype != "",
    !tolower(genotype) %in% c("unknown","hpev_unknown","hpevunknown","na")
  ) %>%
  mutate(
    genotype = trimws(genotype),
    genotype = gsub("^PeV-A", "HPeV", genotype, ignore.case = TRUE),
    genotype = gsub("^PeVA",  "HPeV", genotype, ignore.case = TRUE),
    genotype = gsub("^HPeV-", "HPeV", genotype, ignore.case = TRUE),
    country = recode(country,
                     "Bangladeshi" = "Bangladesh",
                     "Czech" = "Czech Republic",
                     "Ivoire Coast" = "Cote d'Ivoire",
                     "United States" = "USA",
                     "United States of America" = "USA",
                     "UK" = "United Kingdom"
    ),
    year = suppressWarnings(floor(as.numeric(year)))
  )

# Genotypes 
genotypes <- paste0("HPeV", 1:19)
metadata_clean$genotype <- factor(metadata_clean$genotype, levels = genotypes)
genotype_display <- setNames(paste0("PeV-A", 1:19), genotypes)

#Colors 
genotype_colors <- c(
  "HPeV1"="#4C78A8","HPeV2"="#F58518","HPeV3"="#2CA02C","HPeV4"="#E45756",
  "HPeV5"="#72B7B2","HPeV6"="#ECA400","HPeV7"="#B279A2","HPeV8"="#FF9DA6",
  "HPeV9"="#9C755F","HPeV10"="#8C564B","HPeV11"="#59A14F","HPeV12"="#76B7B2",
  "HPeV13"="#EDC948","HPeV14"="#4E9F95","HPeV15"="#E15759","HPeV16"="#F28E2B",
  "HPeV17"="#6B6ECF","HPeV18"="#B07AA1","HPeV19"="#17BECF"
)

#Check A18/A19 

print(metadata_clean %>%
        filter(accession %in% c("OR513951","MH339678")) %>%
        select(accession, genotype, country, year))

#Totals 
genotype_totals <- metadata_clean %>%
  filter(!is.na(genotype)) %>%
  count(genotype, name = "n", .drop = FALSE)

present_genotypes <- genotype_totals %>%
  filter(n > 0) %>% pull(genotype) %>% as.character()
present_genotypes <- genotypes[genotypes %in% present_genotypes]

# legend labels 
genotype_labels <- sapply(genotypes, function(g) {
  n <- genotype_totals$n[as.character(genotype_totals$genotype) == g]
  if (length(n) == 0) n <- 0
  sprintf("%s (%d)", genotype_display[g], n)
})

# country coordinates 
country_coordinates <- data.frame(
  country = c(
    "Australia","Bangladesh","Belarus","Belgium","Bolivia","Brazil","Bulgaria",
    "China","Cote d'Ivoire","Czech Republic","Denmark","Ethiopia","Finland",
    "France","Germany","Ghana","Greece","Hong Kong","Hungary","India","Ireland",
    "Japan","Malawi","Mexico","Netherlands","Norway","Pakistan","South Korea",
    "South Sudan","Sri Lanka","Taiwan","Thailand","United Kingdom","USA",
    "Venezuela","Vietnam"
  ),
  long = c(134.0,90.3,28.0,4.7,-64.7,-51.9,25.5,104.2,-5.5,15.5,9.5,40.5,
           26.0,2.2,10.4,-1.0,22.0,114.2,19.5,78.9,-8.0,138.0,34.3,-102.0,
           5.5,8.5,69.3,128.0,30.2,80.7,121.0,101.0,-3.0,-100.0,-66.6,108.0),
  lat = c(-25.0,23.7,53.7,50.8,-16.3,-14.2,42.7,35.9,7.5,49.8,56.0,9.1,
          64.0,46.2,51.2,7.9,39.0,22.3,47.2,21.0,53.0,37.0,-13.3,23.5,
          52.2,61.0,30.4,36.0,7.0,7.9,23.7,15.5,54.0,38.0,6.4,16.0),
  stringsAsFactors = FALSE
)

country_coordinates <- country_coordinates %>%
  mutate(plot_long = long, plot_lat = lat, external = FALSE)

# external Europe positions
europe_positions <- data.frame(
  country = c("Ireland","United Kingdom","Belgium","Netherlands","France",
              "Germany","Denmark","Norway","Finland","Czech Republic",
              "Hungary","Bulgaria","Greece","Belarus"),
  plot_long_new = c(-30,-20,-18,-8,-23,-6,5,-9,18,14,26,37,31,38),
  plot_lat_new  = c(61,69,46,73,37,41,73,76,76,68,67,60,36,72),
  stringsAsFactors = FALSE
)

country_coordinates <- country_coordinates %>%
  left_join(europe_positions, by = "country") %>%
  mutate(
    external  = !is.na(plot_long_new),
    plot_long = ifelse(external, plot_long_new, plot_long),
    plot_lat  = ifelse(external, plot_lat_new,  plot_lat)
  ) %>%
  select(country, long, lat, plot_long, plot_lat, external)

# Map data
metadata_map <- metadata_clean %>%
  filter(!is.na(country), country != "",
         !tolower(country) %in% c("unknown","na"),
         !is.na(genotype))

country_counts <- metadata_map %>% count(country, genotype, name = "n")

country_counts_complete <- country_counts %>%
  complete(country, genotype = factor(genotypes, levels = genotypes),
           fill = list(n = 0)) %>%
  left_join(country_coordinates, by = "country") %>%
  filter(!is.na(plot_long), !is.na(plot_lat))

country_totals <- country_counts_complete %>%
  group_by(country, long, lat, plot_long, plot_lat, external) %>%
  summarise(total = sum(n), .groups = "drop")

# Temporal data 
metadata_temporal <- metadata_clean %>%
  filter(!is.na(genotype), !is.na(year), year >= 1950, year <= 2026)

temporal_data <- metadata_temporal %>%
  count(year, genotype, name = "n", .drop = FALSE) %>%
  filter(n > 0) %>%
  mutate(
    genotype = factor(genotype, levels = genotypes),
    genotype_position = as.numeric(genotype)
  )

#Size functions 
max_country_total <- max(country_totals$total, na.rm = TRUE)

donut_radius <- function(n) 4.5 + 5.0 * sqrt(n / max_country_total)
bubble_cex   <- function(n) pmin(0.65 + 0.42 * sqrt(n), 5.5)

# Panel A
draw_panel_A <- function() {
  par(mar = c(0.1,0.1,0.1,0.1), xpd = NA, pty = "m")
  
  map("world", fill = TRUE, col = "grey94", border = "white", lwd = 0.30,
      xlim = c(-168,165), ylim = c(-40,80), asp = NA)
  map("world", add = TRUE, fill = FALSE, col = "grey82", lwd = 0.20)
  
  external_data <- country_totals %>% filter(external)
  if (nrow(external_data) > 0) {
    for (i in seq_len(nrow(external_data))) {
      segments(external_data$long[i], external_data$lat[i],
               external_data$plot_long[i], external_data$plot_lat[i],
               col = "grey45", lwd = 0.8)
      points(external_data$long[i], external_data$lat[i],
             pch = 21, bg = "white", col = "grey30", lwd = 0.7, cex = 0.65)
    }
  }
  
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
    r1 <- donut_radius(total_seq); r0 <- r1 * 0.56
    
    ggfree::ringplot(donut_values,
                     x = country_totals$plot_long[i],
                     y = country_totals$plot_lat[i],
                     r0 = r0, r1 = r1,
                     col = donut_colors, border = "white", lwd = 0.45,
                     use.names = FALSE)
    
    text(country_totals$plot_long[i], country_totals$plot_lat[i],
         labels = total_seq, cex = 0.72, font = 2, col = "black")
  }
  
  mtext("A", side = 3, adj = 0, line = -1.2, font = 2, cex = 1.45)
  
  size_candidates <- c(10,50,100,250,500,1000)
  donut_breaks <- size_candidates[size_candidates <= max_country_total]
  if (length(donut_breaks) == 0) donut_breaks <- 1
  legend_values <- tail(donut_breaks, min(4, length(donut_breaks)))
  
  legend_x <- 145; legend_y <- 15
  text(legend_x, legend_y + 15, "Sequences", font = 2, cex = 0.78)
  
  for (j in seq_along(legend_values)) {
    n <- legend_values[j]
    rr <- donut_radius(n) * 0.55
    yy <- legend_y - (j - 1) * 12
    symbols(legend_x, yy, circles = rr, inches = FALSE, add = TRUE,
            bg = "white", fg = "grey35")
    text(legend_x + 12, yy, labels = n, adj = 0, cex = 0.72)
  }
}

# Genotype legend
draw_genotype_legend <- function() {
  par(mar = c(0,0,0,0))
  plot.new()
  legend("center",
         legend = genotype_labels[present_genotypes],
         fill = genotype_colors[present_genotypes],
         border = NA, ncol = 10, cex = 0.82, bty = "n",
         x.intersp = 0.60, y.intersp = 1.20)
}

#Panel B 
draw_panel_B <- function() {
  par(mar = c(4.3,5.5,0.7,7.5), xpd = NA)
  
  year_min <- min(temporal_data$year, na.rm = TRUE)
  year_max <- max(temporal_data$year, na.rm = TRUE)
  
  plot(NA, xlim = c(year_min - 1, year_max + 2), ylim = c(0.5, 19.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  
  abline(h = 1:19, col = "grey90", lwd = 0.60)
  
  year_breaks <- seq(floor(year_min/5)*5, ceiling(year_max/5)*5, by = 5)
  abline(v = year_breaks, col = "grey94", lwd = 0.55)
  
  points(temporal_data$year, temporal_data$genotype_position,
         pch = 21, bg = genotype_colors[as.character(temporal_data$genotype)],
         col = adjustcolor("black", alpha.f = 0.50),
         lwd = 0.45, cex = bubble_cex(temporal_data$n))
  
  rare_data <- temporal_data %>%
    filter(as.character(genotype) %in% c("HPeV18","HPeV19"))
  if (nrow(rare_data) > 0) {
    points(rare_data$year, rare_data$genotype_position,
           pch = 21, bg = genotype_colors[as.character(rare_data$genotype)],
           col = "black", lwd = 0.8,
           cex = pmax(bubble_cex(rare_data$n), 1.25))
  }
  
  axis(1, at = year_breaks, labels = year_breaks, las = 2, cex.axis = 0.78, tck = -0.012)
  mtext("Year", side = 1, line = 3.3, cex = 1)
  
  axis(2, at = 1:19, labels = paste0("PeV-A", 1:19),
       las = 1, cex.axis = 0.78, tick = FALSE)
  
  mtext("B", side = 3, adj = 0, line = -0.35, font = 2, cex = 1.45)
  
  possible_breaks <- c(1,5,10,25,50,100,250)
  size_breaks <- possible_breaks[possible_breaks <= max(temporal_data$n, na.rm = TRUE)]
  
  if (length(size_breaks) > 0) {
    legend_cex <- bubble_cex(size_breaks)
    legend(x = year_max + 4, y = 18.5, legend = size_breaks,
           pch = 21, pt.bg = "grey80", col = "grey35",
           pt.cex = legend_cex, pt.lwd = 0.7,
           title = "Sequences", cex = 0.88, bty = "n",
           xpd = NA, y.intersp = 1.45)
  }
}

#Full figure
draw_complete_figure <- function() {
  layout(matrix(c(1,2,3), nrow = 3), heights = c(4.85, 0.75, 3.65))
  draw_panel_A()
  draw_genotype_legend()
  draw_panel_B()
}

# Final check 

print(country_counts_complete %>%
        filter(as.character(genotype) %in% c("HPeV18","HPeV19"), n > 0) %>%
        select(country, genotype, n))
cat("\nTEMPORAL:\n")
print(temporal_data %>%
        filter(as.character(genotype) %in% c("HPeV18","HPeV19")) %>%
        select(year, genotype, n))
cat("\nLEGEND:\n")
print(genotype_labels[c("HPeV18","HPeV19")])

# Export PDF 
pdf("/Users/hugo/Desktop/vp1_parachovirus/Figure_1_PeV-A.pdf",
    width = 15, height = 8.5, useDingbats = FALSE, family = "Helvetica")
draw_complete_figure()
dev.off()
