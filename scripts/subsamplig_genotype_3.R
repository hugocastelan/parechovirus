# HPeV3 - Stratified Downsampling by Region and Year
# Maximum 10 sequences per region/year

rm(list = ls())

library(Biostrings)
library(dplyr)
library(stringr)
library(readr)

# 1. Files
setwd("/Users/hugo/Downloads/")

fasta_file <- "HPeV3_round3.fasta"

output_fasta  <- "/Users/hugo/HPeV3_downsampled_region_year.fasta"
output_traits <- "/Users/hugo/HPeV3_downsampled_region_traits.tsv"
output_dates  <- "/Users/hugo/HPeV3_downsampled_dates.csv"
output_meta   <- "/Users/hugo/HPeV3_downsampled_metadata.csv"

# 2. Read FASTA
seqs <- readDNAStringSet(fasta_file)
taxa <- names(seqs)
cat("Original sequences:", length(seqs), "\n")

# 3. Extract metadata from FASTA headers
# Format: accession|virus|genotype|country|year
tmp <- str_split_fixed(taxa, "\\|", 5)

metadata <- data.frame(
  taxon     = taxa,
  accession = tmp[, 1],
  virus     = tmp[, 2],
  genotype  = tmp[, 3],
  country   = tmp[, 4],
  year      = tmp[, 5],
  stringsAsFactors = FALSE
)

metadata$year_numeric <- suppressWarnings(as.numeric(metadata$year))

# 4. Check dates
if (any(is.na(metadata$year_numeric))) {
  cat("\nERROR: taxa with invalid dates:\n")
  print(metadata %>% filter(is.na(year_numeric)) %>% select(taxon, country, year))
  stop("Invalid or missing dates detected.")
}

# 5. Assign geographic regions
metadata <- metadata %>%
  mutate(
    region = case_when(
      # Europe
      country %in% c("Netherlands", "Norway", "Finland", "Estonia", "Hungary", 
                     "Greece", "Ireland", "Germany", "Belgium", "Denmark", 
                     "Bulgaria", "Czech", "Czech_Republic", "France", "Italy", 
                     "Spain", "Portugal", "Sweden", "United_Kingdom", "UK", 
                     "England", "Scotland") ~ "Europe",
      
      # East Asia
      country %in% c("Japan", "China", "Hong_Kong", "South_Korea", "Korea", "Taiwan") ~ "East_Asia",
      
      # Southeast Asia
      country %in% c("Vietnam", "Thailand", "Singapore", "Malaysia", "Indonesia", 
                     "Philippines", "Cambodia", "Myanmar") ~ "Southeast_Asia",
      
      # South Asia
      country %in% c("India", "Pakistan", "Sri_Lanka", "Bangladesh", "Nepal") ~ "South_Asia",
      
      # North America
      country %in% c("USA", "United_States", "Canada", "Mexico") ~ "North_America",
      
      # South America
      country %in% c("Bolivia", "Brazil", "Venezuela", "Argentina", "Chile", 
                     "Peru", "Colombia", "Ecuador") ~ "South_America",
      
      # Africa
      country %in% c("Malawi", "Ghana", "Nigeria", "South_Africa", "Kenya", 
                     "Cameroon", "Ivoire_Coast", "Cote_d'Ivoire", "Cote_d_Ivoire", 
                     "Ivory_Coast", "Senegal", "Tanzania", "Uganda", "Ethiopia") ~ "Africa",
      
      # Oceania
      country %in% c("Australia", "New_Zealand") ~ "Oceania",
      
      TRUE ~ "Unknown"
    )
  )

# 6. Check for unknown regions
unknown <- metadata %>% filter(region == "Unknown") %>% count(country, sort = TRUE)
if (nrow(unknown) > 0) {
  cat("\nCOUNTRIES WITHOUT REGION:\n")
  print(as.data.frame(unknown), row.names = FALSE)
  stop("Unknown countries found. Add them to region mapping.")
} else {
  cat("\nAll countries assigned to a region.\n")
}

# 7. Original distribution by region
cat("\nORIGINAL DATASET BY REGION:\n")
metadata %>% count(region, name = "n") %>% arrange(desc(n)) %>% as.data.frame() %>% print(row.names = FALSE)

# 8. Original distribution region × year
cat("\nORIGINAL DATASET BY REGION × YEAR:\n")
metadata %>% count(region, year_numeric, name = "n") %>% arrange(region, year_numeric) %>% as.data.frame() %>% print(row.names = FALSE)

# 9. Stratified downsampling (max 10 per region × year)
set.seed(12345)

metadata_down <- metadata %>%
  group_by(region, year_numeric) %>%
  group_modify(~ {
    if (nrow(.x) > 10) {
      slice_sample(.x, n = 10)
    } else {
      .x
    }
  }) %>%
  ungroup()

# 10. Results summary
cat("\nDOWNSAMPLING RESULTS:\n")
cat("Original sequences:   ", nrow(metadata), "\n")
cat("Downsampled sequences:", nrow(metadata_down), "\n")
cat("Removed sequences:    ", nrow(metadata) - nrow(metadata_down), "\n")

# 11. Final distribution by region
cat("\nFINAL DATASET BY REGION:\n")
metadata_down %>% count(region, name = "n") %>% arrange(desc(n)) %>% as.data.frame() %>% print(row.names = FALSE)

# 12. Final distribution region × year
cat("\nFINAL DATASET BY REGION × YEAR:\n")
metadata_down %>% count(region, year_numeric, name = "n") %>% arrange(region, year_numeric) %>% as.data.frame() %>% print(row.names = FALSE)

# 13. Verify max = 10
check_max <- metadata_down %>% count(region, year_numeric)
cat("\nMaximum sequences in any region × year:", max(check_max$n), "\n")
if (max(check_max$n) > 10) {
  stop("ERROR: At least one region × year has >10 sequences.")
}

# 14. Extract selected FASTA sequences
selected_positions <- match(metadata_down$taxon, names(seqs))
if (any(is.na(selected_positions))) {
  stop("ERROR: selected taxa not found in FASTA.")
}
seqs_down <- seqs[selected_positions]

# 15. Verify order
stopifnot(all(names(seqs_down) == metadata_down$taxon))

# 16. Write FASTA
writeXStringSet(seqs_down, filepath = output_fasta, format = "fasta")

# 17. Write location traits (taxon<TAB>region)
traits <- metadata_down %>% select(taxon, region)
write.table(traits, file = output_traits, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 18. Write dates
dates <- metadata_down %>% transmute(name = taxon, date = year_numeric)
write.csv(dates, output_dates, row.names = FALSE, quote = FALSE)

# 19. Write complete metadata
write.csv(metadata_down, output_meta, row.names = FALSE, quote = FALSE)

# 20. Final validation
cat("\nFINAL VALIDATION:\n")
cat("FASTA sequences:", length(seqs_down), "\n")
cat("Metadata rows:  ", nrow(metadata_down), "\n")
cat("Trait rows:     ", nrow(traits), "\n")
cat("Date rows:      ", nrow(dates), "\n")

stopifnot(
  length(seqs_down) == nrow(metadata_down),
  nrow(metadata_down) == nrow(traits),
  nrow(traits) == nrow(dates),
  all(names(seqs_down) == metadata_down$taxon),
  max(check_max$n) <= 10
)

cat("\nEVERYTHING IS CORRECT.\n")

# 21. Output files
cat("\nGenerated files:\n")
cat(output_fasta, "\n")
cat(output_traits, "\n")
cat(output_dates, "\n")
cat(output_meta, "\n")
