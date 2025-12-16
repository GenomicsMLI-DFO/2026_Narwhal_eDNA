
############################################################
# Project: 2026_Narwhal_eDNA
# Script: Main Analysis
# Author: Geneviève J. Parent
# Date: 2025-12-16
# Description: Data analyses for Narwhal eDNA project
############################################################

# -----------------------------
# 1. Setup Environment
# -----------------------------
# Clear workspace
rm(list = ls())

# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(knitr)
library(seqinr)

# -----------------------------
# 2. Load Data
# -----------------------------
data_160 <- read_csv("Data/ESVtab.corrected_12S160_Narval_Isla_Duporge.csv")
data_MFU <- read_csv("Data/ESVtab.corrected_MiFishU_Narval_Isla_Duporge.csv")
data_taxo <- read_csv("Data/Narval_Isla_Duporge_ESVtab_taxo_report_nt_ONLY.csv")
MOTU_160 <- read_csv("Data/MOTUs.Metabarinfo.corrected_12S160_Narval_Isla_Duporge.csv")
MOTU_MFU <- read_csv("Data/MOTUs.Metabarinfo.corrected_MiFishU_Narval_Isla_Duporge.csv")

# -----------------------------
# 3. Data Exploration
# -----------------------------
# Quick checks
list(data_160 = dim(data_160), data_MFU = dim(data_MFU), data_taxo = dim(data_taxo))
anyNA(data_160); anyNA(data_MFU); anyNA(data_taxo)

# -----------------------------
# 4. Basic Analyses
# -----------------------------
# Function to compute row sums and summary stats
compute_summary <- function(df) {
  df <- df %>%
    rowwise() %>%
    mutate(row_sum = sum(c_across(starts_with("ESV_")), na.rm = TRUE)) %>%
    ungroup()
  
  stats <- list(
    mean = mean(df$row_sum, na.rm = TRUE),
    sd = sd(df$row_sum, na.rm = TRUE),
    total = sum(df$row_sum, na.rm = TRUE)
  )
  
  return(list(data = df, stats = stats))
}

# Apply to both datasets
summary_160 <- compute_summary(data_160)
summary_MFU <- compute_summary(data_MFU)

# Print summary stats
summary_160$stats
summary_MFU$stats

# -----------------------------
# 5. Species-Level Table
# -----------------------------
# Filter taxa: remove NA species, zero counts, and unwanted species
data_taxo <- data_taxo %>%
  filter(!is.na(species), Nreads > 0, species != "Meleagris gallopavo")

# Compute species-level summary per locus
species_reads_table <- data_taxo %>%
  group_by(Loci, species) %>%
  summarise(
    detections = n(),
    total_reads = sum(Nreads, na.rm = TRUE),
    .groups = "drop"
  )

# Compute total reads per locus
total_reads_per_locus <- species_reads_table %>%
  group_by(Loci) %>%
  summarise(locus_total = sum(total_reads), .groups = "drop")

# Join and calculate proportions
species_reads_table <- species_reads_table %>%
  left_join(total_reads_per_locus, by = "Loci") %>%
  mutate(proportion = (total_reads / locus_total) * 100)

# Pivot wider for publication
species_reads_table <- species_reads_table %>%
  pivot_wider(
    id_cols = species,
    names_from = Loci,
    values_from = c(detections, total_reads, proportion),
    values_fill = 0
  )

# Display table
kable(
  species_reads_table,
  caption = "Species Detected per Locus: Counts, Reads, and Proportions",
  align = "lcccccc",
  format = "markdown"
)

# -----------------------------
# 6. Export Table
# -----------------------------

# Create Results folder
if (!dir.exists("Results")) {
  dir.create("Results")
}
write.csv(species_reads_table, file = "Results/species_summary_table_1.csv", row.names = FALSE)

# -----------------------------
# 7. Extract ESVs FASTA for narwhal and unexpected marine mammal species
# -----------------------------
species_list <- c("Monodon monoceros", "Cephalorhynchus heavisidii", "Sagmatias obliquidens")

selected_df <- MOTU_160 %>%
  filter(species_NCBI %in% species_list) %>%
  select(sequence, species_NCBI, ESV) %>%
  distinct(species_NCBI, sequence, .keep_all = TRUE) %>%
  mutate(label = paste0(species_NCBI))

fasta_file <- file.path("Results","selected_species_160.fasta")
write.fasta(sequences = as.list(selected_df$sequence),
            names = gsub(" ", "_", selected_df$label),
            file.out = fasta_file)

selected_df <- MOTU_MFU %>%
  filter(species_NCBI %in% species_list) %>%
  select(sequence, species_NCBI, ESV) %>%
  distinct(species_NCBI, sequence, .keep_all = TRUE) %>%
  mutate(label = paste0(species_NCBI))

fasta_file <- file.path("Results", "selected_species_MFU.fasta")
write.fasta(sequences = as.list(selected_df$sequence),
            names = gsub(" ", "_", selected_df$label),
            file.out = fasta_file)

# -----------------------------
# 8. Session Info
# -----------------------------
sessionInfo()
