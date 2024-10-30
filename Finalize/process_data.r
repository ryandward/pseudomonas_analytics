library(pacman)
p_load(data.table, tidyverse)

# Define output directory
output_dir <- "Finalize"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 1. Read and process condition mapping
mapping <- data.table(
  condition = c(
    "dJMP1", "dJMP2", "dJMP4",
    "Mouse_P1_003", "Mouse_P1_004", "Mouse_P1_006",
    "Mouse_P1_015", "Mouse_P1_016", "Mouse_P1_017",
    "Mouse_P1_018", "Mouse_P1_019", "Mouse_P1_020",
    "Mouse_P1_021", "Mouse_P1_022",
    "P1_mfdpir"
  ),
  specific_condition = c(
    "inoculum_rep2", "inoculum_rep3", "LB_6_generations_rep1",
    "PA14_transconjugant", "LB_6_generations_rep2", "inoculum_rep1",
    "high_inoculum_rep1", "high_inoculum_rep2", "high_inoculum_rep3",
    "low_inoculum_rep1", "low_inoculum_rep2", "low_inoculum_rep3",
    "low_inoculum_rep4", "low_inoculum_rep5",
    "E_coli_mating_strain"
  )
)
fwrite(mapping, file.path(output_dir, "condition_annotations.tsv"))

# 2. Process and filter guides
all_guides <- fread("all_guides.tsv") %>%
  filter(!sequence %in% c(
    "CCTTGATGCTGTTGAGGATC",
    "TGGTCTGGGTGCGCTCGGAG",
    "TAGTCAAAGATACCCCCGAA",
    "GCTCGCGGTTTACTTCGACC"
  )) %>%
  # Keep only non-focused guides
  filter(type != "focused")

fwrite(all_guides, file.path(output_dir, "guide_annotations.tsv"))

# 3. Process count stats with specific conditions
count_stats <- fread("stats.tsv", col.names = c("sequence", "count", "condition"))
count_stats_processed <- count_stats %>%
  left_join(mapping, by = "condition") %>%
  mutate(condition = specific_condition) %>%
  select(-specific_condition) %>%
  # Only keep sequences that are in our filtered guides
  inner_join(all_guides %>% select(sequence), by = "sequence")

fwrite(count_stats_processed, file.path(output_dir, "processed_counts.tsv"))

# 4. Create and save experimental design
exp_design <- count_stats_processed %>%
  select(condition) %>%
  unique() %>%
  mutate(
    group = case_when(
      condition %like% "inoculum" ~ "t0_inoculum",
      condition %like% "LB_6_generations" ~ "6_generations_LB",
      condition %like% "low_inoculum" ~ "10x_inoculum_dilution_mouse",
      condition %like% "high_inoculum" ~ "high_inoculum_mouse",
      TRUE ~ "other"
    )
  )

fwrite(exp_design, file.path(output_dir, "experimental_design.tsv"))

# Print summary of created files
cat("\nCreated processed data files in", output_dir, "folder:\n")
cat("1. condition_annotations.tsv - Mapping between conditions and specific conditions\n")
cat("2. guide_annotations.tsv - Filtered guide sequences (excludes focused guides)\n")
cat("3. processed_counts.tsv - Processed count data with specific conditions\n")
cat("4. experimental_design.tsv - Experimental design with groups\n")
