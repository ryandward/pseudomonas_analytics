library(pacman)

p_load(data.table, tidyverse, edgeR, poolr, scales, pheatmap, viridis)

guides <- fread("Finalize/guides.tsv")
guide_counts <- fread("Finalize/guide_counts.tsv")
guide_map_stats <- fread("Finalize/map_stats.tsv")
experimental_design <- fread("Finalize/experimental_design.tsv")

experimental_design <- experimental_design %>% inner_join(guide_map_stats)


experimental_design <- experimental_design %>%
  # remove counts from low inoculum dilution
  filter(!growth_condition %like% "100x_inoculum_dilution") %>%
  data.table()

# create sample group from gDNA_source, growth_condition, and media
experimental_design <- experimental_design %>%
  mutate(sample_group = paste(gDNA_source, growth_condition, media, sep = "_"))

experimental_design <- experimental_design %>%
  mutate(match_prop = matched / reads)

quality_decision <- guide_counts %>%
  inner_join(guides) %>%
  pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0) %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  group_by(condition) %>%
  mutate(cpm = 1e6 * (count / sum(count))) %>%
  group_by(condition, type) %>%
  summarise(cpm = sum(cpm)) %>%
  pivot_wider(id_cols = condition, names_from = type, values_from = cpm) %>%
  inner_join(experimental_design) %>%
  select(condition, control, knockdown, media, gDNA_source, match_prop, reads)

# Filter out guides that don't have enough reads or match proportion
experimental_design <- experimental_design %>%
  inner_join(quality_decision) %>%
  filter(match_prop > 0.85 & reads > 5e6) %>%
  group_by(sample_group) %>%
  mutate(rep = as.character(seq_len(n())))

# Remove things that don't have replicates
experimental_design <- experimental_design %>%
  filter(sample_group %in% (experimental_design %>%
    group_by(
      sample_group
    ) %>%
    summarise(reps = max(rep)) %>%
    filter(reps > 1) %>%
    pull(sample_group)))

experimental_design <- experimental_design %>% data.table()

################################################################################

guide_counts <- guide_counts %>% filter(condition %in% experimental_design$condition)

setorder(experimental_design, condition)
setorder(guide_counts, condition)

guide_counts.matrix <- guide_counts %>%
  inner_join(guides) %>%
  pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0)

guide_counts <- guide_counts.matrix %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  inner_join(guides) %>%
  unite("guide", c(gene, offset), remove = FALSE)

##########################################################################################

guide_counts.matrix.quality <-
  guide_counts.matrix %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  inner_join(experimental_design) %>%
  filter(type != "focused") %>%
  pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0)

##########################################################################################

guide_counts.quality <-
  guide_counts.matrix.quality %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  inner_join(guides) %>%
  unite("guide", c(gene, offset), remove = FALSE) %>%
  mutate(guide = case_when(type == "control" ~ paste(gene, index), type == "knockdown" ~ paste(gene, offset)))

setorder(guide_counts.quality, condition)
setorder(experimental_design, condition)

conditions.quality <-
  guide_counts.quality %>%
  select(condition) %>%
  unique()

##########################################################################################

data_grid <- guide_counts.quality %>%
  pivot_wider(id_cols = sequence, names_from = condition, values_from = count)

data_grid_matrix <- data_grid %>%
  select(-sequence) %>%
  data.matrix()

row.names(data_grid_matrix) <- data_grid$sequence

colnames(data_grid_matrix) <- data_grid %>%
  colnames() %>%
  data.table(condition = .) %>%
  inner_join(experimental_design) %>%
  pull(condition)

data_group <-
  factor(
    experimental_design[condition %in% conditions.quality$condition, sample_group],
    levels = unique(unique(experimental_design[condition %in% conditions.quality$condition, sample_group]))
  )

data_permut <- model.matrix(~ 0 + data_group)

colnames(data_permut) <- levels(data_group)

rownames(data_permut) <- colnames(data_grid_matrix)

data_permut_check <-
  melt(
    data.table(
      data_permut,
      keep.rownames = "condition"
    ),
    id.vars = "condition"
  ) %>%
  filter(value == 1) %>%
  select(condition) %>%
  inner_join(experimental_design)

print(data_permut_check)

##########################################################################################

data_y <- DGEList(
  counts = data_grid_matrix,
  group = data_group,
  genes = row.names(data_grid_matrix)
)

data_keep <- filterByExpr(data_grid_matrix, data_group)

data_y <- data_y[data_keep, , keep.lib.sizes = FALSE]

data_y <- calcNormFactors(data_y)

data_y <- estimateDisp(data_y, data_permut)

data_fit <- glmQLFit(data_y, data_permut, robust = TRUE)

data_CPM <- cpm(data_y, prior.count = 0)

################################################################################

contrast_levels <-
  c(
    "plated_6_generations_LB - plated_t0_inoculum",
    "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum",
    "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB"
  )

data_contrast <- makeContrasts(
  contrasts = contrast_levels,
  levels = data_permut
)

################################################################################

results_FDR <- data.table(guide_counts.quality)[, .(genes = unique(sequence))]

results_LFC <- data.table(guide_counts.quality)[, .(genes = unique(sequence))]

################################################################################

for (i in seq_len(ncol(data_contrast))) {
  results <- glmQLFTest(data_fit, contrast = data_contrast[, i])

  results <- topTags(results, n = Inf)

  results <- data.table(results$table)

  print(paste("Processing results for", contrast_levels[i], "..."))

  results_FDR <- results[, .(genes, FDR)][results_FDR, on = .(genes)]

  setnames(results_FDR, "FDR", contrast_levels[i])

  results_LFC <-
    results[, .(genes, logFC)][results_LFC, on = .(genes)]

  setnames(results_LFC, "logFC", contrast_levels[i])
}

################################################################################

results_list_FDR <-
  data.table::melt(
    results_FDR,
    id.vars = c("genes"),
    variable.name = "condition",
    value.name = "FDR",
    measure.vars = contrast_levels
  )

results_list_FDR <- results_list_FDR[!is.na(FDR)]

################################################################################

results_list_LFC <-
  data.table::melt(
    results_LFC,
    id.vars = c("genes"),
    variable.name = "condition",
    value.name = "LFC",
    measure.vars = contrast_levels
  )

results_list_LFC[, LFC := round(LFC, 3)]

################################################################################

results_list <-
  results_list_LFC[
    results_list_FDR,
    on = .(genes, condition)
  ] %>%
  rename(sequence = genes) %>%
  inner_join(guides)

results_list <- results_list[!is.na(FDR)]

results_list_by_condition <-
  results_list[
    type == "control",
    .(med_LFC = median(LFC)),
    keyby = .(condition)
  ]

setkey(results_list, condition)

results_list[, LFC := results_list_by_condition[
  results_list, LFC - med_LFC,
  by = .EACHI
]$V1]

################################################################################

median_results_list <-
  results_list[, .(
    medLFC = median(LFC),
    FDR = stouffer(FDR)$p
  ),
  by = .(locus_tag, gene, condition)
  ]

#################################################################################

median_results_list[gene != ".", gene_name_stylized := paste0("italic('", gene, "')")]
median_results_list[gene == ".", gene_name_stylized := paste0("bold('", locus_tag, "')")]
median_results_list[gene == "", gene_name_stylized := paste0("bold('", locus_tag, "')")]
median_results_list[gene == "control", gene_name_stylized := paste0("bold('", locus_tag, "')")]
