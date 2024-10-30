library(pacman)

p_load(data.table, tidyverse, edgeR, poolr, scales, pheatmap, viridis)

count_stats <- fread("stats.tsv", col.names = c("sequence", "count", "condition"))

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

all_guides <- fread("all_guides.tsv") %>%
  filter(!sequence %in% c(
    "CCTTGATGCTGTTGAGGATC",
    "TGGTCTGGGTGCGCTCGGAG",
    "TAGTCAAAGATACCCCCGAA",
    "GCTCGCGGTTTACTTCGACC"
  )) %>%
  group_by(type) %>%
  mutate(
    offset = ifelse(
      is.na(offset),
      NA_real_,
      offset
    ),
    offset = case_when(
      is.na(offset) ~ 0,
      !is.na(offset) ~ offset
    ),
    index = seq_len(n()),
    locus_tag = case_when(
      type == "control" ~ "control",
      type == "knockdown" ~ locus_tag,
      type == "focused" ~ locus_tag
    ),
    gene = case_when(
      type == "control" ~ "control",
      gene == "-" ~ locus_tag,
      type == "knockdown" ~ gene,
      type == "focused" ~ gene
    ),
    guide = case_when(
      type == "control" ~ paste(locus_tag, index),
      type == "focused" ~ paste(locus_tag, index),
      type == "knockdown" ~ paste(locus_tag, offset)
    ),
    target_genome = case_when(
      type == "focused" ~ "focused",
      type == "control" ~ "control",
      type == "knockdown" ~ "essential"
    )
  )

guide_definitions <- fread("guide_definitions.tsv")
# guide names from Neha

first_annotations <- fread("first_annotations.tsv")
# my first attempt to annotate comprehensively

synonyms <- fread("synonyms.tsv", fill = TRUE, na.strings = c(""))

y_genes_neha_found <-
  all_guides %>%
  left_join(guide_definitions) %>%
  left_join(
    first_annotations %>%
      select(gene_name, guide_name)
  ) %>%
  filter(gene != gene_name) %>%
  select(locus_tag, gene, gene_name, locus_tag) %>%
  left_join(synonyms) %>%
  filter(gene %like% "^y") %>%
  unique()

# use the names that Neha found for y-genes
all_guides <- all_guides %>%
  left_join(y_genes_neha_found) %>%
  mutate(gene = case_when(
    !is.na(gene_name) ~ gene_name,
    TRUE ~ gene
  )) %>%
  select(-gene_name, -syn1, -syn2)

all_guides <- all_guides %>%
  mutate(gene = case_when(
    locus_tag == "PA14_07770" ~ "lptD",
    locus_tag == "PA14_23560" ~ "gltX",
    locus_tag == "PA14_12280" ~ "lnt",
    locus_tag == "PA14_14880" ~ "ispG",
    locus_tag == "PA14_17130" ~ "ispC",
    locus_tag == "PA14_31290" ~ "lecA",
    locus_tag == "PA14_32420" ~ "mexS",
    locus_tag == "PA14_33530" ~ "fpvF",
    locus_tag == "PA14_43970" ~ "lpd",
    locus_tag == "PA14_62870" ~ "rlmE",
    locus_tag == "PA14_72970" ~ "tonB1",
    locus_tag == "PA14_65960" ~ "waaA",
    locus_tag == "PA14_66060" ~ "waaE",
    locus_tag == "PA14_66220" ~ "waaP",
    TRUE ~ gene
  ))

# Guide distribution during construction
count_stats %>%
  filter(condition %in% c("P1_mfdpir", "Mouse_P1_003")) %>%
  data.table() %>%
  dcast(sequence ~ condition, value.var = "count") %>%
  right_join(all_guides %>% select(sequence, locus_tag, type) %>% filter(type != "focused")) %>%
  melt(id.vars = c("sequence", "type", "locus_tag"), variable.name = "condition", value.name = "count") %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(condition) %>%
  mutate(cpm = count * 1e6 / sum(count)) %>%
  ggplot(aes(x = cpm, color = type)) +
  geom_density(lwd = 2) +
  facet_grid(condition ~ .) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10^(0:6)),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  theme_minimal()

map_stats <- fread("map_stats.tsv")

exp_design <- fread("unified_counts/exp_design.tsv")

exp_design <- exp_design %>% inner_join(map_stats)

exp_design <- exp_design %>%
  # this one didn't get infected
  filter(!growth_condition %like% "100x_inoculum_dilution") %>%
  filter(!condition %in% c("dJMP3", "dJMP5")) %>%
  # filter(!condition %in% c("dJMP1")) %>%

  data.table()

exp_design <- exp_design %>% mutate(sample_group = paste(gDNA_source, growth_condition, media, sep = "_"))

exp_design <- exp_design %>%
  mutate(match_prop = matched / reads)

quality_decision <- count_stats %>%
  inner_join(all_guides) %>%
  pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0) %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  group_by(condition) %>%
  mutate(cpm = 1e6 * (count / sum(count))) %>%
  group_by(condition, type) %>%
  summarise(cpm = sum(cpm)) %>%
  pivot_wider(id_cols = condition, names_from = type, values_from = cpm) %>%
  arrange(desc(focused)) %>%
  inner_join(exp_design) %>%
  select(condition, control, focused, knockdown, media, gDNA_source, match_prop, reads)

exp_design <- exp_design %>%
  inner_join(quality_decision) %>%
  filter(match_prop > 0.85 & reads > 5e6) %>%
  group_by(sample_group) %>%
  mutate(rep = as.character(seq_len(n())))

exp_design <- exp_design %>% mutate(sample_name = paste(gDNA_source, growth_condition, media, rep, sep = "_"))

exp_design <- exp_design %>%
  filter(sample_group %in% (exp_design %>%
    group_by(
      sample_group
    ) %>%
    summarise(reps = max(rep)) %>%
    filter(reps > 1) %>%
    pull(sample_group)))

exp_design <- exp_design %>% data.table()

################################################################################

count_stats %>%
  inner_join(all_guides) %>%
  filter(type != "focused") %>%
  group_by(sequence, condition) %>%
  inner_join(exp_design) %>%
  pivot_wider(id_cols = sequence, names_from = sample_name, values_fill = 0, values_from = count) %>%
  ungroup() %>%
  select(-sequence) %>%
  data.matrix() %>%
  cor() %>%
  round(4) %>%
  pheatmap(
    cutree_rows = 4,
    cutree_cols = 4,
    display_numbers = TRUE,
    number_format = "%.3f",
    breaks = seq(0.5, 1, length.out = 9999),
    color = c(colorRampPalette(c("navy", "white", "red"))(9997), "gray")
  )

################################################################################

count_stats <- count_stats %>% filter(condition %in% exp_design$condition)

setorder(exp_design, condition)
setorder(count_stats, condition)


count_stats.mat <- count_stats %>%
  inner_join(all_guides) %>%
  pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0)

count_stats.mat %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  inner_join(all_guides)

count_stats <- count_stats.mat %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  inner_join(all_guides) %>%
  unite("guide", c(gene, offset), remove = FALSE)

##########################################################################################
# which kind of gDNA is better? plates.
plot_object <- count_stats.mat %>%
  select(-type, -sequence) %>%
  data.matrix() %>%
  cpm() %>%
  as_tibble() %>%
  cbind(
    count_stats.mat %>%
      select(type, sequence)
  ) %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "cpm") %>%
  inner_join(exp_design) %>%
  unite("growth_condition", c(growth_condition, media)) %>%
  ggplot(aes(x = cpm)) +
  geom_density(aes(fill = rep), alpha = 0.25) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10^(1:6)),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  facet_grid(facets = c("gDNA_source", "growth_condition")) +
  ggtitle("Counts per Million")

print(plot_object)

##########################################################################################
# is there a difference in the composition of the knockdown and control components
# only use plated non-focused

count_stats.mat.quality <-
  count_stats.mat %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  inner_join(exp_design) %>%
  filter(type != "focused") %>%
  pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0)

plot_object <- count_stats.mat.quality %>%
  select(-type, -sequence) %>%
  data.matrix() %>%
  cpm() %>%
  as_tibble() %>%
  cbind(
    count_stats.mat.quality %>%
      select(type, sequence)
  ) %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "cpm") %>%
  inner_join(exp_design) %>%
  ggplot(aes(x = cpm)) +
  geom_density(aes(fill = rep), alpha = 0.35) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10^(1:6)),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  facet_grid(facets = c("type", "media")) +
  ggtitle("Counts per Million") +
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  theme_bw() +
  ggtitle("Distribution of Guides Recovered (Counts per Million)")

print(plot_object)

##########################################################################################

count_stats.quality <-
  count_stats.mat.quality %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  inner_join(all_guides) %>%
  unite("guide", c(gene, offset), remove = FALSE) %>%
  mutate(guide = case_when(type == "control" ~ paste(gene, index), type == "knockdown" ~ paste(gene, offset)))

setorder(count_stats.quality, condition)
setorder(exp_design, condition)

conditions.quality <-
  count_stats.quality %>%
  select(condition) %>%
  unique()

##########################################################################################

data_grid <- count_stats.quality %>%
  pivot_wider(id_cols = sequence, names_from = condition, values_from = count)

data_grid_matrix <- data_grid %>%
  select(-sequence) %>%
  data.matrix()

row.names(data_grid_matrix) <- data_grid$sequence

colnames(data_grid_matrix) <- data_grid %>%
  colnames() %>%
  data.table(condition = .) %>%
  inner_join(exp_design) %>%
  pull(sample_name)

data_group <-
  factor(
    exp_design[condition %in% conditions.quality$condition, sample_group],
    levels = unique(unique(exp_design[condition %in% conditions.quality$condition, sample_group]))
  )

data_permut <- model.matrix(~ 0 + data_group)

colnames(data_permut) <- levels(data_group)

rownames(data_permut) <- colnames(data_grid_matrix)

data_permut_check <-
  melt(
    data.table(
      data_permut,
      keep.rownames = "sample_name"
    ),
    id.vars = "sample_name"
  ) %>%
  filter(value == 1) %>%
  select(sample_name) %>%
  inner_join(exp_design)

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

results_FDR <- data.table(count_stats.quality)[, .(genes = unique(sequence))]

results_LFC <- data.table(count_stats.quality)[, .(genes = unique(sequence))]

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

melted_results_FDR <-
  data.table::melt(
    results_FDR,
    id.vars = c("genes"),
    variable.name = "condition",
    value.name = "FDR",
    measure.vars = contrast_levels
  )

melted_results_FDR <- melted_results_FDR[!is.na(FDR)]

################################################################################

melted_results_LFC <-
  data.table::melt(
    results_LFC,
    id.vars = c("genes"),
    variable.name = "condition",
    value.name = "LFC",
    measure.vars = contrast_levels
  )

melted_results_LFC[, LFC := round(LFC, 3)]

################################################################################

melted_results <-
  melted_results_LFC[
    melted_results_FDR,
    on = .(genes, condition)
  ] %>%
  rename(sequence = genes) %>%
  inner_join(all_guides)

melted_results <- melted_results[!is.na(FDR)]

melted_results_by_condition <-
  melted_results[
    type == "control",
    .(med_LFC = median(LFC)),
    keyby = .(condition)
  ]

setkey(melted_results, condition)

melted_results[, LFC := melted_results_by_condition[
  melted_results, LFC - med_LFC,
  by = .EACHI
]$V1]

################################################################################

median_melted_results <-
  melted_results[, .(
    medLFC = median(LFC),
    FDR = stouffer(FDR)$p
  ),
  by = .(locus_tag, gene, condition)
  ]

#################################################################################

median_melted_results[gene != ".", gene_name_stylized := paste0("italic('", gene, "')")]
median_melted_results[gene == ".", gene_name_stylized := paste0("bold('", locus_tag, "')")]
median_melted_results[gene == "", gene_name_stylized := paste0("bold('", locus_tag, "')")]
median_melted_results[gene == "control", gene_name_stylized := paste0("bold('", locus_tag, "')")]

################################################################################

plot_object <- melted_results %>%
  inner_join(all_guides) %>%
  mutate(condition = case_when(
    condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = LFC)) +
  geom_density(aes(fill = type), alpha = 0.35) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10^(1:6)),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  facet_wrap(facets = c("condition")) +
  ggtitle("Counts per Million") +
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  theme_bw() +
  ggtitle("Knockdown-induced Guide Composition Change (log-fold change)")

print(plot_object)

################################################################################

plot_object <- median_melted_results %>%
  mutate(condition = case_when(
    condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = medLFC)) +
  geom_density(aes(fill = condition), alpha = 0.35) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10^(1:6)),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  ggtitle("Counts per Million") +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  theme_bw() +
  ggtitle("Knockdown-induced Gene Composition Change (log-fold change)")

print(plot_object)

################################################################################
# For pasting results into excel

# median_melted_results %>%
#   mutate(condition = case_when(
#     condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
#     condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
#     condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
#   )) %>%
#   pivot_wider(id_cols = c(locus_tag, gene), values_from = medLFC, names_from = condition) %>%
#   clipr::write_clip()

# median_melted_results %>%
#   mutate(condition = case_when(
#     condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
#     condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
#     condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
#   )) %>%
#   pivot_wider(id_cols = c(locus_tag, gene), values_from = FDR, names_from = condition) %>%
#   clipr::write_clip()
