library(pacman)

p_load(data.table, tidyverse, edgeR, poolr, scales, pheatmap, viridis)

# count_stats <- fread("stats.tsv", col.names = c("sequence", "count", "condition"))

count_stats <- fread("Combined_Sequencing/validated_counts.tsv.gz", col.names = c("condition", "spacer", "count"))

# count_stats <- count_stats %>% mutate(condition = str_replace_all(condition, "\\.counts\\.tsv", ""), sequence = str_replace_all(sequence, "\\*", ""))

all_guides <- fread("annotated_guides_from_neha.tsv")

designed_targets <- fread("designed_targets.tsv") %>%
mutate(type = case_when(
  mismatches == 0 ~ "knockdown",
  mismatches == "None" ~ "control"
))

count_stats <- count_stats %>%
  data.table %>%
  dcast(spacer ~ condition, value.var = "count", fill = 0) %>%
  melt(id.vars = c("spacer"), variable.name = "sample", value.name = "count") %>%
  inner_join(designed_targets %>% select(spacer, name, mismatches) %>% unique())

exp_design <- fread("unified_counts/exp_design.tsv") 

exp_design <- exp_design %>% rename(sample = condition)

exp_design <- exp_design %>%
  filter(gDNA_source == "plated") %>%
  # this one didn't get infected
  filter(!growth_condition %like% "100x_inoculum_dilution") %>%
  filter(!sample %in% c("dJMP3", "dJMP5")) %>%
  filter(!sample %in% c("dJMP1")) %>%
  # filter(!sample %in% c("Mouse_P1_006", "Mouse_P1_004"))

  data.table()

exp_design <- exp_design %>% mutate(sample_group = paste(gDNA_source, growth_condition, media, sep = "_"))

exp_design <-  exp_design %>% group_by(sample_group) %>% mutate(rep = as.character(seq_len(n())))

exp_design <- exp_design %>% mutate(sample_name = paste(gDNA_source, growth_condition, media, rep, sep = "_"))

exp_design <- exp_design %>% filter(sample_group %in% (exp_design %>%
  group_by(
    sample_group
  ) %>%
  summarise(reps = max(rep)) %>%
  filter(reps > 1) %>%
  pull(sample_group)))


exp_design <- exp_design %>% data.table()

################################################################################

count_stats %>%
  select(spacer, sample, count) %>% unique %>%
  inner_join(designed_targets %>% select(spacer, mismatches) %>% unique()) %>%
  inner_join(exp_design) %>%
  group_by(spacer, sample) %>%
  mutate(identifier = paste(sample, sample_name)) %>%
  pivot_wider(id_cols = spacer, names_from = identifier, values_fill = 0, values_from = count) %>%
  ungroup() %>%
  select(-spacer) %>%
  data.matrix() %>%
  cor() %>%
  round(4) %>%
  pheatmap(
    cutree_rows = 3,
    cutree_cols = 3,
    display_numbers = TRUE,
    number_format = "%.3f",
    breaks = seq(0.5, 1, length.out = 9999),
    color = c(colorRampPalette(c("navy", "white", "red"))(9997), "gray")
  )

################################################################################

count_stats <- count_stats %>% filter(sample %in% exp_design$sample)

setorder(exp_design, sample)
setorder(count_stats, sample)


count_stats.mat <- count_stats %>%
  inner_join(designed_targets) %>%
  pivot_wider(id_cols = c(type, spacer), names_from = sample, values_from = count, values_fill = 0)

# count_stats.mat %>%
#   pivot_longer(cols = -c(type, spacer), names_to = "sample", values_to = "count") %>%
#   inner_join(designed_targets)

count_stats <- count_stats.mat %>%
  pivot_longer(cols = -c(type, spacer), names_to = "sample", values_to = "count") %>%
  inner_join(designed_targets) %>%
  unite("guide", c(locus_tag, offset), remove = FALSE)

##########################################################################################
# which kind of gDNA is better? plates.
plot_object <- count_stats.mat %>%
  select(-type, -spacer) %>%
  data.matrix() %>%
  cpm() %>%
  as_tibble() %>%
  cbind(
    count_stats.mat %>%
      select(type, spacer)
  ) %>%
  pivot_longer(cols = -c(type, spacer), names_to = "sample", values_to = "cpm") %>%
  inner_join(exp_design) %>%
  unite("growth_condition", c(growth_condition, media)) %>%
  ggplot(aes(x = cpm)) +
  geom_density(aes(fill = rep), alpha = 0.25) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10^(1:6)),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  facet_grid(type ~ growth_condition) +
  ggtitle("Counts per Million")

print(plot_object)


##########################################################################################


data_grid <- count_stats %>%
  pivot_wider(id_cols = spacer, names_from = sample, values_from = count)

data_grid_matrix <- data_grid %>%
  select(-spacer) %>%
  data.matrix()

row.names(data_grid_matrix) <- data_grid$spacer

colnames(data_grid_matrix) <- data_grid %>%
  colnames() %>%
  data.table(sample = .) %>%
  inner_join(exp_design) %>%
  pull(sample_name)

data_group <-
  factor(
    exp_design[sample %in% exp_design$sample, sample_group],
    levels = unique(unique(exp_design[sample %in% exp_design$sample, sample_group]))
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

# data_keep <- filterByExpr(data_grid_matrix, data_group)

# data_y <- data_y[data_keep, , keep.lib.sizes = TRUE]

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

results_FDR <- data.table(count_stats)[, .(genes = unique(spacer))]

results_LFC <- data.table(count_stats)[, .(genes = unique(spacer))]

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
    variable.name = "sample",
    value.name = "FDR",
    measure.vars = contrast_levels
  )

melted_results_FDR <- melted_results_FDR[!is.na(FDR)]

################################################################################

melted_results_LFC <-
  data.table::melt(
    results_LFC,
    id.vars = c("genes"),
    variable.name = "sample",
    value.name = "LFC",
    measure.vars = contrast_levels
  )

melted_results_LFC[, LFC := round(LFC, 3)]

################################################################################

melted_results <- melted_results_LFC %>% 
inner_join(melted_results_FDR, by = c("genes", "sample")) %>% 
rename(spacer = genes) %>% 
inner_join(designed_targets) %>% 
arrange(LFC) %>% inner_join(all_guides %>% select(spacer, gene))

melted_results <- melted_results[!is.na(FDR)]

melted_results_by_sample <-
  melted_results[
    type == "control",
    .(med_LFC = median(LFC)),
    keyby = .(sample)
  ]

setkey(melted_results, sample)

melted_results[, LFC := melted_results_by_sample[
  melted_results, LFC - med_LFC,
  by = .EACHI
]$V1]

################################################################################

median_melted_results <-
  melted_results[, .(
    medLFC = median(LFC),
    FDR = stouffer(FDR)$p
  ),
  by = .(locus_tag, gene, sample)
  ]

#################################################################################

median_melted_results[gene != ".", gene_name_stylized := paste0("italic('", gene, "')")]
median_melted_results[gene == ".", gene_name_stylized := paste0("bold('", locus_tag, "')")]
median_melted_results[gene == "", gene_name_stylized := paste0("bold('", locus_tag, "')")]
median_melted_results[gene == "control", gene_name_stylized := paste0("bold('", locus_tag, "')")]

################################################################################

plot_object <- melted_results %>%
  mutate(sample = case_when(
    sample == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    sample == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    sample == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = LFC)) +
  geom_density(aes(fill = type), alpha = 0.35) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10^(1:6)),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  facet_wrap(facets = c("sample")) +
  ggtitle("Counts per Million") +
  scale_fill_viridis(discrete = TRUE, option = "mako") +
  theme_bw() +
  ggtitle("Knockdown-induced Guide Composition Change (log-fold change)")

print(plot_object)

################################################################################

plot_object <- median_melted_results %>%
  mutate(sample = case_when(
    sample == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    sample == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    sample == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = medLFC)) +
  geom_density(aes(fill = sample), alpha = 0.35) +
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

genes_of_interest <- c(
  "pgsA",
  "orfN",
  "cysS",
  "purA",
  "purB",
  "purE",
  "purH",
  "purK",
  "purL",
  "purN",
  "lptA",
  "lptH",
  "lptB",
  "lptC",
  "lptD",
  "lptE",
  "lptF",
  "lptG",
  "rpoN",
  "ispD",
  "ispG",
  "mrfp",
  "rep"
)

p_load(tidyverse, ggplot2, data.table, ggrepel, hrbrthemes, viridis, ggallin)
doc_theme <- theme_ipsum(
  base_family = "Arial",
  caption_margin = 12,
  axis_title_size = 12,
  axis_col = "black"
)

volcano_plot <- median_melted_results %>%
  filter(gene != "control") %>%
  mutate(sample = case_when(
    sample == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    sample == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    sample == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = medLFC, y = FDR)) +
  geom_point(
    aes(
      alpha = FDR <= 0.05 & abs(medLFC) >= 1,
      color = gene %in% genes_of_interest
    )
  ) +
  facet_wrap(facets = c("sample")) +
  doc_theme +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_label_repel(
    data = . %>%
      group_by(sample) %>%
      filter(gene != "control") %>%
      arrange(FDR) %>%
      mutate(index = seq_len(n())) %>%
      filter(gene %in% genes_of_interest),
    min.segment.length = 0,
    parse = TRUE,
    max.overlaps = Inf,
    aes(
      label = gene_name_stylized,
      color = FDR < 0.05 & abs(medLFC) > 1
    ),
    alpha = 1
  ) +
  scale_alpha_manual(values = c(0.5, 1), guide = "none") +
  scale_color_manual(values = c("dark gray", "red"), guide = "none") +
  scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans())

print(volcano_plot)



