library(pacman)

p_load(data.table, tidyverse, edgeR, poolr, scales, pheatmap, viridis, ggallin, ggforce)

count_stats <- fread("stats.tsv", col.names = c("spacer", "count", "condition"))

designed_targets <- fread("designed_targets.tsv", na.strings = c("None")) %>%
mutate(type = case_when(
  mismatches == 0 ~ "knockdown",
  is.na(target) ~ "control"
)) %>% filter(is.na(sp_dir) | sp_dir != tar_dir)

# count_stats <- count_stats %>% filter(spacer %in% designed_targets$spacer)

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
  inner_join(designed_targets) %>%
  pivot_wider(id_cols = c(type, spacer), names_from = condition, values_from = count, values_fill = 0) %>%
  pivot_longer(cols = -c(type, spacer), names_to = "condition", values_to = "count") %>%
  group_by(condition) %>%
  mutate(cpm = 1e6 * (count / sum(count))) %>%
  group_by(condition, type) %>%
  summarise(cpm = sum(cpm)) %>%
  pivot_wider(id_cols = condition, names_from = type, values_from = cpm) %>%
  inner_join(exp_design) %>%
  select(condition, control, knockdown, media, gDNA_source, match_prop, reads)

exp_design <- exp_design %>%
  inner_join(quality_decision) %>%
  filter(match_prop > 0.85 & reads > 5e6) %>%
  group_by(sample_group) %>%
  mutate(rep = as.character(seq_len(n())))

exp_design <- exp_design %>% mutate(sample_name = paste(gDNA_source, growth_condition, media, rep, sep = "_"))

exp_design <- exp_design %>% mutate(
  in_vitro = ifelse(media %in% c("LB", "inoculum"), "in_vitro", "in_vivo"),
)

exp_design <- exp_design %>%  mutate(
  outgrowth = ifelse(growth_condition %in% c("6_generations", "10x_inoculum_dilution"), "outgrowth", "inoculum")
)


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
  inner_join(designed_targets) %>%
  group_by(spacer, condition) %>%
  inner_join(exp_design) %>%
  pivot_wider(id_cols = spacer, names_from = sample_name, values_fill = 0, values_from = count) %>%
  ungroup() %>%
  select(-spacer) %>%
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
  inner_join(designed_targets) %>%
  pivot_wider(id_cols = c(type, spacer), names_from = condition, values_from = count, values_fill = 0)

count_stats.mat %>%
  pivot_longer(cols = -c(type, spacer), names_to = "condition", values_to = "count") %>%
  inner_join(designed_targets)

count_stats <- count_stats.mat %>%
  pivot_longer(cols = -c(type, spacer), names_to = "condition", values_to = "count") %>%
  inner_join(designed_targets) %>%
  unite("guide", c(gene, offset), remove = FALSE)

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
  pivot_longer(cols = -c(type, spacer), names_to = "condition", values_to = "cpm") %>%
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
  pivot_longer(cols = -c(type, spacer), names_to = "condition", values_to = "count") %>%
  inner_join(exp_design) %>%
  # filter(gDNA_source == "plated") %>%
  filter(type != "focused") %>%
  pivot_wider(id_cols = c(type, spacer), names_from = condition, values_from = count, values_fill = 0)

plot_object <- count_stats.mat.quality %>%
  select(-type, -spacer) %>%
  data.matrix() %>%
  cpm() %>%
  as_tibble() %>%
  cbind(
    count_stats.mat.quality %>%
      select(type, spacer)
  ) %>%
  pivot_longer(cols = -c(type, spacer), names_to = "condition", values_to = "cpm") %>%
  inner_join(exp_design) %>%
  # unite("growth_condition", c(growth_condition, media)) %>%
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
  pivot_longer(cols = -c(type, spacer), names_to = "condition", values_to = "count") %>%
  inner_join(designed_targets) %>%
  unite("guide", c(gene, offset), remove = FALSE) %>%
  mutate(guide = case_when(type == "control" ~ paste(spacer), type == "knockdown" ~ paste(gene, offset)))

setorder(count_stats.quality, condition)
setorder(exp_design, condition)

conditions.quality <-
  count_stats.quality %>%
  select(condition) %>%
  unique()

##########################################################################################

data_grid <- count_stats.quality %>%
  pivot_wider(id_cols = spacer, names_from = condition, values_from = count)

data_grid_matrix <- data_grid %>%
  select(-spacer) %>%
  data.matrix()

row.names(data_grid_matrix) <- data_grid$spacer

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

rownames(data_permut) <- colnames(data_grid_matrix)

# get rid of "timing" and "drug_dose" in names of columns of aba_permut
colnames(data_permut) <- gsub("data_group", "", colnames(data_permut))

##########################################################################################

data_y <- DGEList(
  counts = data_grid_matrix,
  group = data_group,
  genes = row.names(data_grid_matrix)
)

# data_keep <- filterByExpr(data_grid_matrix, data_group)

# data_y <- data_y[data_keep, , keep.lib.sizes = FALSE]

data_y <- calcNormFactors(data_y)

data_y <- estimateGLMRobustDisp(data_y, data_permut)

data_fit <- glmQLFit(data_y, data_permut, robust = TRUE)

data_CPM <- cpm(data_y, prior.count = 0)

################################################################################

data_contrast <- makeContrasts(
  in_vivo = plated_10x_inoculum_dilution_mouse - plated_t0_inoculum,
  in_vitro = plated_6_generations_LB - plated_t0_inoculum,
  in_vivo_vs_in_vitro = plated_10x_inoculum_dilution_mouse - plated_6_generations_LB,
  levels = data_permut
)

contrast_levels <- colnames(data_contrast)

################################################################################

results_FDR <- data.table(count_stats.quality)[, .(genes = unique(spacer))]

results_LFC <- data.table(count_stats.quality)[, .(genes = unique(spacer))]

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
  rename(spacer = genes) %>%
  inner_join(designed_targets)

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



all_string <- fread("STRG0A01FJP.protein.enrichment.terms.v12.0.txt.gz") %>%
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) %>%
  unique()



gene_groups <- all_string %>%
      # filter(term %in% (all_string %>% group_by(term) %>% tally() %>% pull(unique(term)))) %>%
      group_by(category, term, description) %>%
      summarise(gene_count = n(), locus_tag = list(sort(unique(locus_tag)))) %>%
      mutate(locus_tag_group = vapply(locus_tag, paste, collapse = ",", FUN.VALUE = character(1)))

targets <- fread("designed_targets.tsv")

targets <- targets %>% filter(spacer %in% count_stats$spacer)

term_stats <- gene_groups %>% 
  unnest(locus_tag) %>% 
  inner_join(targets %>% 
  select(locus_tag) %>% unique()) %>% 
  group_by(term, gene_count, description) %>% 
  summarize(genes_targeted = n()) %>%
  select(term, description, gene_count, genes_targeted)

complete_terms <- term_stats %>% 
  filter(gene_count == genes_targeted)

# only perform enrichments where all genes are available
# gene_groups <- complete_terms %>% inner_join(gene_groups)

repeated_gene_groups <- gene_groups %>% 
  group_by(locus_tag) %>% 
  mutate(times_listed = n()) %>% 
  arrange(locus_tag) %>% 
  ungroup ()


# pick the best annotation for each locus_tag_group, i.e., highest in term, and the lowest in the category_rank
ranked_annotations <- repeated_gene_groups %>%
      group_by(locus_tag_group, category) %>%
      arrange(versionsort::ver_sort(term)) %>%
      slice(n()) %>%
      ungroup() %>%
      mutate(category_rank = case_when(
        category == "Biological Process (Gene Ontology)" ~ 1,
        category == "Molecular Function (Gene Ontology)" ~ 2,
        category == "Cellular Component (Gene Ontology)" ~ 3,
        category == "Protein Domains and Features (InterPro)" ~ 4,
        category == "Protein Domains (SMART)" ~ 5,
        category == "Protein Domains (Pfam)" ~ 6,
        category == "Annotated Keywords (UniProt)" ~ 7,
        category == "Reactome Pathways" ~ 8,
        category == "Subcellular localization (COMPARTMENTS)" ~ 9,
        category == "Local Network Cluster (STRING)" ~ 10,
        TRUE ~ NA_integer_
      )) %>%
      group_by(locus_tag_group) %>%
      filter(category_rank == min(category_rank))

enrichments <- ranked_annotations %>%
  ungroup() %>%
  distinct(locus_tag_group, .keep_all = TRUE) %>%
  select(-locus_tag_group) %>%
  unnest(locus_tag) %>%
  inner_join(term_stats) 


# Get the unique terms
unique_terms <- unique(enrichments$term)

target_spacers_for_terms <- term_stats %>%
  inner_join(enrichments, relationship = "many-to-many") %>%
  inner_join(targets, relationship = "many-to-many")

####################


#########################################################################################

# Split the spacer column by term
locus_tags_list <- split(target_spacers_for_terms$spacer, target_spacers_for_terms$term)

dge <- data_y

# Find the indices of each set of locus tags in rownames(dge)
gene_indices <- lapply(locus_tags_list, function(locus_tags) which(rownames(dge) %in% locus_tags))


v <- voomWithQualityWeights(dge, data_permut, plot = TRUE)


v_targets <- v$E %>%
  data.table(keep.rownames = "spacer") %>%
  select(spacer) %>%
  left_join(
    targets %>%
      filter(locus_tag %in% all_string$locus_tag) %>%
      group_by(spacer) %>%
      filter(is.na(target) | target == "None" | (sp_dir != tar_dir & abs(as.numeric(offset)) == min(abs(as.numeric(offset))) & overlap == max(overlap)))
      %>%
      group_by(target)
  )

# Perform the competitive gene set test for all gene sets
all_sets <- lapply(colnames(data_contrast), function(contrast_name) {
  contrast_column <- data_contrast[, contrast_name]
  result <- camera(
    v,
    index = gene_indices, design = data_permut,
    # weights = v_targets$weight,
    inter.gene.cor = 0.05,
    contrast = contrast_column
  ) %>%
    data.table(keep.rownames = "term") %>%
    mutate(term = factor(term, levels = unique_terms), contrast = contrast_name)
  result
}) %>%
  do.call(rbind, .)

  #######################


best_sets <- all_sets %>% 
  inner_join(enrichments) %>% 
  inner_join(v_targets) %>%
  inner_join(term_stats) %>% 
  group_by(contrast, term, description) %>% 
  nest(locus_tags = locus_tag) %>% 
  group_by(locus_tags, contrast) %>% 
  mutate(missing_genes = gene_count - genes_targeted) %>%
  arrange(FDR, missing_genes) %>% 
  slice(1) %>% 
  ungroup %>% 
  rename(guide_count = NGenes) %>%
  select(term, guide_count, Direction, PValue, FDR, contrast, description, genes_targeted, gene_count)  %>% unique()  %>% 
  data.table()


################################################################################

best_sets %>% 
  filter(contrast == "in_vitro") %>%
  filter(FDR<=0.05) %>% arrange(FDR) %>% head(12) %>% 
  inner_join(enrichments) %>% 
  inner_join(count_stats %>% group_by(condition) %>% mutate(cpm = cpm(count))) %>% 
  ungroup %>%
  inner_join(exp_design) %>%
  mutate(facet_title = paste(term, paste0("(", guide_count, " guides ", toupper(Direction), ": ", signif(FDR,3), ")"))) %>%
  mutate(facet_title = factor(facet_title, levels = facet_title %>% unique())) %>%
  ggplot(aes(y = (cpm), x = factor(sample_group))) +
  geom_sina(aes(color = factor(sample_group)), scale = "count") +
  geom_violin(aes(fill = factor(sample_group), color = factor(sample_group)), alpha = 0.25, scale = "count", draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~ facet_title + paste(contrast, description, sep = "\n"), ncol = 6, scales = "free_y") +
  theme_minimal() +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(10^(0:5)),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) 

























