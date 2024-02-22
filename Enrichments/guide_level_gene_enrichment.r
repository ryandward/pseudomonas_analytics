source("cleanup_analysis.R")

p_load("ggforce")

dge <- data_y

colnames(dge$design) <- gsub("[^[:alnum:]]", "_", colnames(dge$design))
colnames(dge$design) <- gsub("___", " - ", colnames(dge$design))

dge$samples$group <- gsub("[^[:alnum:]]", "_", dge$samples$group)
dge$samples$group <- gsub("___", " - ", dge$samples$group)

contrast_levels <- gsub("[^[:alnum:]]", "_", contrast_levels)
contrast_levels <- gsub("___", " - ", contrast_levels)

contrast_list <- setNames(as.list(contrast_levels), contrast_levels)

contrast_list$levels <- dge$design

contrasts <- do.call(makeContrasts, contrast_list)

all_string <- fread("Enrichments/STRG0A01FJP.protein.enrichment.terms.v12.0.txt.gz") %>%
  mutate(locus_tag = str_replace(`#string_protein_id`, ".*\\.", "")) %>%
  unique()

targets <- all_guides %>% rename(spacer = sequence)

gene_groups <- all_string %>%
  # filter(term %in% (all_string %>% group_by(term) %>% tally() %>% pull(unique(term)))) %>%
  group_by(category, term, description) %>%
  summarise(gene_count = n(), locus_tag = list(sort(unique(locus_tag)))) %>%
  mutate(locus_tag_group = vapply(locus_tag, paste, collapse = ",", FUN.VALUE = character(1)))

term_stats <- gene_groups %>%
  unnest(locus_tag) %>%
  inner_join(targets %>%
    select(locus_tag) %>% unique()) %>%
  group_by(term, gene_count, description) %>%
  summarize(genes_targeted = n())

repeated_gene_groups <- gene_groups %>%
  group_by(locus_tag) %>%
  mutate(times_listed = n()) %>%
  arrange(locus_tag) %>%
  ungroup()


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

#########################################################################################


# Split the spacer column by term
sets_to_locus_tags <- split(target_spacers_for_terms$spacer, target_spacers_for_terms$term)

# Find the indices of each set of locus tags in rownames(dge)
sets_to_locus_tags_indices <- lapply(sets_to_locus_tags, function(locus_tags) which(rownames(dge) %in% locus_tags))

v <- voomWithQualityWeights(dge, dge$design, plot = TRUE)



v_targets <- v$E %>%
  data.table(keep.rownames = "spacer") %>%
  select(spacer) %>%
  left_join(
    targets
  )


# Perform the competitive gene set test for all gene sets
all_sets <- lapply(colnames(contrasts), function(contrast_name) {
  contrast_column <- contrasts[, contrast_name]
  result <- camera(
    v,
    index = sets_to_locus_tags_indices,
    design = dge$design,
    # inter.gene.cor = 0.05,
    contrast = contrast_column
  ) %>%
    data.table(keep.rownames = "term") %>%
    mutate(term = factor(term, levels = unique_terms), contrast = contrast_name)
  result
}) %>%
  do.call(rbind, .)

all_sets <- all_sets %>%
  inner_join(enrichments) %>%
  inner_join(v_targets) %>%
  inner_join(term_stats) %>%
  group_by(contrast, term, description) %>%
  nest(locus_tags = locus_tag) %>%
  group_by(locus_tags, contrast) %>%
  mutate(missing_genes = gene_count - genes_targeted) %>%
  arrange(FDR, missing_genes) %>%
  ungroup() %>%
  rename(guide_count = NGenes) %>%
  select(term, guide_count, Direction, PValue, FDR, contrast, description, genes_targeted, gene_count) %>%
  unique() %>%
  data.table()

contrast_assignments <- contrasts %>%
  data.table(keep.rownames = "group") %>%
  melt(
    id.vars = "group",
    variable.name = "contrast",
    value.name = "assignment"
  ) %>%
  filter(assignment != 0)

group_assignments <- dge$design %>%
  data.table(keep.rownames = "sample") %>%
  melt(id.vars = "sample", variable.name = "group") %>%
  filter(value != 0) %>%
  select(-value)

original_data <- dge$counts %>%
  data.table(keep.rownames = "spacer") %>%
  melt(
    value.name = "count",
    id.vars = "spacer",
    variable.name = "sample"
  )

annotated_data <- dge$samples %>%
  data.table(keep.rownames = "sample") %>%
  inner_join(original_data) %>%
  group_by(sample) %>%
  mutate(cpm = 1e6 * count / sum(count))

# Export the results #################################################################


all_sets_for_export <- all_sets %>%
  filter(FDR <= 0.05) %>%
  select(term, description, gene_count, genes_targeted, guide_count, FDR, contrast) %>%
  mutate(contrast = case_when(
    contrast == "plated_6_generations_LB - plated_t0_inoculum" ~ "LB vs. INOCULUM",
    contrast == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "mouse vs. INOCULUM",
    contrast == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "mouse vs. LB",
  ))

fwrite(all_sets_for_export, "Enrichments/annotated_gene_set_enrichments.tsv", sep = "\t")


all_sets_for_export_with_genes <- enrichments %>%
  inner_join(targets) %>%
  arrange(gene) %>%
  group_by(term) %>%
  select(term, description, gene) %>%
  unique() %>%
  mutate(genes = paste(gene, collapse = ", ")) %>%
  select(term, description, genes) %>%
  unique() %>%
  inner_join(all_sets_for_export)

fwrite(all_sets_for_export_with_genes, "Enrichments/annotated_gene_set_enrichments_with_genes.tsv", sep = "\t")

### Plot Crafting Area #################################################################
### you could create a function out of this

# this_term <- "CL:2702"
# this_term <- "CL:2730"
# this_term <- "CL:2388"
this_term <- "GO:0006720"

title <- term_stats %>%
  filter(term == this_term) %>%
  pull(description)

title <- paste(title, " (", this_term, ")", sep = "")

enrichment_plot <- contrast_assignments %>%
  inner_join(
    group_assignments,
    relationship = "many-to-many"
  ) %>%
  inner_join(
    all_sets %>%
      filter(term == this_term) %>%
      inner_join(contrast_assignments) %>%
      mutate(contrast = factor(contrast, levels = unique(contrast)))
  ) %>%
  mutate(
    contrast = case_when(
      contrast == "plated_6_generations_LB - plated_t0_inoculum" ~ "LB vs. INOCULUM",
      contrast == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "mouse vs. INOCULUM",
      contrast == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "mouse vs. LB",
    )
  ) %>%
  mutate(
    label = paste(
      contrast,
      paste(
        "FDR:",
        signif(FDR, 2)
      ),
      paste(
        Direction,
        paste(
          paste(
            genes_targeted, gene_count,
            sep = "/"
          ), " genes present",
          sep = ""
        )
      ),
      sep = "\n"
    )
  ) %>%
  mutate(
    label = factor(label, levels = unique(label))
  ) %>%
  inner_join(
    enrichments
  ) %>%
  inner_join(
    targets
  ) %>%
  inner_join(
    annotated_data
  ) %>%
  inner_join(
    v_targets
  ) %>%
  arrange(
    assignment
  ) %>%
  mutate(
    assignment = case_when(
      assignment == 1 ~ "Treatment",
      assignment == -1 ~ "Control"
    )
  ) %>%
  mutate(
    group = case_when(
      group == "plated_6_generations_LB" ~ "LB",
      group == "plated_t0_inoculum" ~ "INOCULUM",
      group == "plated_10x_inoculum_dilution_mouse" ~ "mouse"
    )
  ) %>%
  mutate(
    group = factor(group, levels = unique(group))
  ) %>%
  ggplot(
    aes(x = as.character(assignment), y = cpm)
  ) +
  geom_sina(
    aes(color = group),
    # scale = "width"
  ) +
  geom_violin(
    alpha = 0.35,
    draw_quantiles = c(0.25, 0.5, 0.75),
    # scale = "width",
    lwd = 1.25
  ) +
  geom_tile(
    aes(alpha = factor(ifelse(FDR <= 0.05, "highlight", "no_highlight"))),
    width = Inf, height = Inf, fill = "light grey"
  ) +
  scale_alpha_manual(
    values = c("highlight" = 0.00, "no_highlight" = 0.025), guide = FALSE
  ) +
  scale_size(
    range = c(0.1, 3)
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(10^(0:5)),
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  facet_wrap(~label) +
  # label x axis "Assignment" and y axis "Counts per million"
  labs(
    x = NULL,
    y = "Counts per Million"
  ) +
  ggtitle(title) +
  theme_minimal()


plot(enrichment_plot)
