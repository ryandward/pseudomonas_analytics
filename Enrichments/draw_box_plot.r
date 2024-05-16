source("Enrichments/guide_level_gene_enrichment.r")
source("final_analytical_questions.r")

# this_term <- "GO:0045229" # with orfN 
# this_term <- "CL:2706" # with orfN
# this_term <- "GO:0046486" #with pgs
this_term <- "KW-0658" # purine biogenesis
# this_term <- "CL:2388" # large
# this_term <- "GO:0046474"
# this_term <- "GO:0016780" # with orfN and pgsA
# this_term <- "GO:0006629" #lipid metabolic process with orfN
title <- term_stats %>%
  filter(term == this_term) %>%
  pull(description)


title <- paste(title, " (", this_term, ")", sep = "")

title <- title %>% stringr::str_wrap(width = 50)

this_genes_targeted <- all_sets %>%
  filter(term == this_term) %>%
  pull(genes_targeted)

this_gene_count <- all_sets %>%
  filter(term == this_term) %>%
  pull(gene_count)

title <- paste(
  title, "\n",
  paste(
    this_genes_targeted, this_gene_count,
    sep = "/"
  ),
  " genes present",
  sep = " "
)

# non_normalized_melted_results <- non_normalized_melted_results_LFC %>% 
#   inner_join(non_normalized_melted_results_FDR)

enrichment_plot <- contrast_assignments %>%
  filter(assignment == 1) %>%
  inner_join(
    non_normalized_melted_results %>%
      rename(contrast = condition, `Guide-level\nFDR` = FDR)
  ) %>%
  inner_join(
    group_assignments,
    relationship = "many-to-many"
  ) %>%
  mutate(`Guide-level\nLog-fold change` = case_when(
    assignment == -1 ~ 0,
    assignment == 1 ~ LFC
  )) %>%
  inner_join(
    all_sets %>%
      filter(term == this_term) %>%
      inner_join(contrast_assignments) %>%
      mutate(contrast = factor(contrast, levels = unique(contrast)))
  ) %>%
  mutate(Direction = case_when(
    FDR <= 0.05 ~ Direction,
    TRUE ~ "No change"
  )) %>%
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
      Direction,
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
  arrange(abs(`Guide-level\nLog-fold change`)) %>%
  mutate(
    group = factor(group, levels = unique(group))
  ) %>%
  select(
    assignment, `Guide-level\nLog-fold change`, `Guide-level\nFDR`, label, FDR, spacer
  ) %>%
  unique() %>%
  ggplot(
    aes(x = label, y = `Guide-level\nLog-fold change`)
  ) +
  # draw a horizontal line at 0
  geom_hline(yintercept = 0, color = "darkgrey", lwd = 1.5, lty = "dashed") +
  # geom_tile(
  #   aes(alpha = factor(ifelse(FDR <= 0.05, "highlight", "no_highlight"))),
  #   width = Inf, height = Inf, fill = "light grey"
  # ) +
  geom_sina(
    aes(
      fill = `Guide-level\nLog-fold change`,
      size = `Guide-level\nFDR`,
      # weight = -log10(`Guide-level\nFDR`)
    ),
    shape = 21,
    color = "black",
    # lwd = 0.1,
    scale = "count"
  ) +
  geom_boxplot(
    aes(
      # weight = -log10(`Guide-level\nFDR`)
    ),
    alpha = 0.0,
    # draw_quantiles = c(0.25, 0.5, 0.75),
    # scale = "width",
    lwd = 0.5
  ) +
  scale_alpha_manual(
    values = c("highlight" = 0.00, "no_highlight" = 0.025), guide = FALSE
  ) +
  scale_size(
    range = c(0.1, 3)
  ) +
  # scale_y_continuous(
  #   trans = scales::pseudo_log_trans(base = 10),
  #   breaks = c(10^(0:5)),
  #   labels = scales::label_number(scale_cut = scales::cut_short_scale())
  # ) +
  # facet_wrap(~label) +
  # label x axis "Assignment" and y axis "Counts per million"
  labs(
    x = NULL,
    y = "Guide-level log-fold change"
  ) +
  ggtitle(title) +
scale_fill_gradient2(
  low = rgb(1, 0, 0, alpha=0.5), 
  mid = rgb(1, 1, 1, alpha=0.5), 
  high = rgb(0, 0, 1, alpha=0.5), 
  midpoint = 0
) +
  scale_size_continuous(range = c(10, 2), trans = "log10") +
  theme_minimal() +
  # turn angle to 45 degrees
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5)
  )

plot(enrichment_plot)