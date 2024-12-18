source("Enrichments/guide_level_gene_enrichment.r")
source("final_analytical_questions.r")

require(ggbeeswarm)


# this_term <- "KW-0658" # Figure 3
# this_term <- "CL:957" # Figure 3
this_term <- "CL:2388" # Supplementary Figure
# this_term <- "CL:924" # Supplementary Figure
# this_term <- "KW-0289" # Supplementary Figure

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

# Function to calculate weighted median
weighted_median <- function(y, w) {
  order <- order(y)
  y <- y[order]
  w <- w[order]
  cumw <- cumsum(w)
  cutoff <- sum(w) / 2
  median_index <- which(cumw >= cutoff)[1]
  return(y[median_index])
}

# Data preparation
processed_data <- contrast_assignments %>%
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
    FDR > 0.05 ~ "No Response",
    Direction == "Up" ~ "Resilient",
    Direction == "Down" ~ "Vulnerable"
  )) %>%
  mutate(
    contrast = case_when(
      contrast == "plated_6_generations_LB - plated_t0_inoculum" ~ "in-vitro vs. initial",
      contrast == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "in-vivo vs. initial",
      contrast == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "in-vivo vs. in-vitro"
    )
  ) %>%
  filter(contrast != "in-vivo vs. initial") %>%
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
    assignment, `Guide-level\nLog-fold change`, `Guide-level\nFDR`, label, FDR, spacer, gene
  ) %>%
  unique()

# Calculate weighted median for each label
weighted_medians <- processed_data %>%
  filter(FDR <= 0.05) %>%
  group_by(label) %>%
  summarize(
    weighted_median = weighted_median(`Guide-level\nLog-fold change`, -log10(`Guide-level\nFDR`))
  )

# Define the enrichment plot
enrichment_plot <- ggplot(
  processed_data,
  aes(x = label, y = `Guide-level\nLog-fold change`)
) +
  # draw a horizontal line at 0
  geom_hline(yintercept = 0, color = "dark grey", lwd = 1.5, lty = "solid") +
  geom_boxplot(
    color = "black",
    alpha = 0,
    lwd = 0.15,
    outlier.shape = NA
  ) +
  # geom_violin(
  #   fill = "grey",
  #   alpha = 0.25
  # ) +
  # geom_quasirandom(
  #   data = processed_data,
  #   aes(
  #     size = `Guide-level\nFDR`,
  #   ),
  #   shape = 21,
  #   color = "#000000bb",
  #   fill = "#000000bb",
  #   alpha = 0.75
  # ) +
  geom_segment(
    data = weighted_medians,
    aes(x = as.numeric(label) - 0.375, xend = as.numeric(label) + 0.375, y = weighted_median, yend = weighted_median),
    color = "red",
    lwd = 2,
    lty = "solid",
    alpha = 0.5
  ) +
  geom_sina(
    data = processed_data,
    aes(
      size = `Guide-level\nFDR`,
    ),
    fill = "grey",
    shape = 21,
    alpha = 0.5
  ) +
  # geom_label(
  #   data = weighted_medians,
  #   aes(x = as.numeric(label), y = weighted_median, label = round(weighted_median, 2)),
  #   # vjust = -0.5,
  #   color = "red"
  # ) +
  # scale_size(
  #   range = c(0, 4)
  # ) +
  labs(
    x = NULL,
    y = "Guide-level log-fold change"
  ) +
  ggtitle(title) +
  scale_size_continuous(
    range = c(10, 0),
    trans = "log10",
    breaks = c(1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10),
    limits = c(min(processed_data$`Guide-level\nFDR`) / 10, 1e0),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, size = 8),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.title = element_text(hjust = 0.5, size = 8)
  )

# Plot the enrichment plot
plot(enrichment_plot)

###

# Calculate weighted median for each label

control_set <- control_set %>%
  mutate(Direction = case_when(
    PValue <= 0.05 ~ Direction,
    TRUE ~ "No Response"
  ))

control_labels <- control_set %>%
  rename(group = contrast) %>%
  mutate(
    group = case_when(
      group == "plated_6_generations_LB - plated_t0_inoculum" ~ "in-vitro vs. initial",
      group == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "in-vivo vs. initial",
      group == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "in-vivo vs. in-vitro"
    )
  ) %>%
  mutate(
    label = paste(
      group,
      paste(
        "p-value:",
        signif(PValue, 2)
      ),
      Direction,
      sep = "\n"
    )
  ) %>%
  mutate(
    group = factor(group, levels = unique(group))
  ) %>%
  mutate(
    label = factor(label, levels = unique(label))
  )

control_weighted_medians <- non_normalized_melted_results %>%
  filter(type == "control") %>%
  group_by(condition) %>%
  summarize(
    weighted_median = weighted_median(LFC, -log10(FDR))
  ) %>%
  mutate(
    group = case_when(
      condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "in-vitro vs. initial",
      condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "in-vivo vs. initial",
      condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "in-vivo vs. in-vitro"
    )
  ) %>%
  mutate(group = factor(group, levels = unique(group))) %>%
  inner_join(control_labels)

non_normalized_melted_results %>%
  rename(group = condition) %>%
  inner_join(
    control_set %>%
      rename(group = contrast)
  ) %>%
  mutate(
    group = case_when(
      group == "plated_6_generations_LB - plated_t0_inoculum" ~ "in-vitro vs. initial",
      group == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "n-vivo vs. initial",
      group == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "in-vivo vs. in-vitro"
    )
  ) %>%
  mutate(
    group = factor(group, levels = unique(group))
  ) %>%
  inner_join(control_labels) %>%
  filter(type == "control") %>%
  ggplot(aes(x = label, y = LFC)) +
  geom_hline(yintercept = 0, color = "dark grey", lwd = 2, lty = "solid") +
  geom_boxplot(outlier.shape = NA, alpha = 0.25) +
  geom_quasirandom(aes(size = FDR),
    color = "NA",
    fill = "black",
    alpha = 0.3,
    shape = 21,
    alpha = 0.5
  ) +
  scale_size_continuous(range = c(10, 0), trans = "log10") +
  geom_segment(
    data = control_weighted_medians,
    aes(x = as.numeric(group) - 0.4, xend = as.numeric(group) + 0.4, y = weighted_median, yend = weighted_median),
    color = "red",
    lwd = 3,
    lty = "solid",
    alpha = 0.5
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = "Guide-level log-fold change"
  )
