library(pacman)
p_load(data.table, tidyverse, edgeR, poolr, scales, pheatmap, viridis)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

count_stats <- fread("stats.tsv", col.names = c("sequence", "count", "condition"))

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

map_stats <- fread("map_stats.tsv")

exp_design <- fread("unified_counts/exp_design.tsv")

exp_design <- exp_design %>% inner_join(map_stats)

exp_design <- exp_design %>%
  # this one didn't get infected
  # filter(!growth_condition %like% "100x_inoculum_dilution") %>%
  filter(!condition %in% c("dJMP3", "dJMP5")) %>%
  # filter(!condition %in% c("dJMP1")) %>%

  data.table()

exp_design_density_plots <- exp_design %>% mutate(sample_group = paste(gDNA_source, growth_condition, media, sep = "_"))

exp_design <- exp_design_density_plots %>%
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



count_stats.mat.quality <-
  count_stats.mat %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  inner_join(exp_design_density_plots, by = "condition") %>%
  filter(type != "focused") %>%
  pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0)

plot_data <- count_stats.mat.quality %>%
  select(-type, -sequence) %>%
  data.matrix() %>%
  cpm() %>%
  as_tibble() %>%
  cbind(
    count_stats.mat.quality %>%
      select(type, sequence)
  ) %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "cpm") %>%
  inner_join(exp_design, by = "condition")

# Reorder media levels
plot_data <- plot_data %>%
  mutate(growth_condition = factor(growth_condition, levels = c("t0", "6_generations", "10x_inoculum_dilution", "100x_inoculum_dilution")))

# Calculate densities and summarize results
density_data <- plot_data %>%
  mutate(log_cpm = log1p(cpm)) %>% # Apply log transformation
  group_by(type, growth_condition, condition) %>%
  summarise(
    density = list(density(log_cpm, from = 0, to = log1p(10000))),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(density, ~ expm1(.x$x)), # Transform x-values back to original scale
    y = map(density, ~ .x$y)
  ) %>%
  select(-density) %>%
  unnest(c(x, y)) %>%
  group_by(type, growth_condition, x) %>%
  summarise(
    ymax = max(y),
    ymin = min(y),
    ymean = mean(y), # Calculate mean for the middle line
    .groups = "drop"
  )

# Custom labeling function to include 0
custom_log_labels <- function(base = 10) {
  function(x) {
    labels <- scales::label_log(base)(x)
    labels[x == 0] <- "0"
    labels
  }
}

# Define custom colors for media using specific colors from the Paired palette
ribbon_colors <- brewer.pal(n = 12, name = "Paired")[c(1, 3, 5, 7)]
line_colors <- brewer.pal(n = 12, name = "Paired")[c(2, 4, 6, 8)]

# Plot with geom_ribbon for density range and geom_line for the middle line
plot_object <- ggplot() +
  geom_ribbon(data = density_data, aes(x = x, ymin = ymin, ymax = ymax, fill = growth_condition)) +
  geom_line(data = density_data, aes(x = x, y = ymean, color = growth_condition), lwd = 0.75) +
  geom_line(data = density_data, aes(x = x, y = ymin, color = growth_condition), linetype = "dashed", lwd = 0.35) +
  geom_line(data = density_data, aes(x = x, y = ymax, color = growth_condition), linetype = "dashed", lwd = 0.35) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 100, 10000, 100000, 1000000),
    labels = custom_log_labels(10)
  ) +
  facet_grid(facets = c("type", "growth_condition"), scales = "free_y") +
  # ggtitle("Distribution of Guides Recovered (Counts per Million) with SEM") +
  scale_fill_manual(values = ribbon_colors) +
  scale_color_manual(values = line_colors) +
  # dont show legend for fill or color
  theme_minimal() +
  theme(legend.position = "none")

print(plot_object)

# immediately source cleanup_analysis.R to prevent code duplication
source("cleanup_analysis.R")
