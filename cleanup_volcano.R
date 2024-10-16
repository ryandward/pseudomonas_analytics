source("cleanup_analysis.R")
source("final_analytical_questions.r")

p_load(tidyverse, ggplot2, data.table, ggrepel, hrbrthemes, viridis, ggallin)
doc_theme <- theme_ipsum(
  base_family = "Arial",
  caption_margin = 12,
  axis_title_size = 12,
  axis_col = "black"
)
# top_ten <- median_melted_results %>% arrange(FDR) %>% select(gene) %>% filter(gene != "control") %>% unique %>% head(10)
median_melted_results %>%
  filter(gene != "control") %>%
  mutate(condition = case_when(
    condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = medLFC, y = FDR)) +
  geom_point(aes(colour = FDR < 0.05 & abs(medLFC) > 1)) +
  facet_wrap(facets = c("condition")) +
  doc_theme +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("gray", "red"), na.value = "grey") +
  geom_label_repel(
    data = . %>%
      group_by(condition) %>%
      filter(gene != "control") %>%
      arrange(FDR) %>%
      mutate(index = seq_len(n())) %>%
      filter(gene %in% c("pgsA", "orfN", "cysS")),
    min.segment.length = 0,
    parse = TRUE,
    max.overlaps = Inf,
    aes(label = gene)
  ) +
  scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans())

################################################################################

non_normalized_melted_results %>%
  filter(gene != "control") %>%
  mutate(condition = case_when(
    condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = LFC, y = FDR)) +
  geom_point(aes(colour = FDR < 0.05 & abs(LFC) > 1)) +
  facet_wrap(facets = c("condition")) +
  doc_theme +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("gray", "red"), na.value = "grey") +
  geom_label_repel(
    data = . %>%
      group_by(condition) %>%
      filter(gene != "control") %>%
      arrange(FDR) %>%
      mutate(index = seq_len(n())) %>%
      filter(index <= 10 & FDR < 0.05),
    # filter(gene %in% top_ten$gene & FDR < 0.05),
    min.segment.length = 0,
    max.overlaps = Inf,
    aes(label = guide)
  ) +
  scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans())

################################################################################

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
  "mrfp"
)

median_melted_results <- median_melted_results %>%
  group_by(condition) %>%
  arrange(FDR) %>%
  filter(locus_tag != "control") %>%
  filter(FDR <= 0.05) %>%
  slice(1:25) %>%
  select(locus_tag, condition) %>%
  mutate(top_five = TRUE) %>%
  right_join(median_melted_results) %>%
  mutate(top_five = case_when(is.na(top_five) ~ FALSE, TRUE ~ TRUE))

median_melted_results <- median_melted_results %>%
  mutate(top_five = case_when((FDR <= 0.05 & gene %in% c("pgsA")) | top_five == TRUE ~ TRUE, TRUE ~ FALSE))
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggpointdensity) # Required for coloring points by density

# Define a common theme
common_theme <- theme_minimal() +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Calculate x-axis limits
x_limits <- range(median_melted_results$medLFC, na.rm = TRUE)

# Define modulus transformation and its inverse
p_value <- 0.25
trans <- modulus_trans(p = p_value)

# Generate evenly spaced points in the original range from 0 to 300
n_breaks <- 5
original_space <- seq(0, 300, length.out = 1000)

# Apply modulus transformation to the original points
transformed_values <- trans$transform(original_space)

# Select evenly spaced transformed points (between min and max)
selected_transformed_breaks <- seq(min(transformed_values), max(transformed_values), length.out = n_breaks)

# Apply the inverse transformation to get original values for breaks
selected_original_breaks <- trans$inverse(selected_transformed_breaks)

# Ensure that 1 is included in the breaks
selected_original_breaks <- unique(round(c(1, selected_original_breaks)))

# Create the main plot with density-based point coloring
volcano_plot_density <- median_melted_results %>%
  filter(gene != "control") %>%
  filter(condition != "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum") %>%
  mutate(condition = case_when(
    condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = medLFC, y = FDR)) +
  # Add points colored by density, adjust point shape and size
  geom_pointdensity(aes(alpha = FDR <= 0.05 & abs(medLFC) >= 1), shape = 20, size = 3.5) + # Filled circle, larger size
  scale_color_viridis_c(
    option = "plasma", # Use a different color palette for better contrast
    end = 0.9,
    name = "Local Density",
    trans = trans,
    breaks = selected_original_breaks # Use calculated breaks that include 1 and are rounded for clarity
  ) +
  scale_alpha_manual(
    values = c("TRUE" = 0.5, "FALSE" = 0.15), # Increase opacity for non-significant points
    guide = "none" # Hide the significance legend
  ) +
  facet_wrap(facets = c("condition")) +
  common_theme +
  # Updated threshold lines to be black, thicker, and dashed
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", lwd = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", lwd = 1) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", lwd = 1) +
  geom_label_repel(
    data = . %>%
      group_by(condition) %>%
      filter(gene != "control") %>%
      arrange(FDR) %>%
      mutate(index = seq_len(n())) %>%
      filter(top_five == TRUE),
    min.segment.length = 0,
    parse = TRUE,
    max.overlaps = Inf,
    aes(
      label = gene_name_stylized
    ),
    alpha = 0.8
  ) +
  scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans()) +
  coord_cartesian(xlim = x_limits)

# Create the marginal density plot
marginal_plot <- median_melted_results %>%
  filter(gene != "control") %>%
  mutate(condition = case_when(
    condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = medLFC)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  facet_wrap(facets = c("condition")) +
  common_theme +
  theme(legend.position = "none") +
  coord_cartesian(xlim = x_limits)

# Combine the main plot and the marginal plot using cowplot
combined_plot <- plot_grid(
  volcano_plot_density,
  marginal_plot,
  ncol = 1,
  align = "v",
  rel_heights = c(3, 1) # Adjust the relative heights
)

# print(combined_plot)

print(volcano_plot_density + theme(panel.spacing = unit(5, "lines")))
