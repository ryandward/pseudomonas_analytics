source("cleanup_analysis.R")
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

count_stats.mat.quality <-
  count_stats.mat %>%
  pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
  inner_join(exp_design_botneck, by = "condition") %>%
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
  mutate(media = factor(media, levels = c("inoculum", "LB", "mouse")))

# Calculate densities and summarize results
density_data <- plot_data %>%
  mutate(log_cpm = log1p(cpm)) %>% # Apply log transformation
  group_by(type, media, condition) %>%
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
  group_by(type, media, x) %>%
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
ribbon_colors <- brewer.pal(n = 12, name = "Paired")[c(1, 3, 5)]
line_colors <- brewer.pal(n = 12, name = "Paired")[c(2, 4, 6)]

# Plot with geom_ribbon for density range and geom_line for the middle line
plot_object <- ggplot() +
  geom_ribbon(data = density_data, aes(x = x, ymin = ymin, ymax = ymax, fill = media)) +
  geom_line(data = density_data, aes(x = x, y = ymean, color = media), lwd = 0.75) +
  geom_line(data = density_data, aes(x = x, y = ymin, color = media), linetype = "dashed", lwd = 0.35) +
  geom_line(data = density_data, aes(x = x, y = ymax, color = media), linetype = "dashed", lwd = 0.35) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 100, 10000, 100000, 1000000),
    labels = custom_log_labels(10)
  ) +
  facet_grid(facets = c("type", "media"), scales = "free_y") +
  # ggtitle("Distribution of Guides Recovered (Counts per Million) with SEM") +
  scale_fill_manual(values = ribbon_colors) +
  scale_color_manual(values = line_colors) +
  theme_minimal()

print(plot_object)
