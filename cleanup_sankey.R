source("cleanup_analysis.R")

p_load(ggsankey, tidyverse, viridis)

sankey.levels <- c(
  "Resistant",
  "No Response",
  "Vulnerable"
)

sankey.df <- median_melted_results %>%
  mutate(
    Response = case_when(
      medLFC < -1 & FDR < 0.05 ~ "Vulnerable",
      medLFC > 1 & FDR < 0.05 ~ "Resistant",
      TRUE ~ "No Response"
    )
  ) %>%
  select(locus_tag, condition, Response) %>%
  mutate(condition = case_when(
    condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  pivot_longer(!c(condition, locus_tag)) %>%
  pivot_wider(id_cols = c(locus_tag), names_from = condition, values_from = value) %>%
  make_long(`In vitro`, `In vivo v. in vitro`) %>%
  mutate(
    node = factor(node, levels = sankey.levels),
    next_node = factor(next_node, levels = sankey.levels)
  )

sankey.tally <- sankey.df %>%
  group_by(node, x) %>%
  tally()

sankey.df %>%
  inner_join(sankey.tally) %>%
  ggplot(aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = paste(node, n, sep = "\n")
  )) +
  geom_sankey(
    flow.colour = "black",
    flow.alpha = 0.25,
    node.color = "black"
  ) +
  scale_fill_manual(values = c(
    "light blue", "grey", "red"
  )) +
  geom_sankey_label(
    aes(colour = "node"),
    size = 3.5,
    color = 1
  ) +
  theme_sankey(base_size = 16) +
  guides(
    fill = guide_legend(title = "Relative Response to Knockdown")
  ) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
