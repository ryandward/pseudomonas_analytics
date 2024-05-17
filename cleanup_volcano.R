source("cleanup_analysis.R")

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
      slice(1:25) %>% select(locus_tag, condition) %>%
      mutate(top_five = TRUE) %>%
      right_join(median_melted_results) %>%
      mutate(top_five = case_when(is.na(top_five) ~ FALSE, TRUE ~ TRUE))

median_melted_results <- median_melted_results %>% 
  mutate(top_five = case_when((FDR <= 0.05 & gene %in% c("pgsA")) | top_five == TRUE ~ TRUE, TRUE ~ FALSE))

# median_melted_results <- median_melted_results %>%
# #modify to bold if the name is a locus tag and matches the format PA14_*, by replacing the word "italic" with "bold"
#   mutate(gene_name_stylized = case_when(
#     grepl("^PA14_", gene_name) ~ paste0("**", gene_name, "**"),
#     TRUE ~ gene_name
#   ))


volcano_plot <- median_melted_results %>%
  filter(gene != "control") %>%
  mutate(condition = case_when(
    condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
    condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
    condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
  )) %>%
  ggplot(aes(x = medLFC, y = FDR)) +
  geom_point(
    alpha = 0.5,
    shape = 20,
    size = 4,
    aes(
      alpha = FDR <= 0.05 & abs(medLFC) >= 1,
      color = top_five
    )
  ) +
  facet_wrap(facets = c("condition")) +
  doc_theme +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", lwd = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", lwd = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", lwd = 0.75) +
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
      label = gene_name_stylized,
      color = FDR < 0.05 & abs(medLFC) > 1
    ),
    alpha = 1
  ) +
  scale_alpha_manual(values = c(0.5, 1), guide = "none") +
  scale_color_manual(values = c("dark gray", "black"), guide = "none") +
  scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans())

print(volcano_plot)




# designed_targets <- fread("designed_targets.tsv", na.strings = c("None")) %>%
# mutate(type = case_when(
#   mismatches == 0 ~ "knockdown",
#   is.na(target) ~ "control"
# )) %>% filter(is.na(sp_dir) | sp_dir != tar_dir)



# median_spacers <- count_stats %>%
# filter(count != 0) %>%
# select(condition, sequence, count) %>%
#   # rbind(freezer_stock) %>% inner_join(targets %>% select(spacer, target) %>% unique()) %>% 
#   group_by(condition) %>%
#   mutate(cpm = cpm(count), lcpm = log(cpm)) %>% 
#   group_by(sequence) %>%
#   mutate(LMT_count = log(mean(cpm)), t_deviance = lcpm - LMT_count, t_sd = sd(cpm)) %>% 
#   filter(abs(t_deviance) == min(abs(t_deviance)) & t_sd == min(t_sd)) %>%
#   inner_join(designed_targets %>% rename(sequence = spacer) %>% filter(sp_dir != tar_dir & overlap == 20 & offset >= 0)) %>%
#   group_by(locus_tag) %>%
#   mutate(LMG_count = log(mean(cpm)), g_deviance = t_deviance - LMG_count, g_sd = sd(cpm)) %>%
#   filter(abs(g_deviance) == min(abs(g_deviance)) & g_sd == min(g_sd)) %>%
#   ungroup %>%
#   select(locus_tag, sequence) %>% unique()

# median_spacers %>% inner_join(non_normalized_melted_results) %>%
#   filter(gene != "control") %>%
#   mutate(condition = case_when(
#     condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
#     condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
#     condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
#   )) %>%
#   ggplot(aes(x = FDR, y = LFC)) +
#   geom_point(
#     aes(
#       alpha = FDR <= 0.05 & abs(LFC) >= 1,
#       color = gene %in% genes_of_interest
#     )
#   ) +  
#   ggrepel::geom_label_repel(
#     data = . %>%
#       group_by(condition) %>%
#       filter(gene != "control") %>%
#       arrange(FDR) %>%
#       mutate(index = seq_len(n())) %>%
#       filter(gene %in% genes_of_interest | index <= 30),
#     min.segment.length = 0,
#     parse = TRUE,
#     max.overlaps = Inf,
#     aes(
#       label = gene,
#       color = FDR < 0.05 & abs(LFC) > 1
#     ),
#     alpha = 1
#   ) +

#   facet_wrap(facets = c("condition")) +
#   theme(legend.position = "bottom") +
#   geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +
#   geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "blue") +
#   scale_alpha_manual(values = c(0.5, 1), guide = "none") +
#   scale_color_manual(values = c("dark gray", "red"), guide = "none") +
#   scale_x_continuous(trans = scales::reverse_trans() %of% scales::log10_trans())

