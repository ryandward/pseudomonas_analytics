# https://www.nature.com/articles/nmeth.3253

source("cleanup_analysis.R")

botneck_stats <- fread("stats.tsv", col.names = c("sequence", "count", "condition"))

botneck_stats <- botneck_stats %>%
  inner_join(all_guides) %>%
  filter(type != "focused")

exp_design_botneck <- fread("unified_counts/exp_design.tsv") %>% mutate(sample_group = paste(gDNA_source, growth_condition, media, sep = "_"))


botneck.t0 <- botneck_stats %>%
  inner_join(all_guides) %>%
  inner_join(exp_design_botneck) %>%
  filter(condition %in% c("dJMP1", "dJMP2")) %>%
  group_by(sequence, type) %>%
  summarise(count = median(count)) %>%
  group_by(type) %>%
  mutate(
    fi0 = count / sum(count),
    count0 = count
  ) %>%
  select(type, sequence, fi0, count0) %>%
  nest() %>%
  rename(data0 = data) %>%
  mutate(s0 = map_dbl(data0, ~ sum(.$count0)))

botneck <- botneck_stats %>%
  inner_join(all_guides) %>%
  inner_join(exp_design_botneck) %>%
  filter(generations != 0) %>%
  group_by(type, condition, generations) %>%
  mutate(
    fis = count / sum(count)
  ) %>%
  nest() %>%
  mutate(
    ss = map_dbl(data, ~ sum(.$count))
  ) %>%
  full_join(botneck.t0) %>%
  mutate(data = map2(data, data0, inner_join))

botneck <- botneck %>%
  mutate(
    data = map(
      data, ~ .x %>%
        mutate(
          ratio = ((fis - fi0)^2) / (fi0 * (1 - fi0)^2)
        )
    )
  )

botneck <- botneck %>%
  mutate(
    f_hat = map_dbl(
      data,
      ~ sum(.$ratio)
    ) * (1 / map_dbl(data, ~ n_distinct(.$sequence))),
    Nb = generations / (f_hat - 1 / s0 - 1 / ss)
  )

botneck.stats <- botneck %>%
  inner_join(exp_design_botneck) %>%
  group_by(sample_group, type) %>%
  summarise(
    Nb.med = median(Nb),
    Nb.range = max(Nb) - min(Nb),
    Nb.mean = mean(Nb),
    Nb.sd = sd(Nb)
  )

botneck.stats <- botneck.stats %>% mutate(Nb.sd = case_when(is.na(Nb.sd) ~ 0, TRUE ~ Nb.sd))

botneck %>%
  inner_join(exp_design_botneck) %>%
  filter(!media %like% "Gent") %>%
  filter(gDNA_source == "plated") %>%
  filter(!condition %in% c("dJMP3", "dJMP5")) %>%
  # change the name plated_100x_inoculum_dilution_mouse to "100x dilution (mouse)"
  # change the name plated_10x_inoculum_dilution_mouse to "10x dilution (mouse)"
  # change the name plated_6_generations_LB to "6 generations (LB)"
  mutate(sample_group = case_when(
    sample_group == "plated_100x_inoculum_dilution_mouse" ~ "100x dilution (mouse)",
    sample_group == "plated_10x_inoculum_dilution_mouse" ~ "10x dilution (mouse)",
    sample_group == "plated_6_generations_LB" ~ "6 generations (LB)"
  )) %>%
  ggplot(aes(x = sample_group, y = Nb, fill = type)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Population Complexity: Bottleneck Metric") +
  scale_y_continuous(
    # trans = scales::pseudo_log_trans(base = 10),
    # breaks = c(0, 10^seq(1, 7)),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  theme(axis.title.x = element_blank()) +
  # make control guides grey, and knockdown guides the "red" from the paired palette
  scale_fill_manual(values = c("grey", "red"))


# botneck_stats %>%
#   filter(count != 0) %>%
#   filter(type == "knockdown") %>%
#   group_by(locus_tag, condition) %>%
#   summarize(unique_count = n()) %>%
#   data.table() %>%
#   dcast(locus_tag ~ condition, value.var = "unique_count", fill = 0) %>%
#   melt(value.name = "unique_count", variable.name = "condition") %>%
#   group_by(
#     condition, unique_count
#   ) %>%
#   summarize(count_count = n()) %>%
#   data.table() %>%
#   dcast(condition ~ unique_count, value.var = "count_count", fill = 0) %>%
#   i()
# nner_join(exp_design_botneck) %>%
#   filter(!condition %in% c("dJMP3", "dJMP5")) %>%
#   filter(!media %like% "Gent") %>%
#   filter(gDNA_source == "plated") %>%
#   clipr::write_clip()

botneck_stats %>%
  filter(count != 0) %>%
  filter(type == "knockdown") %>%
  rbind(
    botneck_stats %>%
      filter(count != 0) %>%
      filter(type == "control") %>%
      filter(condition == "P1_mfdpir")
  ) %>%
  group_by(locus_tag, condition) %>%
  summarize(unique_count = n()) %>%
  data.table() %>%
  dcast(locus_tag ~ condition, value.var = "unique_count", fill = 0) %>%
  melt(value.name = "unique_count", variable.name = "condition") %>%
  group_by(
    condition, unique_count
  ) %>%
  summarize(count_count = n()) %>%
  data.table() %>%
  dcast(condition ~ unique_count, value.var = "count_count", fill = 0) %>%
  inner_join(exp_design_botneck) %>%
  filter(!condition %in% c("dJMP3", "dJMP5")) %>%
  filter(!media %like% "Gent") %>%
  filter(gDNA_source == "plated" | condition == "P1_mfdpir") %>%
  mutate(growth_condition = case_when(
    condition == "P1_mfdpir" ~ "mating strain",
    TRUE ~ growth_condition
  )) %>%
  select(growth_condition, media, `0`, `1`, `2`, `3`, `4`) %>%
  data.table() %>%
  melt(id.vars = c("growth_condition", "media"), variable.name = "guides_per_gene", value.name = "times_seen") %>%
  group_by(growth_condition, media, guides_per_gene) %>%
  mutate(rep = row_number()) %>%
  ggplot(aes(fill = guides_per_gene, y = times_seen, x = rep)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.9),
    # position = "stack",
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  facet_grid(. ~ growth_condition)



stat_summary(
  fun = mean,
  geom = "bar",
  position = position_dodge(width = 0.9),
  # position = "stack",
  color = "black"
) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.9),
    width = 0.2
  )


range_fun <- function(y) {
  ymin <- min(y)
  ymax <- max(y)
  return(c(ymin = ymin, ymax = ymax))
}

botneck_stats %>%
  filter(count != 0) %>%
  filter(type == "knockdown") %>%
  rbind(
    botneck_stats %>%
      filter(count != 0) %>%
      filter(type == "control") %>%
      filter(condition == "P1_mfdpir")
  ) %>%
  group_by(locus_tag, condition) %>%
  summarize(unique_count = n()) %>%
  data.table() %>%
  dcast(locus_tag ~ condition, value.var = "unique_count", fill = 0) %>%
  melt(value.name = "unique_count", variable.name = "condition") %>%
  group_by(
    condition, unique_count
  ) %>%
  summarize(count_count = n()) %>%
  data.table() %>%
  dcast(condition ~ unique_count, value.var = "count_count", fill = 0) %>%
  inner_join(exp_design_botneck) %>%
  filter(!condition %in% c("dJMP3", "dJMP5")) %>%
  filter(!media %like% "Gent") %>%
  filter(gDNA_source == "plated" | condition == "P1_mfdpir") %>%
  mutate(growth_condition = case_when(
    condition == "P1_mfdpir" ~ "mating strain",
    TRUE ~ growth_condition
  )) %>%
  select(growth_condition, media, `0`, `1`, `2`, `3`, `4`) %>%
  data.table() %>%
  melt(id.vars = c("growth_condition", "media"), variable.name = "guides_per_gene", value.name = "times_seen") %>%
  filter(media == "inoculum") %>%
  group_by(guides_per_gene) %>%
  ggplot(aes(x = guides_per_gene, y = times_seen)) +
  stat_summary(
    fun = median,
    geom = "bar",
    position = position_dodge(width = 0.9),
    # position = "stack",
    color = "NA",
    alpha = 0.5
  ) +
  stat_summary(
    fun.data = range_fun,
    geom = "errorbar",
    position = position_dodge(width = 0.9),
    width = 0.2,
    # lwd = 1
  ) +
  # label x axis as "Guides per Gene"
  xlab("Guides per Gene") +
  ylab("Times Seen (Median and Range)") +
  theme_minimal()
