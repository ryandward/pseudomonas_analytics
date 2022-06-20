all_counts %>% inner_join(exp_design) %>% inner_join(annotations) %>%
	mutate(sample_group = case_when(
		sample_group %like% "mating" ~ "Mating Strain",
		sample_group %like% "^inoculum_pellet" ~ "PA14 Library")) %>% arrange(promoter) %>% pivot_wider(id_cols = c(gene, locus_tag, sequence, sample_name), names_from = c(promoter), values_from = count) %>% View

all_counts %>% inner_join(exp_design) %>% inner_join(annotations) %>% 	filter(sample_group %like% "mating" | condition == "P1_PA14") %>% pivot_wider(id_cols = c(sequence, gene, promoter), names_from = condition, values_from = count, values_fill = 0) %>% mutate(diff = P1_mfdpir - P1_PA14) %>% arrange(diff) %>% View

