median_melted_results %>% 
	mutate(
		Response = case_when(
			medLFC < 0 & FDR < 0.05 ~ "Vulnerable",
			medLFC > 0 & FDR < 0.05 ~ "Resistant",
			TRUE ~ "No Response")) %>%
	select(locus_tag, condition, Response) %>% 
	mutate(condition = case_when(
		condition == "mouse vs inoc" ~ "In vivo",
		condition == "mouse vs plate" ~ "In vivo vs. in vitro",
		condition == "plate vs inoc" ~ "In vitro"))	%>%
	pivot_longer(!c(condition, locus_tag)) %>% 
	pivot_wider(id_cols = c(locus_tag), names_from = condition, values_from = value) %>%
	make_long(`In vitro`, `In vivo`, `In vivo vs. in vitro` ) %>% 
	ggplot(aes(
		x = x, 
		next_x = next_x,
		node = node,
		next_node = next_node,
		fill = factor(node),
		label = node)) +
	geom_sankey(
		flow.colour = "black",
		flow.alpha = 0.25, 
		node.color = "black") +
	scale_fill_viridis(discrete = T, direction = 1, option = "inferno", alpha = 0.55) +
	geom_sankey_label(size = 3.5, color = 1, fill = "white") +
	theme_sankey(base_size = 16) +
	guides(fill = guide_legend(title = "Relative Response to Knockdown")) +
	theme(axis.title.x = element_blank()) 

