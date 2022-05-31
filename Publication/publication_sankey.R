p_load(ggsankey, tidyverse, viridis)

median_melted_results %>% 
	mutate(
		Response = case_when(
			medLFC < -1 & FDR < 0.05 ~ "Vulnerable",
			medLFC < 1 & FDR < 0.05 ~ "Resistant",
			
			TRUE ~ "No Response")) %>%
	select(locus_tag, condition, Response) %>% 
	mutate(condition = case_when(
		condition == "mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0" ~ "In vivo",
		condition == "mouse_plated_10x_inoculum_dilution - LB_plated_6_generations" ~ "In vivo v. in vitro",
		condition == "LB_plated_6_generations - inoculum_pellet_t0" ~ "In vitro"))	%>%
	pivot_longer(!c(condition, locus_tag)) %>% 
	pivot_wider(id_cols = c(locus_tag), names_from = condition, values_from = value) %>%
	make_long(`In vitro`, `In vivo`, `In vivo v. in vitro` ) %>% 
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
	# scale_fill_viridis(
	# 	discrete = T, direction = 1, option = "inferno", alpha = 0.55) +
	scale_fill_manual(values = c(
		"grey", "#63ccff", "#ff6f60")) +
	geom_sankey_label(aes(colour = "node"),
										size = 3.5, color = 1) +
	theme_sankey(base_size = 16) +
	guides(
		fill = guide_legend(title = "Relative Response to Knockdown")) +
	theme(
		axis.title.x = element_blank(),
		legend.position = "bottom")

