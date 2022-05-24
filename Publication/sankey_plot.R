source("Publication/publication_counter.R")

median_melted_results %>% 
	mutate(
	Response = case_when(
		medLFC < -4 & FDR < 0.05 ~ "LFC < -4",
		medLFC < -2 & FDR < 0.05 ~ "LFC < -2",
		medLFC < -1 & FDR < 0.05 ~ "LFC < -1",
		TRUE ~ "No Response")) %>%
	mutate(condition = case_when(
		condition == "mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0" ~ "Lung (v. Inoculum)",
		condition == "mouse_plated_10x_inoculum_dilution - LB_plated_6_generations" ~ "Lung (v. Plate)",
		condition == "LB_plated_6_generations - inoculum_pellet_t0" ~ "Plate (v. Inoculum)"))	%>%
	select(locus_tag, condition, Response) %>% 
	pivot_longer(!c(condition, locus_tag)) %>% 
	pivot_wider(id_cols = c(locus_tag), names_from = condition, values_from = value) %>%
	make_long( `Plate (v. Inoculum)`, `Lung (v. Plate)`, `Lung (v. Inoculum)`) %>% 
	ggplot(aes(
		x = x, 
		next_x = next_x,
		node = node,
		next_node = next_node,
		fill = factor(node),
		label = node)) +
	geom_sankey(
		flow.colour = "black",
		flow.alpha = 0.5, 
		node.color = "black") +
	scale_fill_viridis(discrete = T) +
	geom_sankey_label(size = 3.5, color = 1, fill = "white") +
	theme_sankey(base_size = 16) +
	guides(fill = guide_legend(title = "Response to Knockdown")) +
	theme(axis.title.x = element_blank()) 

