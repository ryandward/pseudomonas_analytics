median_melted_results %>% filter(gene != "control") %>% 	
	mutate(condition = case_when(
		condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
		condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
		condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro")) %>% 
	ggplot(aes(x = medLFC, y = FDR)) + 
	geom_point(aes(colour = FDR < 0.01 & abs(medLFC) > 1)) + 
	facet_wrap(facets = c("condition")) + 
	theme_bw() + 
	theme(legend.position = "bottom") +
	scale_color_manual(values = c("gray", "red"), na.value = "grey") +
	geom_label_repel(
		data = . %>% filter(FDR < 1e-30),
		aes(label = gene)) +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans())