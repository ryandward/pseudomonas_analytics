source("publication_counter.R")

mixed <- median_melted_results %>% mutate(
	depleted = (medLFC < -1 & FDR < 0.05), 
	resistant = (medLFC > 1 & FDR < 0.05)) %>%
	mutate(condition = case_when(
		condition == "mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0" ~ "Lung v Inoculum",
		condition == "mouse_plated_10x_inoculum_dilution - LB_plated_6_generations" ~ "Lung v Plate",
		condition == "LB_plated_6_generations - inoculum_pellet_t0" ~ "Plate v Inoculum"))	%>%
select(locus_tag, condition, depleted, resistant) %>% 
	pivot_longer(!c(condition, locus_tag)) %>% 
	unite("condition", c("condition", "name"), sep = " ") %>% 
	pivot_wider(id_cols = c(locus_tag), names_from = condition, values_from = value)

mixed.matrix <- mixed %>%
	make_comb_mat(mode = "intersect")

p <- UpSet(
	mixed.matrix, 
	top_annotation = HeatmapAnnotation(
		"Intersection" = anno_barplot(
			comb_size(mixed.matrix),
			border = FALSE, 
			height = unit(8, "cm"),
			add_numbers = T,
			gp = gpar(fill = c("grey", rep("grey", 2), rep("grey", 6)), lty = "solid")), 
		show_annotation_name = FALSE),
	right_annotation = rowAnnotation(
		"Genetic Phenotypes" = anno_barplot(
			set_size(mixed.matrix),
			border = FALSE,
			gp = gpar(
				# fill = viridis(set_size(mixed.matrix)/2 %>% length, direction = -1, alpha = 0.5),
				# fill = c(rep("light blue", 2), rep("coral", 2)),
				fill = viridis(3, direction = -1, alpha = 0.35) %>% rep(each = 2), 
				col = c("red", "black", "red", "black"),
				lwd = 3,
				lty = c("dotted", "solid", "dotted", "solid")),
			width = unit(6, "cm"),
			add_numbers = T),
		annotation_name_gp = gpar(fontsize = 8)),
	set_order = set_name(mixed.matrix),
	row_names_gp = grid::gpar(fontsize = 10),
	comb_col = c("black", rep("black", 2), rep("black", 6)))

print(p)