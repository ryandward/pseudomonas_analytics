source("cleanup_analysis.R")

p_load(tidyverse, ggplot2, data.table, ggrepel, hrbrthemes, viridis, ggallin)

doc_theme <- theme_ipsum(
	base_family = "Arial", 
	caption_margin = 12,
	axis_title_size = 12,
	axis_col = "black")


# top_ten <- median_melted_results %>% arrange(FDR) %>% select(gene) %>% filter(gene != "control") %>% unique %>% head(10)

median_melted_results %>% filter(gene != "control") %>% 	
	mutate(condition = case_when(
		condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
		condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
		condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro")) %>% 
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
			mutate(index = 1:n()) %>% 
			# filter((index <=10 & FDR < 0.05) | gene %in% c("pgsA", "orfN", "cysS")),
			filter(gene %in% c("pgsA", "orfN", "cysS")),
		# filter(gene %in% top_ten$gene & FDR < 0.05),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf,
		aes(label = gene)) +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans())

################################################################################

melted_results %>% filter(gene != "control") %>% 	
	mutate(condition = case_when(
		condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
		condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
		condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro")) %>% 
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
			mutate(index = 1:n()) %>% 
			filter(index <=10 & FDR < 0.05),
		# filter(gene %in% top_ten$gene & FDR < 0.05),
		min.segment.length = 0,
		max.overlaps = Inf,
		aes(label = guide)) +
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

median_melted_results %>% 
	filter(gene != "control") %>% 	
	mutate(condition = case_when(
		condition == "plated_10x_inoculum_dilution_mouse - plated_t0_inoculum" ~ "In vivo",
		condition == "plated_10x_inoculum_dilution_mouse - plated_6_generations_LB" ~ "In vivo v. in vitro",
		condition == "plated_6_generations_LB - plated_t0_inoculum" ~ "In vitro"
	)) %>% 
	ggplot(aes(x = medLFC, y = FDR)) + 
	geom_point(
		aes(
			alpha = FDR <= 0.05 & abs(medLFC) >= 1,
			color = gene %in% genes_of_interest)
		) + 
	facet_wrap(facets = c("condition")) + 
	doc_theme + 
	theme(legend.position = "bottom") +
	geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +
	geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
	geom_label_repel(
		data = . %>% 
			group_by(condition) %>% 
			filter(gene != "control") %>% 
			arrange(FDR) %>% 
			mutate(index = 1:n()) %>% 
			filter(gene %in% genes_of_interest),
		min.segment.length = 0,
		parse = TRUE,
		max.overlaps = Inf,
		aes(
			label = gene_name_stylized,
			color = FDR < 0.05 & abs(medLFC) > 1),
		alpha =1
	) +
	scale_alpha_manual(values = c(0.5, 1), guide = "none") +
	scale_color_manual(values = c("dark gray", "red"), guide = "none") +
	scale_y_continuous(trans = scales::reverse_trans() %of% scales::log10_trans())


