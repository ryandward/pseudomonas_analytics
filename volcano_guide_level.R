source("deep_integration_analysis.R")

purine_enrichment <- fread("KW-0658.tsv", header = FALSE, col.names = c("gene_name"))

#############################################################################
# prepare data for volcano plots

this.contrast <- "mouse vs plate0"

to_plot <- melted_results[condition == this.contrast]

setorder(to_plot, FDR)

to_plot[FDR < 0.01, Significance := "FDR < 0.01"]
to_plot[FDR < 0.001, Significance := "FDR < 0.001"]
to_plot[FDR < 0.0001, Significance := "FDR < 0.0001"]

to_plot[gene_name %in% purine_enrichment[, gene_name], Pathway := "Purine Biosynthesis"]

to_plot[!gene_name %in% purine_enrichment[, gene_name], index := .I]

to_plot[index <= 10, `Hits` := "Top ten (Non-purine)"]

################################################################################
# First plot. Colors for FALSE DISCOVERY RATE GRADIENTS
# 
# 
# plot_object <- 
# 	ggplot(
# 		data = to_plot, 
# 		aes(x = LFC, y = -log10(FDR))) +
# 	geom_point(
# 		aes(color = Significance), size = 2) +
# 	xlim(
# 		median_melted_results[type != "control", min(LFC)], median_melted_results[type != "control", max(LFC)]) +
# 	ylim(
# 		median_melted_results[type != "control", min(-log10(FDR))], median_melted_results[type != "control", max(-log10(FDR))]) +
# 	theme_bw(
# 		base_size = 12) +
# 	ggtitle(this.contrast) +
# 	theme(
# 		plot.title = element_text(size = 16),
# 		axis.text = element_text(size = 14, color = "black"),
# 		axis.title = element_text(size = 14, color = "black"),
# 		legend.text = element_text(size = 14, color = "black"),
# 		legend.title = element_text(size = 14, color = "black"),
# 		legend.position = "bottom") + 
# 	geom_hline(
# 		yintercept = 2, 
# 		linetype = "dashed", 
# 		color = "red")
# 
# print(plot_object)
# 
# ################################################################################
# # Second plot. Pathway i.e. purine enrichment
# 
# plot_object <- 
# 	ggplot(
# 		data = to_plot, 
# 		aes(x = LFC, y = -log10(FDR))) +
# 	geom_point(
# 		aes(color = Pathway), 
# 		size = 2) +
# 	xlim(
# 		median_melted_results[type != "control", min(LFC)], median_melted_results[type != "control", max(LFC)]) +
# 	ylim(
# 		median_melted_results[type != "control", min(-log10(FDR))], median_melted_results[type != "control", max(-log10(FDR))]) +
# 	theme_bw(
# 		base_size = 12) +
# 		ggtitle(this.contrast) +
# 			theme(
# 		plot.title = element_text(size = 16),
# 		axis.text = element_text(size = 14, color = "black"),
# 		axis.title = element_text(size = 14, color = "black"),
# 		legend.text = element_text(size = 14, color = "black"),
# 		legend.title = element_text(size = 14, color = "black"),
# 		legend.position = "bottom") + 
# 	geom_hline(
# 		yintercept = 1.30103, 
# 		linetype = "dashed", 
# 		color = "red") +
# 	geom_label_repel(
# 		data = to_plot[!is.na(Pathway)],
# 		aes(label = gene_name_stylized),
# 		size = 5,
# 		box.padding = unit(0.5, "lines"),
# 		point.padding = unit(0.5, "lines"),
# 		max.iter = 500,
# 		max.overlaps = 100,
# 		parse = TRUE)
# 
# print(plot_object)

################################################################################
# Third plot. Top ten non-purine hits.

plot_object <- 
	ggplot(
		data = to_plot, 
		aes(x = LFC, y = -log10(FDR))) +
	geom_point(
		aes(color = `Hits`), 
		size = 2) +
	xlim(
		melted_results[type != "control", min(LFC)], melted_results[type != "control", max(LFC)]) +
	ylim(
		melted_results[type != "control", min(-log10(FDR))], melted_results[type != "control", max(-log10(FDR))]) +
	theme_bw(
		base_size = 12) +
	ggtitle(this.contrast) +
	theme(
		plot.title = element_text(size = 16),
		axis.text = element_text(size = 14, color = "black"),
		axis.title = element_text(size = 14, color = "black"),
		legend.text = element_text(size = 14, color = "black"),
		legend.title = element_text(size = 14, color = "black"),
		legend.position = "bottom") + 
	geom_hline(
		yintercept = 1.30103, 
		linetype = "dashed", 
		color = "red") +
	geom_label_repel(
		data = to_plot[gene_name %like% "isp"],
		aes(label = gene_name),
		size = 5,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.5, "lines"),
		max.iter = 500,
		max.overlaps = 100,
		parse = TRUE)

print(plot_object)
