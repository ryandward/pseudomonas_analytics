# source("Publication//publication_counter.R")

library(pacman)

p_load(
	data.table,
	scales,
	edgeR,
	pheatmap,
	svglite,
	ggplot2,
	ggrepel,
	colourpicker,
	RColorBrewer,
	poolr,
	statmod
)

purine_enrichment <- fread("KW-0658.tsv", header = FALSE, col.names = c("gene"))


#############################################################################
# prepare data for volcano plots
for (this.contrast in unique(median_melted_results$condition)) {
	

	to_plot <- median_melted_results[condition == this.contrast]
	
	setorder(to_plot, FDR)
	
	to_plot[FDR < 0.05, Significance := "FDR < 0.05"]
	to_plot[FDR < 0.005, Significance := "FDR < 0.005"]
	to_plot[FDR < 0.0005, Significance := "FDR < 0.0005"]
	
	plot_object <-
		ggplot(
			data = to_plot,
			aes(x = medLFC, y = -log10(FDR))) +
		geom_point(
			aes(color = `Significance`),
			size = 2) +
		xlim(
			median_melted_results[type != "control", min(medLFC)], median_melted_results[type != "control", max(medLFC)]) +
		ylim(
			median_melted_results[type != "control", min(-log10(FDR))], median_melted_results[type != "control", max(-log10(FDR))]) +
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
			data = to_plot[gene %like% "isp" | gene %like% "orf" | medLFC < -12],
			aes(label = gene),
			size = 5,
			box.padding = unit(0.5, "lines"),
			point.padding = unit(0.5, "lines"),
			max.iter = 500,
			max.overlaps = 100,
			parse = TRUE)

	print(plot_object)
}

