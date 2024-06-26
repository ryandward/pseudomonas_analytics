# source("deep_integration_analysis.R")

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
			aes(color = Significance),
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
			data = to_plot[gene_name %like% "isp" | gene_name %like% "orf" | (medLFC > 0 & FDR < 0.05)],
			aes(label = gene_name_stylized),
			size = 5,
			box.padding = unit(0.5, "lines"),
			point.padding = unit(0.5, "lines"),
			max.iter = 50000,
			max.overlaps = 100,
			parse = TRUE)
	
	print(plot_object)
}

