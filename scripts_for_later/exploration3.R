require('pacman'); 
p_load(data.table, edgeR, reshape2, digest, forcats)

CPM <- cpm(y, prior.count = 0)

CPMdt <- data.table(CPM)

invitroCPM <- CPMdt[, .SD, .SDcols = patterns("mouse")]

invitroCPM <- data.matrix(invitroCPM)

rownames(invitroCPM) <- rownames(CPM)

plot_matrix <- log10(invitroCPM + 1)

break_halves <- length(unique(as.vector(plot_matrix)))

breaks <- c(
	seq(min(plot_matrix), 
			median(plot_matrix), 
			length.out = break_halves)[-break_halves], 
	seq(median(plot_matrix), 
			max(plot_matrix), 
			length.out = break_halves))

plot_colors <- c(
	colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves], 
	colorRampPalette(c("white", "#007ac1"))(break_halves)[-1])

to_plot_title <- paste("Log10 Counts per Million")

to_plot <- pheatmap(plot_matrix,
										col = plot_colors,
										breaks = breaks,
										border_color = NA,
										# cellwidth = 15,
										# cellheight = 12,
										cutree_rows = 4,
										cutree_cols = 4,
										main = to_plot_title,
										angle_col = 315,
										# fontsize_col = 10,
										# fontsize_row = 10,
										# cluster_cols = FALSE,
										# cutree_cols = 10,
										clustering_method = "ward.D2",
										clustering_distance_rows = "maximum",
										clustering_distance_cols = "maximum"
)

