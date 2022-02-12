require('pacman')

p_load(
	ggthemr,
	data.table,
	scales,
	edgeR,
	pheatmap,
	svglite,
	ggplot2,
	ggrepel,
	colourpicker,
	RColorBrewer,
	poolr
)


annotated_key <- fread("annotated_key.tsv")

exp_design <- fread("exp_design.tsv")
exp_design <- exp_design[condition != "Undetermined"]

inoculum_exp_design <- 
	data.table(
		condition = "Inoculum",
		media = "inoculum",
		gDNA_source = "pellet",
		growth_condition = "t0",
		rep = 2,
		generations = 0)

exp_design <- rbind(exp_design, inoculum_exp_design)

exp_design[, verbose := paste(media, gDNA_source, growth_condition, rep, sep = "_")]

setorder(exp_design, condition)

all_counts <- 
	fread(
		"all_counts_seal.tsv", 
		header = FALSE, 
		col.names = c(
			"promoter", 
			"spacer", 
			"count", 
			"condition"))

inoculum_counts <- 
	fread(
		"inoculum_counts.tsv", 
		header = FALSE, 
		col.names = c(
			"promoter", 
			"spacer", 
			"count"))

inoculum_counts[, condition := "Inoculum"]

all_counts <- rbind(all_counts, inoculum_counts)

all_counts <- all_counts[promoter == "P1"][, .(spacer, condition, count)]

setorder(all_counts, condition)

setorder(exp_design, condition)

################################################################################
# Check for Data Integrity
################################################################################

data_grid <- 
	data.table::dcast(
		all_counts,
		spacer ~ condition,
		value.var = "count",
		fill = 0,
		fun.aggregate = sum)

data_grid_remelted <-
	melt(
		data_grid,
		variable.name = "condition",
		value.name = "count",
		id.vars = c('spacer'))

# print(
# 	data_grid_remelted[, .(
# 		median_count = median(count)),
# 		by = .(condition)][exp_design, 
# 											 on = .(condition)])

################################################################################

data_grid_matrix <- data.matrix(data_grid[, -c("spacer")])

row.names(data_grid_matrix) <- data_grid$spacer

crossjoin_correlation_grid <- cor(data_grid_matrix)

# Create a square matrix from the list of pairwise condition correlations.

################################################################################

all_counts <-
	all_counts[condition != "inoculum"][, .(spacer, condition, count)]

plot_matrix <- crossjoin_correlation_grid

break_halves <- length(unique(as.vector(plot_matrix)))

breaks <- c(
	seq(min(plot_matrix),
			median(plot_matrix),
			length.out = break_halves)[-break_halves],
	seq(median(plot_matrix),
			max(plot_matrix),
			length.out = break_halves))

breaks <- breaks[-length(breaks)]

breaks <- c(breaks, 0.99999999)

plot_colors <-
	c(colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves],
		colorRampPalette(c("white", "#007ac1"))(break_halves)[-c(1, break_halves)])

plot_colors <-
	colorRampPalette(c("white", "#007ac1"))(break_halves * 2 - 1)[-c(1, break_halves)]

plot_colors <- c(plot_colors, "dark grey")

to_plot_title <- paste("Raw Count Condition Correlations")

to_plot <- pheatmap(
	plot_matrix,
	col = plot_colors,
	breaks = breaks,
	border_color = NA,
	cellwidth = 20,
	cellheight = 20,
	main = to_plot_title,
	angle_col = 315,
	# fontsize_col = 10,
	# fontsize_row = 10,
	# cluster_cols = FALSE,
	show_rownames = TRUE,
	show_colnames = TRUE,
	clustering_method = "ward.D2",
	# clustering_distance_rows = "maximum",
	# clustering_distance_cols = "maximum"
)

print(to_plot)

################################################################################

data_group <-
	factor(
		exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")],
		levels = unique(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")]))

data_permut <- model.matrix(~ 0 + data_group)

colnames(data_permut) <- levels(data_group)

rownames(data_permut) <- colnames(data_grid_matrix)

data_permut_check <-
	melt(
		data.table(
			data_permut, 
			keep.rownames = "condition"), 
		id.vars = "condition")[value == 1][, .(condition, variable)]

data_permut_check <-
	data_permut_check[exp_design, on = .(condition == condition)]

print(data_permut_check)

data_y <- DGEList(
	counts = data_grid_matrix,
	group = data_group,
	genes = row.names(data_grid_matrix)
)

data_keep <- filterByExpr(data_grid_matrix, data_group, large.n = 1000)

data_y <- data_y[data_keep, , keep.lib.sizes = FALSE]

data_y <- calcNormFactors(data_y)

data_y <- estimateDisp(data_y, data_permut)

data_fit <- glmQLFit(data_y, data_permut, robust = TRUE)

data_CPM <- cpm(data_y, prior.count = 0)

################################################################################
# CPM Heatmap

plot_grid <- cor(data_CPM)

plot_matrix <- data.matrix(plot_grid)

break_halves <- length(unique(as.vector(plot_matrix)))

breaks <- c(
	seq(min(plot_matrix),
			median(plot_matrix),
			length.out = break_halves)[-break_halves],
	seq(median(plot_matrix),
			max(plot_matrix),
			length.out = break_halves)
)

breaks <- breaks[-length(breaks)]
breaks <- c(breaks, 0.99999999)

plot_colors <-
	colorRampPalette(c("white", "#007ac1"))(break_halves * 2 - 1)[-c(1, break_halves)]

plot_colors <- c(plot_colors, "dark grey")

to_plot_title <-
	paste("Clustering of Conditional Correlations (CPM)")

to_plot <- pheatmap(
	plot_matrix,
	col = plot_colors,
	breaks = breaks,
	border_color = NA,
	cellwidth = 20,
	cellheight = 20,
	main = to_plot_title,
	angle_col = 315,
	# fontsize_col = 10,
	# fontsize_row = 10,
	# cluster_cols = FALSE,
	show_rownames = TRUE,
	show_colnames = TRUE,
	clustering_method = "ward.D2",
	clustering_distance_rows = "maximum",
	clustering_distance_cols = "maximum"
)

print(to_plot)

################################################################################
# label with what the sequencing experiment actually is

to_plot <- pheatmap(
	plot_matrix,
	col = plot_colors,
	breaks = breaks,
	border_color = NA,
	cellwidth = 20,
	cellheight = 20,
	main = to_plot_title,
	angle_col = 315,
	# fontsize_col = 10,
	# fontsize_row = 10,
	# cluster_cols = FALSE,
	show_rownames = TRUE,
	show_colnames = TRUE,
	clustering_method = "ward.D2",
	labels_row = factor(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")]),
	labels_col = factor(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")]),
	clustering_distance_rows = "maximum",
	clustering_distance_cols = "maximum"
)

print(to_plot)

################################################################################

data_CPM_by_group <- 
	copy(data_CPM)

colnames(data_CPM_by_group) <- 
	factor(exp_design[,  paste(
		media, gDNA_source, growth_condition, sep = "_")])

print(to_plot)

################################################################################
# other diagnostic plots

plotMDS(log2(data_CPM_by_group))
plotQLDisp(data_fit)
plotBCV(data_y)

################################################################################

contrast_levels <-
	c(
		"mouse_plated_10x_inoculum_dilution - inoculum_plated_t0",
		"mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0",
		"inoculum_plated_t0 - inoculum_pellet_t0"
	)

data_contrast <- makeContrasts(
	contrasts = contrast_levels,
	levels = data_permut)

################################################################################

results_FDR <- all_counts[, .(genes = unique(spacer))]

results_LFC <- all_counts[, .(genes = unique(spacer))]

################################################################################

for (i in 1:ncol(data_contrast)) {
	
	results <- glmQLFTest(data_fit, contrast = data_contrast[, i])
	
	results <- topTags(results, n = Inf)
	
	results <- data.table(results$table)
	
	print(paste("Processing results for", contrast_levels[i], "..."))
	
	results_FDR <- results[, .(genes, FDR)][results_FDR, on = .(genes)]
	
	setnames(results_FDR, "FDR", contrast_levels[i])
	
	results_LFC <-
		results[, .(genes, logFC)][results_LFC, on = .(genes)]
	
	setnames(results_LFC, "logFC", contrast_levels[i])
}

################################################################################

melted_results_FDR <-
	data.table::melt(
		results_FDR,
		id.vars = c("genes"),
		variable.name = "condition",
		value.name = "FDR",
		measure.vars = contrast_levels
	)

melted_results_FDR <- melted_results_FDR[!is.na(FDR)]

################################################################################

melted_results_LFC <-
	data.table::melt(
		results_LFC,
		id.vars = c("genes"),
		variable.name = "condition",
		value.name = "LFC",
		measure.vars = contrast_levels
	)

################################################################################

melted_results <-
	melted_results_LFC[
		melted_results_FDR,
		on = .(genes, condition)]

melted_results <- melted_results[!is.na(FDR)]

melted_results[
	, LFC := melted_results[
		i  = genes %like% "Ctrl",
		j  = .(ctrl_medLFC = median(LFC, na.rm = TRUE)),
		by = .(condition)][
			i  = .SD,
			on = .(condition),
			j  = .(adj_medLFC = LFC - ctrl_medLFC),
			by = .EACHI]$adj_medLFC]

################################################################################

melted_results <-
	annotated_key[melted_results, on = .(name == genes)]

################################################################################

median_melted_results <-
	melted_results[, .(
		medLFC = median(LFC),
		FDR = stouffer(FDR)$p),
		by = .(locus_tag, gene_name, type, condition)]

################################################################################

setorder(median_melted_results, FDR)

################################################################################
# add fancy names to median melted results

median_melted_results[gene_name != ".", gene_name_stylized := paste0("italic('", gene_name, "')")]
median_melted_results[gene_name == ".", gene_name_stylized := paste0("bold('", locus_tag, "')")]
median_melted_results[gene_name == "control", gene_name_stylized := paste0("bold('", locus_tag, "')")]

################################################################################
################################################################################
################################################################################
# USE DESCRIPTIVE TERMINOLOGY TO DESCRIBE THE GROUPS OF SAMPLES
# EVENTUALLY THIS NEEDS TO BE MOVED INTO A NEW FILE
# PERHAPS STICK IT INTO A SIMPLE DEFINITIONS FILE
# PERHAPS 2 COLUMNS DEFINING WHICH SAMPLES BELONG TO GROUPS
################################################################################
################################################################################
################################################################################
# melt CPM for later use

CPM_melted <- melt(
	data.table(
		data_CPM, 
		keep.rownames = "spacer"), 
	id.vars = "spacer", 
	variable.name = "condition", 
	value.name = "CPM")

grouped_CPM <- copy(CPM_melted)

grouped_CPM[
	condition %in% c("Inoculum", "Mouse_P1_003"),  
	Condition := 'Pelleted inoculum t_0']

# grouped_CPM[
# 	condition == "Mouse_P1_006",
# 	Condition := 'Plated inoculum t_0']

grouped_CPM[
	condition %in% c("Mouse_P1_015", "Mouse_P1_016", "Mouse_P1_017"),  
	Condition := 'Plated ex-vivo 10× dilution']

grouped_CPM[
	condition %in% c("Mouse_P1_018", "Mouse_P1_019", "Mouse_P1_020", "Mouse_P1_021", "Mouse_P1_022"),  
	Condition := 'Plated ex-vivo 100× dilution']

grouped_CPM[
	condition %in% c("Mouse_P1_007", "Mouse_P1_008"),  
	Condition := 'Pelleted ex-vivo 10× dilution']

grouped_CPM <- exp_design[grouped_CPM, on = .(condition)]

grouped_CPM <- annotated_key[grouped_CPM, on = .(name == spacer)]

grouped_CPM[!is.na(rep), verbose := paste(media, gDNA_source, growth_condition, rep, sep = "_")]

grouped_CPM[is.na(rep), verbose := paste(media, gDNA_source, growth_condition, sep = "_")]


################################################################################