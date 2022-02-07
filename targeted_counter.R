require('pacman')
p_load(ggthemr, Rtsne, data.table, scales, edgeR, statmod, poolr, pheatmap, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, RColorBrewer)

annotated_key <- fread("annotated_key.tsv")
exp_design <- fread("exp_design.tsv")

all_counts <- fread("all_counts_seal.tsv", header = FALSE, col.names = c("promoter", "spacer", "count", "condition"))

# get rid of all the other promoters
all_counts <- all_counts[promoter == "P1"][, .(spacer, condition, count)]

# # these conditions suck
rejected_conditions = c("Mouse_P1_018")
# 18?

all_counts <- all_counts[!(condition %in% rejected_conditions)]
exp_design <- exp_design[!(condition %in% rejected_conditions)]

setorder(all_counts, condition)
setorder(exp_design, condition)


# exp_design[media == "Gent", media := "LB"]

################################################################################
# Check for Data Integrity
################################################################################
data_grid <- data.table::dcast(
	all_counts, 
	spacer ~ condition,
	value.var = "count", 
	fill = 0)

data_grid_remelted <- melt(data_grid, variable.name = "condition", value.name = "count", id.vars = c('spacer'))

print(data_grid_remelted[, .(median_count = median(count)), by = .(condition)][exp_design, on = .(condition)])
################################################################################
# Create a correlation matrix for each condition in the remelted matrix.
# It has to be dcasted, then remelted to have zeroes put into their appropriate spots!
rm(crossjoin_correlation)

for(i in unique(data_grid_remelted$condition)){
	for(j in unique(data_grid_remelted$condition)){
		if(!(exists("crossjoin_correlation"))){
			crossjoin_correlation <- data.table(i, j, cor(data_grid_remelted[condition == i]$count, data_grid_remelted[condition == j]$count))
		}
		else{
			crossjoin_correlation <-
				rbind(crossjoin_correlation, data.table(i, j, cor(data_grid_remelted[condition == i]$count, data_grid_remelted[condition == j]$count)))
		}
	}
}
crossjoin_correlation_grid <- dcast(crossjoin_correlation, i~j, value.var = "V3")
# Create a square matrix from the list of pairwise condition correlations.
################################################################################
all_counts <- all_counts[condition != "inoculum"][, .(spacer, condition, count)]

plot_grid <- crossjoin_correlation_grid
plot_matrix <- data.matrix(plot_grid[,-1])
rownames(plot_matrix) <- plot_grid$i

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

plot_colors <- c(
	colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves],
	colorRampPalette(c("white", "#007ac1"))(break_halves)[-c(1, break_halves)])

plot_colors <- colorRampPalette(c("white", "#007ac1"))(break_halves*2-1)[-c(1,break_halves)]
plot_colors <- c(plot_colors, "dark grey")

to_plot_title <- paste("Raw Count Condition Correlations")

to_plot <- pheatmap(plot_matrix,
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
#/Check for Data Integrity in naming scheme
################################################################################

data_grid_matrix <- data.matrix(data_grid[, -c("spacer")])

row.names(data_grid_matrix) <- data_grid$spacer

data_group <- factor(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")],
										 levels = unique(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")]))

data_permut <-  model.matrix( ~ 0 + data_group)

colnames(data_permut) <- levels(data_group)

rownames(data_permut) <- colnames(data_grid_matrix)

data_permut_check <- melt(data.table(data_permut,keep.rownames = "condition"), id.vars = "condition")[value==1][, .(condition, variable)]
data_permut_check <- data_permut_check[exp_design, on = .(condition==condition)]

print(data_permut_check)

data_y <- DGEList(counts = data_grid_matrix, 
									group = data_group, 
									genes = row.names(data_grid_matrix))

# data_keep <- filterByExpr(data_y, data_permut)

data_keep <- rowSums(cpm(data_y) > 10) >= 5

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

plot_colors <- c(
	colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves],
	colorRampPalette(c("white", "#007ac1"))(break_halves)[-c(1, break_halves)])

plot_colors <- colorRampPalette(c("white", "#007ac1"))(break_halves*2-1)[-c(1,break_halves)]
plot_colors <- c(plot_colors, "dark grey")

to_plot_title <- paste("Clustering of Conditional Correlations (CPM)")

to_plot <- pheatmap(plot_matrix,
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
										clustering_method = "ward.D",
										clustering_distance_rows = "maximum",
										clustering_distance_cols = "maximum")

print(to_plot)

# label with what the sequencing experiment actually is
to_plot <- pheatmap(plot_matrix,
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
										clustering_method = "ward.D",
										labels_row = factor(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")]),
										labels_col = factor(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")]),
										clustering_distance_rows = "maximum",
										clustering_distance_cols = "maximum")

print(to_plot)


# colnames(data_CPM) <- factor(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")])

# / CPM heatmap

plotQLDisp(data_fit)
plotMDS(data_CPM)
plotBCV(data_y)

################################################################################

contrast_levels <- c("mouse_plated_10x_inoculum_dilution - inoculum_plated_t0",
										 "inoculum_plated_t0 - inoculum_pellet_t0")

data_contrast <- makeContrasts(contrasts = contrast_levels, 
															 levels = data_permut)

################################################################################

results_FDR <- all_counts[, .(genes = unique(spacer))]
results_LFC <- all_counts[, .(genes = unique(spacer))]

################################################################################
################################################################################

for (i in 1:ncol(data_contrast)){
	
	results <- glmQLFTest(data_fit, contrast = data_contrast[,i])
	results <- topTags(results, n = Inf)
	results <- data.table(results$table)
	
	print(paste("Processing results for", contrast_levels[i], "..."))
	
	results_FDR <- results[, .(genes, FDR)][results_FDR, on = .(genes)]
	setnames(results_FDR, "FDR", contrast_levels[i])
	
	results_LFC <- results[, .(genes, logFC)][results_LFC, on = .(genes)]
	setnames(results_LFC, "logFC", contrast_levels[i])}

################################################################################

melted_results_FDR <- 
	data.table::melt(
		results_FDR, 
		id.vars = c(
			"genes"),
		variable.name = "condition", 
		value.name = "FDR",
		measure.vars = contrast_levels)

melted_results_FDR <- melted_results_FDR[!is.na(FDR)]

################################################################################

melted_results_LFC <- 
	data.table::melt(
		results_LFC, 
		id.vars = c(
			"genes"),
		variable.name = "condition", 
		value.name = "LFC",
		measure.vars = contrast_levels)

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

melted_results <- annotated_key[melted_results, on = .(name == genes)]

################################################################################

guide_toss <- melted_results[, .(sumFDR = sum(FDR)), by = .(name, gene_name, locus_tag)]

setorder(guide_toss, sumFDR)

guide_rejected <- guide_toss[, .(maxsumFDR = max(sumFDR)), by = .(gene_name, locus_tag)]

guide_rejected <- guide_toss[guide_rejected, on = .(locus_tag, gene_name)][sumFDR == maxsumFDR, .(name)]

melted_results_adjusted <- melted_results[!(name %in% guide_rejected$name)]

################################################################################

median_melted_results <- melted_results_adjusted[, .(
	medLFC = median(LFC), 
	FDR = stouffer(FDR)$p), 
	by = .(locus_tag, gene_name, type, condition)]

################################################################################

setorder(median_melted_results, FDR)

################################################################################

plot(-log10(FDR)~LFC, data = melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0"], pch = 20, main = "Guide-level: mouse_plated_10x_inoculum_dilution - inoculum_plated_t0")
points(-log10(FDR)~LFC, data = melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0" & type == "control"], pch = 20, col = "green")

plot(-log10(FDR)~medLFC, data = median_melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0"], pch = 20, main = "Gene-level (Median guide): mouse_plated_10x_inoculum_dilution - inoculum_plated_t0")
abline(v = median_melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0" & type == "control"]$medLFC)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# control_CPM <- melt(data.table(data_CPM, keep.rownames = TRUE), id.vars = "rn")[rn %like% "Ctrl"]
# control_CPM <- control_CPM[variable == "mouse_plated_10x_inoculum_dilution"]
# setorder(control_CPM, rn)
# 
# all_counts[spacer %like% "Ctrl", .(median_ctrl = median(count)), by = .(condition)][data_permut_check, on = .(condition==rn)][, .(condition, rep, variable, median_ctrl)]
# 
# ################################################################################
# ################################################################################
# 
# mouse_grid_FDR <- dcast(median_melted_results, locus_tag + gene_name + type ~ condition, value.var = "FDR")
# mouse_grid_LFC <- dcast(median_melted_results, locus_tag + gene_name + type ~ condition, value.var = "medLFC")
# 
# ################################################################################
# ################################################################################
# 
# FDR_order <- median_melted_results[, .(mean_order = mean(FDR)), by = .(locus_tag, gene_name, type)]
# setorder(FDR_order, mean_order, locus_tag)
# FDR_order[, FDR := .I]
# FDR_order[, medLFC := .I]
# FDR_order <- FDR_order[, .(locus_tag, gene_name, type, FDR, medLFC)]
# FDR_order[, condition := "order"]
# 
# # glue the dummy variable onto the results matrix
# # median_melted_results <- rbind(FDR_order, median_melted_results)
# 
# ################################################################################
# # construct the results then copy them to clipboard
# 
# to_clip <- dcast(rbind(FDR_order, median_melted_results), locus_tag + gene_name + type ~ condition, value.var = "FDR")
# setorder(to_clip, order)
# to_clip <- to_clip[, -c("order")]
# clipr::write_clip(to_clip)
# 
# to_clip <- dcast(rbind(FDR_order, median_melted_results), locus_tag + gene_name + type ~ condition, value.var = "medLFC")
# setorder(to_clip, order)
# to_clip <- to_clip[, -c("order")]
# clipr::write_clip(to_clip)
# 
# plot(-log10(FDR) ~ medLFC, data = median_melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0"])

##################################
# Chi Square Test, did it work?
##################################

did_it_work <- data.table(data_CPM, keep.rownames = "spacer")[, .(spacer, Mouse_P1_003, Mouse_P1_015, Mouse_P1_016, Mouse_P1_017)]

did_it_work[, t0 := Mouse_P1_003]
did_it_work[, tf := Mouse_P1_015 + Mouse_P1_016 + Mouse_P1_017]
did_it_work <- did_it_work[, .(spacer, t0, tf)]
did_it_work <- melt(did_it_work, id.vars = "spacer", variable.name = "timing", value.name = "CPM")
did_it_work[!spacer %like% "Ctrl", spacer:= "knockdown"]
did_it_work[spacer %like% "Ctrl", spacer:= "control"]
did_it_work <- did_it_work[, .(sumCPM = sum(CPM)), by = .(spacer, timing)]
did_it_work <- dcast(did_it_work, spacer~timing, value.var = "sumCPM")

did_it_work_matrix <- data.matrix(did_it_work[, -1])
rownames(did_it_work_matrix) <- did_it_work[, spacer]
chisq.test(did_it_work_matrix)

did_it_work <- data.table(did_it_work_matrix, keep.rownames = "type")
did_it_work <- melt(did_it_work, variable.name = "timing", value.name = "sum_CPM", id.vars = "type")

did_it_work[
	, sum_CPM := did_it_work[
		i  = type == "control",
		j  = .(ctrl_sum_CPM = sum_CPM),
		by = .(timing)][
			i  = .SD, 
			on = .(timing),
			j  = .(sum_CPM = sum_CPM/ctrl_sum_CPM),
			by = .EACHI]$sum_CPM]

this_plot <- ggplot(data = did_it_work, aes(x = timing , y = sum_CPM, fill = type)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	scale_fill_brewer(palette = "Paired") +
	ggtitle(paste("Barplot. Did it work? Controls, Knockdowns.")) +
	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))

ggthemr("dust")
print(this_plot)
