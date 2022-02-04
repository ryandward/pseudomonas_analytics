require('pacman'); 
p_load(data.table, edgeR, reshape2, digest, forcats, pheatmap)

#######################################################

setup <- fread("/home/ryandward/Mobile_CRISPRi_Neha_Ryan/setup.tsv")

counts <- setup[counts, on = .(condition)]

#######################################################
#######################################################

FL_counts <- FL[counts, on = .(guide), nomatch = 0]

setorder(FL_counts, identifier, treatment, timing, rep, run)

count_matrix <- data.table::dcast(
	FL_counts, 
	forcats::as_factor(promoter) + 
		forcats::as_factor(guide) ~ 
		forcats::as_factor(identifier) + 
		forcats::as_factor(treatment) + 
		forcats::as_factor(timing) + 
		forcats::as_factor(gdna) +
		forcats::as_factor(rep),
	value.var = "count", 
	fill = 0)

groups <- factor(gsub("_[0-9]+$", "", colnames(count_matrix[, -1:-2])),
								 levels = unique(gsub("_[0-9]+$", "", colnames(count_matrix[, -1:-2]))))

permut <- model.matrix( ~ 0 + groups)

colnames(permut) <- unique(groups)

count_data <- data.matrix(count_matrix[, -1:-2])

rownames(count_data) <- paste(count_matrix$promoter, count_matrix$guide)


y <- DGEList(
	counts = count_data,
	group = groups,
	genes = row.names(count_data))

keep <- filterByExpr(y, permut)

y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y)

y <- estimateDisp(y, permut)

plotBCV(y)

fit <- glmQLFit(y, permut, robust = TRUE)

plotQLDisp(fit)


CPM <- cpm(y, prior.count = 0)

plot_matrix <- log10(CPM + 1)

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



#######################################################
#######################################################

contrast <- makeContrasts(
	contrasts = c(
		# "invitro_LB_t1_pellet - invitro_inoculum_t0_pellet",
		# "invitro_LB_t2_pellet - invitro_inoculum_t0_pellet",
		"invitro_gent_t1_pellet - invitro_LB_t1_pellet",
		"invitro_gent_t2_pellet - invitro_LB_t2_pellet"
	), 
	levels = permut)

results <- glmQLFTest(
	fit, 
	contrast = contrast)

results <- topTags(
	results, 
	n = Inf)

results <- data.table(
	results$table)

melted_results <- 
	data.table::melt(
		results, 
		id.vars = c(
			"genes",
			"logCPM",
			"F",
			"FDR",
			"PValue"),
		variable.name = "condition", 
		value.name = "logFC")

melted_results[, condition := gsub("logFC\\.", "", condition)]
melted_results[, condition := gsub("_0", "", condition)]
melted_results[, condition := gsub("\\.\\.\\.", "__", condition)]
melted_results[, c("promoter", "guide") := tstrsplit(genes, " ", type.convert = TRUE, fixed = TRUE)]

melted_results[, genes := NULL]

plot(-log10(FDR) ~ logFC, data = melted_results)


structured_results <- dcast(melted_results, FDR + guide ~ promoter + condition, value.var = "logFC")

plot_matrix <- data.matrix(structured_results[, -1:-2])

rownames(plot_matrix) <- structured_results$guide

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

to_plot_title <- paste("logFC for Pelleted Libraries")

to_plot <- pheatmap(plot_matrix,
										col = plot_colors,
										breaks = breaks,
										border_color = NA,
										cellwidth = 25,
										cellheight = 15,
										# cutree_rows = 4,
										# cutree_cols = 4,
										main = to_plot_title,
										angle_col = 315,
										# fontsize_col = 10,
										fontsize_row = 10,
										cluster_cols = FALSE,
										# cutree_cols = 10,
										clustering_method = "ward.D2",
										clustering_distance_rows = "maximum",
										clustering_distance_cols = "maximum"
)


