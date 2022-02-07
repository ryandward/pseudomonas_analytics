require('pacman'); 
p_load(data.table, edgeR, reshape2, digest, forcats, pheatmap)

clipboard <- function(
	x , 
	sep = "\t" , 
	row.names = FALSE , 
	col.names = TRUE 
) {
	con <- pipe(
		"xclip -selection clipboard -i" , 
		open = "w"
	)
	
	write.table(
		x , 
		con , 
		sep = sep , 
		row.names = row.names , 
		col.names = col.names
	)
	close(con)}

RowVar <- function(x , ...) {
	rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

FL <- fread("awk 'NR%2==1' /home/ryandward/Mobile_CRISPRi_Neha_Ryan/FLGuidesOnly.fasta",
						col.names = c("guide"))

FL[, guide := gsub(">", "", guide)]

promoter_matches <- fread("/home/ryandward/Mobile_CRISPRi_Neha_Ryan/promoter_stats/promoter_matches.tsv",
													col.names = c(
														"ID",
														"promoter",
														"condition"))

promoter_matches[, paste0("ID", 1:2) := tstrsplit(ID, " ", type.convert = TRUE, fixed = TRUE)]
promoter_matches[, ID := NULL]
promoter_matches[, ID2 := NULL]

guide_matches <- fread("/home/ryandward/Mobile_CRISPRi_Neha_Ryan/guide_stats/guide_matches.tsv",
													col.names = c(
														"ID",
														"guide",
														"condition"))

guide_matches[, paste0("ID", 1:2) := tstrsplit(ID, " ", type.convert = TRUE, fixed = TRUE)]
guide_matches[, ID := NULL]
guide_matches[, ID2 := NULL]

setindex(guide_matches, condition, ID1)
setindex(promoter_matches, condition, ID1)


# matches <- promoter_matches[guide_matches, on = .(ID1, condition), nomatch = 0]
matches <- merge(
	promoter_matches, 
	guide_matches, 
	on = .(ID1, condition), 
	all = TRUE)

matches[is.na(promoter), promoter := "unknown"]
matches[is.na(guide), guide := "unknown"]


counts <- matches[
	, 
	.(count = .N),
	by = .(guide, 
				 promoter, 
				 condition)]

rm(promoter_matches)
rm(guide_matches)
rm(matches)

# Simplify names of groups
counts[, paste0("ID", 1:4) := tstrsplit(condition, "_", type.convert = TRUE, fixed = TRUE)]

counts[ID1 == "FL", identifier := "invitro"]
counts[ID1 == "Mouse", identifier := "mouse"]
counts[ID1 == "Undetermined", identifier := "undetermined"]

counts[grep("[0-9]+", ID2), run := ID2]
counts[grep("[0-9]+", ID3), run := ID3]
counts[grep("[0-9]+", ID4), run := ID4]

counts[, run := as.numeric(run)]

# counts[, condition := NULL]
counts[, ID1 := NULL]
counts[, ID2 := NULL]
counts[, ID3 := NULL]
counts[, ID4 := NULL]

############### to next document

# FL_counts <- FL[counts, on = .(guide), nomatch = 0]
# 
# identifier_order <- factor(unique(FL_counts$identifier), levels = unique(FL_counts$identifier))
# 
# # how to widen the data
# count_matrix <- data.table::dcast(FL_counts, promoter + guide ~ identifier +  run, value.var = "count", fill = 0)
# 
# high_levels <- gsub("_.*", "", colnames(count_matrix[,-1:-2]))
# 
# groups <- factor(sort(high_levels), levels = unique(sort(high_levels)))
# 
# permut <- model.matrix( ~ 0 + groups)
# 
# colnames(permut) <- levels(groups)
# 
# count_data <- data.matrix(count_matrix[, -1:-2])
# 
# rownames(count_data) <- paste(count_matrix$promoter, count_matrix$guide)
# 
# 
# y <- DGEList(
# 	counts = count_data,
# 	group = groups,
# 	genes = row.names(count_data))
# 
# keep <- filterByExpr(y, permut)
# 
# y <- y[keep, , keep.lib.sizes = FALSE]
# 
# y <- calcNormFactors(y)
# 
# y <- estimateDisp(y, permut)
# 
# plotBCV(y)
# 
# fit <- glmQLFit(y, permut, robust = TRUE)
# 
# plotQLDisp(fit)
# 
# CPM <- cpm(y, prior.count = 0)
# 
# 
# plot_matrix <- log10(CPM + 1)
# 
# break_halves <- length(unique(as.vector(plot_matrix)))
# 
# breaks <- c(
# 	seq(min(plot_matrix), 
# 			median(plot_matrix), 
# 			length.out = break_halves)[-break_halves], 
# 	seq(median(plot_matrix), 
# 			max(plot_matrix), 
# 			length.out = break_halves))
# 
# plot_colors <- c(
# 	colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves], 
# 	colorRampPalette(c("white", "#007ac1"))(break_halves)[-1])
# 
# to_plot_title <- paste("Log10 Counts per Million")
# 
# to_plot <- pheatmap(plot_matrix,
# 										col = plot_colors,
# 										breaks = breaks,
# 										border_color = NA,
# 										cellwidth = 12,
# 										# cellheight = 12,
# 										# cutree_rows = 4,
# 										# cutree_cols = 4,
# 										main = to_plot_title,
# 										angle_col = 45,
# 										# fontsize_col = 10,
# 										# fontsize_row = 10,
# 										# cluster_cols = FALSE,
# 										# cutree_cols = 10,
# 										clustering_method = "ward.D2",
# 										clustering_distance_rows = "maximum",
# 										clustering_distance_cols = "maximum"
# 										)
# 
# 
# 
# #######################################################
# 
# contrast <- makeContrasts(
# 	contrasts = c(
# 		"mouse - invitro"
# 		), 
# 	levels = permut)
# 
# results <- glmQLFTest(
# 	fit, 
# 	contrast = contrast)
# 
# results <- topTags(
# 	results, 
# 	n = Inf)
# 
# results <- data.table(
# 	results$table)
# 
# melted_results <- 
# 	data.table::melt(
# 		results, 
# 		id.vars = c(
# 			"genes",
# 			"logCPM",
# 			"F",
# 			"FDR",
# 			"PValue"),
# 		variable.name = "condition", 
# 		value.name = "logFC")
# 
# melted_results[, condition := gsub("logFC\\.", "", condition)]
# melted_results[, condition := gsub("_0", "", condition)]
# melted_results[, condition := gsub("\\.\\.\\.", "__", condition)]
# 
# plot(-log10(FDR) ~ logFC, data = melted_results[condition == "logFC"])
# 
