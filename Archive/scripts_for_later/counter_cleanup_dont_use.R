require('pacman')

p_load(data.table, 
			 ggthemr, 
			 edgeR, 
			 pheatmap, 
			 ggplot2, ggrepel,
			 RColorBrewer)

##### READ THINGS FROM DISK #####
annotated_key <- fread("annotated_key.tsv")

exp_design <- fread("exp_design.tsv")

all_counts <- fread(
	"all_counts_seal.tsv", 
	header = FALSE, 
	col.names = c("promoter", "spacer", "count", "condition"))

##### ONLY KEEP PROMOTER 1 #####
all_counts <- all_counts[promoter == "P1"][, .(spacer, condition, count)]

##### ORDER THINGS IN THE SAME ORDER #####
setorder(exp_design, condition)
setorder(all_counts, condition)

##### WHAT ARE WE INTERESTED IN? #####
main_question <- c(
	"Mouse_P1_003",
	"Mouse_P1_006",
	"Mouse_P1_015", 
	"Mouse_P1_016", 
	"Mouse_P1_017")

##### FILTER OUR TABLES TO OUR INTERESTING QUESTION #####
exp_design <- exp_design[condition %in% main_question]
all_counts <- all_counts[condition %in% main_question]

##### WE NEED ZEROES WHERE DATA IS MISSING #####
##### MAKE LIST INTO A GRID #####
data_grid <- data.table::dcast(
	all_counts, 
	spacer ~ condition,
	value.var = "count",
	fill = 0)

##### TURN GRID BACK INTO LIST #####
all_counts <- melt(data_grid, variable.name = "condition", value.name = "count", id.vars = c('spacer'))

##### EDGER USES MATRIXES AND NOT DATA.TABLES #####
data_grid_matrix <- data.matrix(data_grid[, -c("spacer")])

##### MATRICES USE ROWNAMES> #####
row.names(data_grid_matrix) <- data_grid$spacer

#### EXPERIMENTAL GROUPS WITH LEVELS #####
data_group <- factor(
	exp_design[, paste(media, gDNA_source, growth_condition, sep = "_")])

##### CREATE A MATRIX AS DEFINED BY THE EXPERIMENTAL GROUPS WITH LEVELS #####
data_permut <- model.matrix(~ 0 + data_group)

##### GIVE THE APPROPRIATE ROW AND COLUMN NAMES #####
colnames(data_permut) <- levels(data_group)
rownames(data_permut) <- colnames(data_grid_matrix)

##### CREATE A CHECK DATA.TABLE TO COMPARE IF ALL THIS WORKED ######
data_permut_check <- melt(
	data.table(
		data_permut,
		keep.rownames = "condition"),
	id.vars = "condition",
	variable.name = "group" )[value == 1][, .(condition, group)]

data_permut_check <- data_permut_check[exp_design, on = .(condition == condition)]

##### CHECK THE TABLE VERY CLOSELY FOR IRREGULARITIES BETWEEN COLUMNS #####
print(data_permut_check)

##### CREATE EDGER OBJECT CONTAINING COUNTS, GROUPS, AND GENES
data_y <- DGEList(
	counts = data_grid_matrix,
	group = data_group,
	genes = row.names(data_grid_matrix))

##### KEEP ONLY THINGS WITH 10+ COUNTS IN 2+ CONDITIONS
data_keep <- rowSums(cpm(data_y) > 10) >= 2
data_y <- data_y[data_keep, , keep.lib.sizes = FALSE]

##### NORMALIZE THE COUNT DATA #####
data_y <- calcNormFactors(data_y)

##### CALCULATE DISPERSION #####
data_y <- estimateDisp(data_y, data_permut)

##### FIT A LINE THROUGH THE DISPERSION #####
data_fit <- glmQLFit(data_y, data_permut, robust = TRUE)

##### CALCULATE COUNTS PER MILLION #####
data_CPM <- cpm(data_y, prior.count = 0)

##### WHAT CONTRASTS BETWEEN GROUPS ARE WE INTERESTED IN? #####
contrast_levels <- c("mouse_plated_10x_inoculum_dilution - inoculum_plated_t0",
										 "inoculum_plated_t0 - inoculum_pellet_t0")

##### MAKE THE CONTRAST FROM THE MATRICES #####
data_contrast <- makeContrasts(
	contrasts = contrast_levels,
	levels = data_permut)


##### EMPTY FDR AND LFC MATRICES TO POPULATE #####
results_FDR <- all_counts[, .(genes = unique(spacer))]
results_LFC <- all_counts[, .(genes = unique(spacer))]

##### POPULATE THE MATRICES BY ITERATING THROUGH THE GROUP CONTRASTS #####
for (i in 1:ncol(data_contrast)){

	results <- glmQLFTest(data_fit, contrast = data_contrast[,i])
	results <- topTags(results, n = Inf)
	results <- data.table(results$table)

	print(paste("Processing results for", contrast_levels[i], "..."))

	results_FDR <- results[, .(genes, FDR)][results_FDR, on = .(genes)]
	setnames(results_FDR, "FDR", contrast_levels[i])

	results_LFC <- results[, .(genes, logFC)][results_LFC, on = .(genes)]
	setnames(results_LFC, "logFC", contrast_levels[i])}


##### CLEAN UP THE FDR, LFC, THEN MERGE THEM #####
melted_results_FDR <-
	data.table::melt(
		results_FDR,
		id.vars = c(
			"genes"),
		variable.name = "condition",
		value.name = "FDR",
		measure.vars = contrast_levels)

melted_results_FDR <- melted_results_FDR[!is.na(FDR)]

melted_results_LFC <-
	data.table::melt(
		results_LFC,
		id.vars = c(
			"genes"),
		variable.name = "condition",
		value.name = "LFC",
		measure.vars = contrast_levels)

melted_results <-
	melted_results_LFC[
		melted_results_FDR,
		on = .(genes, condition)]
# 
# melted_results <- melted_results[!is.na(FDR)]
# 
# melted_results[
# 	, LFC := melted_results[
# 		i  = genes %like% "Ctrl",
# 		j  = .(ctrl_medLFC = median(LFC, na.rm = TRUE)),
# 		by = .(condition)][
# 			i  = .SD, 
# 			on = .(condition),
# 			j  = .(adj_medLFC = LFC - ctrl_medLFC),
# 			by = .EACHI]$adj_medLFC]
# 
# ################################################################################
# 
# melted_results <- annotated_key[melted_results, on = .(name == genes)]
# 
# ################################################################################
# 
# 
# guide_toss <- melted_results[, .(sumFDR = sum(FDR)), by = .(name, gene_name, locus_tag)]
# 
# setorder(guide_toss, sumFDR)
# 
# guide_rejected <- guide_toss[, .(maxsumFDR = max(sumFDR)), by = .(gene_name, locus_tag)]
# 
# guide_rejected <- guide_toss[guide_rejected, on = .(locus_tag, gene_name)][sumFDR == maxsumFDR, .(name)]
# 
# melted_results_adjusted <- melted_results[!(name %in% guide_rejected$name)]
# 
# ################################################################################
# 
# median_melted_results <- melted_results_adjusted[, .(
# 	medLFC = median(LFC), 
# 	FDR = stouffer(FDR)$p), 
# 	by = .(locus_tag, gene_name, type, condition)]
# 
# ################################################################################
# 
# setorder(median_melted_results, FDR)
# 
# ################################################################################
# 
# plot(-log10(FDR)~LFC, data = melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0"], pch = 20, main = "Guide-level: mouse_plated_10x_inoculum_dilution - inoculum_plated_t0")
# points(-log10(FDR)~LFC, data = melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0" & type == "control"], pch = 20, col = "green")
# 
# plot(-log10(FDR)~medLFC, data = median_melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0"], pch = 20, main = "Gene-level (Median guide): mouse_plated_10x_inoculum_dilution - inoculum_plated_t0")
# abline(v = median_melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0" & type == "control"]$medLFC)
# 
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# # control_CPM <- melt(data.table(data_CPM, keep.rownames = TRUE), id.vars = "rn")[rn %like% "Ctrl"]
# # control_CPM <- control_CPM[variable == "mouse_plated_10x_inoculum_dilution"]
# # setorder(control_CPM, rn)
# # 
# # all_counts[spacer %like% "Ctrl", .(median_ctrl = median(count)), by = .(condition)][data_permut_check, on = .(condition==rn)][, .(condition, rep, variable, median_ctrl)]
# # 
# # ################################################################################
# # ################################################################################
# # 
# # mouse_grid_FDR <- dcast(median_melted_results, locus_tag + gene_name + type ~ condition, value.var = "FDR")
# # mouse_grid_LFC <- dcast(median_melted_results, locus_tag + gene_name + type ~ condition, value.var = "medLFC")
# # 
# # ################################################################################
# # ################################################################################
# # 
# # FDR_order <- median_melted_results[, .(mean_order = mean(FDR)), by = .(locus_tag, gene_name, type)]
# # setorder(FDR_order, mean_order, locus_tag)
# # FDR_order[, FDR := .I]
# # FDR_order[, medLFC := .I]
# # FDR_order <- FDR_order[, .(locus_tag, gene_name, type, FDR, medLFC)]
# # FDR_order[, condition := "order"]
# # 
# # # glue the dummy variable onto the results matrix
# # # median_melted_results <- rbind(FDR_order, median_melted_results)
# # 
# # ################################################################################
# # # construct the results then copy them to clipboard
# # 
# # to_clip <- dcast(rbind(FDR_order, median_melted_results), locus_tag + gene_name + type ~ condition, value.var = "FDR")
# # setorder(to_clip, order)
# # to_clip <- to_clip[, -c("order")]
# # clipr::write_clip(to_clip)
# # 
# # to_clip <- dcast(rbind(FDR_order, median_melted_results), locus_tag + gene_name + type ~ condition, value.var = "medLFC")
# # setorder(to_clip, order)
# # to_clip <- to_clip[, -c("order")]
# # clipr::write_clip(to_clip)
# # 
# # plot(-log10(FDR) ~ medLFC, data = median_melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0"])
# 
# ##################################
# # Chi Square Test, did it work?
# ##################################
# 
# did_it_work <- data.table(data_CPM, keep.rownames = "spacer")[, .(spacer, Mouse_P1_003, Mouse_P1_015, Mouse_P1_016, Mouse_P1_017)]
# 
# did_it_work[, t0 := Mouse_P1_003]
# did_it_work[, tf := Mouse_P1_015 + Mouse_P1_016 + Mouse_P1_017]
# did_it_work <- did_it_work[, .(spacer, t0, tf)]
# did_it_work <- melt(did_it_work, id.vars = "spacer", variable.name = "timing", value.name = "CPM")
# did_it_work[!spacer %like% "Ctrl", spacer:= "knockdown"]
# did_it_work[spacer %like% "Ctrl", spacer:= "control"]
# did_it_work <- did_it_work[, .(sumCPM = sum(CPM)), by = .(spacer, timing)]
# did_it_work <- dcast(did_it_work, spacer~timing, value.var = "sumCPM")
# 
# did_it_work_matrix <- data.matrix(did_it_work[, -1])
# rownames(did_it_work_matrix) <- did_it_work[, spacer]
# chisq.test(did_it_work_matrix)
# 
# did_it_work <- data.table(did_it_work_matrix, keep.rownames = "type")
# did_it_work <- melt(did_it_work, variable.name = "timing", value.name = "sum_CPM", id.vars = "type")
# 
# did_it_work[
# 	, sum_CPM := did_it_work[
# 		i  = type == "control",
# 		j  = .(ctrl_sum_CPM = sum_CPM),
# 		by = .(timing)][
# 			i  = .SD, 
# 			on = .(timing),
# 			j  = .(sum_CPM = sum_CPM/ctrl_sum_CPM),
# 			by = .EACHI]$sum_CPM]
# 
# this_plot <- ggplot(data = did_it_work, aes(x = timing , y = sum_CPM, fill = type)) +
# 	geom_bar(stat = "identity", position = position_dodge()) +
# 	scale_fill_brewer(palette = "Paired") +
# 	ggtitle(paste("Barplot. Did it work? Controls, Knockdowns.")) +
# 	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))
# 
# ggthemr("dust")
# print(this_plot)

