require('pacman')

p_load(vegan, data.table, scales, edgeR, statmod, poolr, pheatmap, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, RColorBrewer)


all_counts <- fread("overall_stats.tsv",
										header = FALSE,
						 col.names = c("spacer",
						 							"count",
						 							"condition"))

exp_design <- fread("exp_design.tsv", 
										 na.strings = c("NA"))

exp_design <- exp_design[condition != "Undetermined"]


################################################################################
# Fix the naming conventions

guide_key <- fread("guide_sequences.tsv", 
									 header = FALSE,
									 col.names = c("name", 
									 							"sequence"))
guide_key[, paste0("ID", 1:4) := tstrsplit(name, "_", type.convert = TRUE, fixed = TRUE)]

NC_002516.2 <- fread("NC_002516.2.bed",
										 header = FALSE,
												 col.names = c(
												 	"chromosome",
												 	"left",
												 	"right",
												 	"locus_tag",
												 	"gene_name",
												 	"strand",
												 	"feature",
												 	"completeness"
												 ))

NC_002516.2 <- NC_002516.2[feature %like% "gene"]
NC_002516.2[, strain := "PAO1"]

annotated_key <- merge(guide_key,
											 header = FALSE,
											 NC_002516.2, 
											 by.x = "ID2", 
											 by.y = "locus_tag", 
											 all.x = FALSE, 
											 all.y = FALSE)

NC_008463.1 <- fread("NC_008463.1.bed",
										 col.names = c(
										 	"chromosome",
										 	"left",
										 	"right",
										 	"locus_tag",
										 	"gene_name",
										 	"strand",
										 	"feature",
										 	"completeness"
										 ))
NC_008463.1 <- NC_008463.1[feature %like% "gene"]
NC_008463.1[, c("strain", "locus_tag") := tstrsplit(locus_tag, "_", type.convert = TRUE, fixed = TRUE)]
NC_008463.1[, strain := "PA14"]

annotated_key <- rbind(annotated_key, 
											 merge(guide_key, 
											 			NC_008463.1, 
											 			by.x = "ID2", 
											 			by.y = "locus_tag", 
											 			all.x = FALSE, 
											 			all.y = FALSE))

lost_guides <- guide_key[!(name %in% annotated_key$name)]

# fix a the "unk" PA3145 that seem to have been messed up... missing a column
lost_guides[ID2 == "unk", `:=` (ID2 = ID1, ID3 = NA, ID4 = ID3)]

# put it back into the annotated guides list
annotated_key <- rbind(annotated_key, 
											 merge(lost_guides, 
											 			NC_002516.2, 
											 			by.x = "ID2", 
											 			by.y = "locus_tag", 
											 			all.x = FALSE, 
											 			all.y = FALSE))

# one gene RS22570 that is still lost
lost_guides <- guide_key[!(name %in% annotated_key$name)]
lost_guides[ID3 == "-", `:=` (ID3 = "hypothetical")]

# add the appropriate columns to just shove it onto the list
lost_guides <- merge(lost_guides, 
			NC_002516.2, 
			by.x = "ID2", 
			by.y = "locus_tag", 
			all.x = TRUE, 
			all.y = FALSE)

# just shove it onto the list
annotated_key <- rbind(annotated_key, lost_guides)

# fix some of the strain names
annotated_key[name %like% "PA14_", strain := "PA14" ]
annotated_key[!name %like% "PA14_", strain := "PAO1" ]
annotated_key[ID2 %like% "RS" & is.na(strain),  `:=` (strain = "PA14", chromosome = "NC_008463.1") ]

if(nrow(annotated_key[ID1 != strain])==0)
{
	annotated_key[, strain := ID1]
}

# guide_key[!(name %in% annotated_key$name)] should now be zero

unfound_guide_message <- paste("The number of guides with information not found:",nrow(guide_key[!(name %in% annotated_key$name)]), ".")
print(paste(unfound_guide_message, "If this number is not 0, something bad happened. Look carefully at the input."))

annotated_key[, type := "unknown"]

annotated_key[gene_name == "." & !(name %like% "Ctrl"), gene_name := ID3]
annotated_key[gene_name == "-", gene_name := "."]

annotated_key[, locus_tag := ID2]
annotated_key[, offset := ID4]

annotated_key[name %like% "Ctrl", `:=` (type = "control", locus_tag = "control", gene_name = "control")]

annotated_key <- annotated_key[, .(chromosome,
																	 left,
																	 right,
																	 locus_tag,
																	 gene_name,
																	 strand,
																	 feature,
																	 completeness,
																	 name,
																	 offset,
																	 type)]
#/fix
################################################################################


rejected_conditions = c("Mouse_P1_020", "Mouse_P1_005", "Mouse_P1_018")
# rejected_conditions = c("Mouse_P1_020")

all_counts <- all_counts[!(condition %in% rejected_conditions)]
exp_design <- exp_design[!(condition %in% rejected_conditions)]

exp_design[media == "Gent", media := "LB"]

################################################################################
# Check for Data Integrity
################################################################################
data_grid <- data.table::dcast(
	all_counts[condition %like% "Mouse"], 
	spacer ~ condition,
	value.var = "count", 
	fill = 0)

data_grid_remelted <- melt(data_grid, variable.name = "condition", value.name = "count", id.vars = c('spacer'))

################################################################################
# Create a correlation matrix for each condition in the remelted matrix.
# It has to be dcasted, then remelted to have zeroes put into their appropriate spots!
#                                        
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
#                                       
# Create a square matrix from the list of pairwise condition correlations.
################################################################################

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

data_permut_check <- melt(data.table(data_permut,keep.rownames = TRUE), id.vars = "rn")[value==1][, .(rn, variable)]
data_permut_check <- data_permut_check[exp_design, on = .(rn==condition)]


data_y <- DGEList(counts = data_grid_matrix, 
									group = data_group, 
									genes = row.names(data_grid_matrix))

data_keep <- filterByExpr(data_y, 
													data_permut)

data_y <- data_y[data_keep, , keep.lib.sizes = FALSE]
data_y <- calcNormFactors(data_y)
data_y <- estimateDisp(data_y, data_permut)

plotBCV(data_y)

data_fit <- glmQLFit(data_y, data_permut, robust = TRUE)

plotQLDisp(data_fit)

data_CPM <- cpm(data_y, prior.count = 0)

colnames(data_CPM) <- factor(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")])

################################################################################
# MDS PLOT

plotMDS(data_y)

# 
# colors <- rep(
# 	c(alpha("black", 0.25),
# 		alpha("red", 0.5)), 2)
# 
# pch <- c(
# 	rep(19, 2), 
# 	rep(21, 2))
# 
# plotMDS(aba_y, 
# 				col = colors[aba_group], 
# 				pch = pch[aba_group], 
# 				cex = 3)
# 
# legend("bottomleft", 
# 			 legend = levels(aba_group), 
# 			 pch = pch, 
# 			 col = colors, 
# 			 ncol = 2)
# 
# title("MDS: Guide Count")
# /MDS
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
			length.out = break_halves))

breaks <- breaks[-length(breaks)]
breaks <- c(breaks, 0.99999999)

plot_colors <- c(
	colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves], 
	colorRampPalette(c("white", "#007ac1"))(break_halves)[-c(1, break_halves)])

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
# / CPM heatmap
################################################################################
# 
# contrast_levels <- c("Gent_pellet_6_generations - LB_pellet_6_generations",
# 										 "Gent_plated_6_generations - LB_plated_6_generations",
# 										 "LB_plated_6_generations - LB_pellet_6_generations",
# 										 "Gent_plated_6_generations - Gent_pellet_6_generations",
# 										 "inoculum_plated_t0 - inoculum_pellet_t0",
# 										 "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0",
# 										 "mouse_plated_100x_inoculum_dilution - inoculum_plated_t0",
# 										 "mouse_pellet_10x_inoculum_dilution - inoculum_pellet_t0",
# 										 "inoculum_pellet_t0 - inoculum_plated_t0")

# There doesn't seem to be much difference when adding gent to the plates, so I modified
# the design matrix ~ line 20.

contrast_levels <- c("LB_plated_6_generations - LB_pellet_6_generations",
										 "inoculum_plated_t0 - inoculum_pellet_t0",
										 "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0",
										 "mouse_plated_100x_inoculum_dilution - inoculum_plated_t0",
										 "mouse_pellet_10x_inoculum_dilution - inoculum_pellet_t0",
										 "inoculum_pellet_t0 - inoculum_plated_t0")

# contrast_levels <- c("inoculum_pellet_t0 - inoculum_plated_t0")

contrast_levels <- c("mouse_plated_10x_inoculum_dilution - inoculum_plated_t0",
										 "mouse_plated_10x_inoculum_dilution - LB_plated_6_generations",
										 "mouse_pellet_10x_inoculum_dilution - inoculum_pellet_t0",
										 "mouse_pellet_10x_inoculum_dilution - LB_pellet_6_generations")


data_contrast <- makeContrasts(contrasts = contrast_levels, 
															levels = data_permut)

################################################################################

results_FDR <- all_counts[, .(genes = unique(spacer))]
results_LFC <- all_counts[, .(genes = unique(spacer))]

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

melted_results_LFC[
	, LFC := melted_results_LFC[
		i  = genes %like% "Ctrl",
		j  = .(ctrl_medLFC = median(LFC, na.rm = TRUE)),
		by = .(condition)][
			i  = .SD, # In this case, .SD is multiple in nature -- it refers to each of these sub-data.tables, one-at-a-time
			on = .(condition),
			j  = .(adj_medLFC = LFC - ctrl_medLFC),
			by = .EACHI]$adj_medLFC]

melted_results <- 
	melted_results_LFC[
		melted_results_FDR, 
		on = .(genes, condition)]

################################################################################

melted_results <- annotated_key[melted_results, on = .(name == genes)]

################################################################################

median_melted_results <- melted_results[, .(
	medLFC = median(LFC), 
	FDR = stouffer(FDR)$p), 
	by = .(locus_tag, gene_name, type, condition)]

################################################################################

setorder(median_melted_results, FDR)

################################################################################

plot(-log10(FDR)~LFC, data = melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0"], pch = 20)
points(-log10(FDR)~LFC, data = melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0" & type == "control"], pch = 20, col = "green")

control_CPM <- melt(data.table(data_CPM, keep.rownames = TRUE), id.vars = "rn")[rn %like% "Ctrl"]
control_CPM <- control_CPM[variable == "mouse_plated_10x_inoculum_dilution"]
setorder(control_CPM, rn)

all_counts[spacer %like% "Ctrl", .(median_ctrl = median(count)), by = .(condition)][data_permut_check, on = .(condition==rn)][, .(condition, rep, variable, median_ctrl)]

################################################################################
################################################################################
CPM_reps <- data.table(data_CPM[, c(8,9,10)], keep.rownames = "spacer")
setnames(CPM_reps, c("spacer", "rep_1", "rep_2", "rep_3"))
CPM_reps <- melt(CPM_reps, variable.name = "rep", value.name = "CPM", id.vars = "spacer")
CPM_reps[spacer %like% "Ctrl", .(median_ctrl = median(CPM)), by = .(rep)]
################################################################################
################################################################################
CPM_reps <- data.table(data_CPM[, c(8,9,10)], keep.rownames = "spacer")
setnames(CPM_reps, c("spacer", "rep_1", "rep_2", "rep_3"))
CPM_reps <- melt(CPM_reps, variable.name = "rep", value.name = "CPM", id.vars = "spacer")
plot(hist(CPM_reps[spacer %like% "Ctrl", .(std_CPM = std(CPM)), by = .(spacer)]$std_CPM), main = "Standard Deviation of Controls")
################################################################################
################################################################################

mouse_grid_FDR <- dcast(median_melted_results, locus_tag + gene_name + type ~ condition, value.var = "FDR")
mouse_grid_LFC <- dcast(median_melted_results, locus_tag + gene_name + type ~ condition, value.var = "medLFC")

plot(-log10(FDR)~medLFC, data = median_melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0"], pch = 20, main = "Median: mouse_plated_10x_inoculum_dilution - inoculum_plated_t0")
abline(v = median_melted_results[condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0" & type == "control"]$medLFC)

FDR_order <- median_melted_results[, .(mean_order = mean(FDR)), by = .(locus_tag, gene_name, type)]
setorder(FDR_order, mean_order, locus_tag)
FDR_order[, FDR := .I]
FDR_order[, medLFC := .I]
FDR_order <- FDR_order[, .(locus_tag, gene_name, type, FDR, medLFC)]
FDR_order[, condition := "order"]

# glue the dummy variable onto the results matrix
median_melted_results <- rbind(FDR_order, median_melted_results)

################################################################################
# construct the results then copy them to clipboard

# construct the results, with the dummy variable, under the condition = "order"
to_clip <- dcast(median_melted_results, locus_tag + gene_name + type ~ condition, value.var = "FDR")
# order based on the "order" condition
setorder(to_clip, order)
# remove the "order" condition
to_clip <- to_clip[, -c("order")]
# copy it to clipboard
clipr::write_clip(to_clip)

## Repeat the same thing for medLFC
to_clip <- dcast(median_melted_results, locus_tag + gene_name + type ~ condition, value.var = "medLFC")
setorder(to_clip, order)
to_clip <- to_clip[, -c("order")]
clipr::write_clip(to_clip)
#/construct
################################################################################
