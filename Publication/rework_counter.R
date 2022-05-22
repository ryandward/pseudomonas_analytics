# counting guidelines: fastq_operations/link_promoters_guides_then_count.sh 
# GLBRC work directory: /home/GLBRCORG/ryan.d.ward/MouseLib
# GLBRC work directory, replicates: /home/GLBRCORG/ryan.d.ward/experiment_35362

require('pacman')

p_load(
	data.table,
	tidyverse,
	scales,
	edgeR,
	pheatmap,
	svglite,
	ggplot2,
	ggrepel,
	colourpicker,
	RColorBrewer,
	poolr,
	statmod,
	viridis,
	hrbrthemes
)

annotated_key <- fread("Publication/annotated_key.tsv")

focused <- 
	fread(
		"Publication/targeted_library.txt",
		col.names = c("name"))

exp_design <- fread("Publication/working_experimental_design.tsv")

setorder(exp_design, condition)

##########################################################################################

all_counts <- 
	fread(
		"Publication/all_counts_seal.tsv.gz", 
		header = FALSE, 
		col.names = c(
			"promoter", 
			"spacer", 
			"count", 
			"condition"))

inoculum_counts <- 
	fread(
		"Publication/inoculum_counts.tsv.gz", 
		header = FALSE, 
		col.names = c(
			"promoter", 
			"spacer", 
			"count"))

all_counts2 <- fread(
	"Publication/all_counts_seal2.tsv.gz",
	header = FALSE,
	col.names = c(
		"count",
		"condition",
		"spacer"))

# For all the results we have from the sequencing run performed at UWM, 
# the data are presumed to be P1.

all_counts2 <- all_counts2 %>%
	mutate(promoter = "P1")

inoculum_counts[, condition := "Inoculum"]

all_counts <- rbind(all_counts, inoculum_counts)

all_counts <- rbind(all_counts, all_counts2)

mating_counts <- fread(
	"additional_sequences/additional.tsv.gz",
	header = FALSE, 
	col.names = c(
		"count",
		"condition",
		"promoter",
		"spacer")) %>% 
	filter(condition == "P1_mfdpir")

all_counts <- rbind(mating_counts, all_counts)


# all_counts <- all_counts[promoter == "P1"][, .(spacer, condition, count)]
# exp_design <- exp_design %>% filter(condition_group %in% c("inoculum_pellet_t0", "mouse_plated_10x_inoculum_dilution", "LB_plated_6_generations", "P1_mating_strain"))

exp_design <- exp_design %>% filter(!growth_condition %like% "CLEAR") 

all_counts <-
	all_counts[condition %in% exp_design$condition] [, .(spacer, condition, count)]


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

################################################################################

data_grid_matrix <- data.matrix(data_grid[, -c("spacer")])

row.names(data_grid_matrix) <- data_grid$spacer

data_grid_matrix %>% cpm %>% 
	as.data.table(keep.rownames = "spacer") %>% 
	melt(id.vars = "spacer", variable.name = "condition", value.name = "CPM") %>% 
	inner_join(annotated_key, by = c("spacer" = "name")) %>% select(spacer, condition, type, CPM) %>% 
	inner_join(exp_design) %>%
	mutate(rep = factor(as.character(rep))) %>%
	ggplot(aes(x = CPM, fill = rep)) + 
	geom_density() + 
	scale_x_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^(1:6)),
		labels = label_number_si()) + 
	scale_fill_viridis(alpha = 0.35, discrete = T, direction = -1) +
	theme_ipsum_rc() +
	facet_wrap(facets = c("condition_group", "type"), ncol = 6) -> p

print(p)

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
		id.vars = "condition",
		variable.name = "condition_group")[value == 1][, .(condition, condition_group)]

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

contrast_levels <-
	c("mouse_plated_10x_inoculum_dilution - LB_plated_6_generations",
		"mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0")

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

melted_results_LFC[, LFC := round(LFC, 3)]

################################################################################

melted_results <-
	melted_results_LFC[
		melted_results_FDR,
		on = .(genes, condition)]

melted_results <- melted_results[!is.na(FDR)]

melted_results_by_condition <- 
	melted_results[
		genes %like% "Ctrl",
		.(med_LFC = median(LFC)),
		keyby = .(condition)]

setkey(melted_results, condition)

melted_results[, LFC := melted_results_by_condition[
	melted_results, LFC - med_LFC, by = .EACHI]$V1]


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
median_melted_results[gene_name == "", gene_name_stylized := paste0("bold('", locus_tag, "')")]
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

# grouped_CPM <- copy(CPM_melted)

grouped_CPM <- CPM_melted[data_permut_check[, .(condition, condition_group)], on = . (condition)]

grouped_CPM[
	condition_group == "inoculum_pellet_t0", Condition := 'Inoculum (t0)']

grouped_CPM[
	condition_group == "mouse_plated_10x_inoculum_dilution",	Condition := 'Lung (10x dilution, plate)']

grouped_CPM[
	condition_group == "mouse_plated_100x_inoculum_dilution",	Condition := 'Lung (100x dilution, plate)']

grouped_CPM[
	condition_group == "mouse_pellet_10x_inoculum_dilution",	Condition := 'Lung (10x dilution, pellet)']

grouped_CPM[
	condition_group == "inoculum_plated_t0",	Condition := 'Inoculum grown on plates']

grouped_CPM[is.na(Condition), Condition := condition_group]

grouped_CPM <- exp_design[grouped_CPM, on = .(condition)]

grouped_CPM <- annotated_key[grouped_CPM, on = .(name == spacer)]

grouped_CPM[!is.na(rep), verbose := paste(media, gDNA_source, growth_condition, rep, sep = "_")]

grouped_CPM[is.na(rep), verbose := paste(media, gDNA_source, growth_condition, sep = "_")]

setorder(median_melted_results, locus_tag)

results_summary <- melted_results[FDR < 0.05, .N, by = .(condition)]
median_results_summary <- median_melted_results[FDR < 0.05, .N, by = .(condition)]

##########################################################################################
# function to create boxplots from CPM


box_CPM <- function(this_gene) {
	print(
		grouped_CPM %>% filter(
			condition %in% c(
				"Mouse_P1_015", 
				"Mouse_P1_016", 
				"Mouse_P1_017", 
				"Inoculum", 
				"Mouse_P1_003",
				"Mouse_P1_004",
				"dJMP4"))  %>%
			filter(gene_name == this_gene | locus_tag == this_gene) %>%
			mutate(Guide = factor(
				as.character(offset),
				levels = as.character(sort(unique(offset))))) %>%
			ggplot(aes(x = Condition, y = CPM, fill = Guide, color = Guide)) +
			geom_boxplot(outlier.colour = NA, alpha = 0.5) +
			geom_point(position = position_jitterdodge()) +
			ylab("Counts per Million") +
			xlab("Condition") +
			theme_ipsum() +
			scale_fill_ipsum() +
			scale_color_ipsum() +
			scale_y_continuous(
				trans = scales::pseudo_log_trans(base = 10),
				breaks = c(0, 10^(1:5)),
				labels = label_number_si(),
				limits = c(0, cpm.max)) +
			ggtitle(bquote(bold(Guides ~ Recovered ~ `for` ~ bolditalic(.(this_gene)) ~ "(CPM)"))))
}

