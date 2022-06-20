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
	statmod
)

##########################################################################################
# girl, fix your annotations

annotations <- fread("unified_counts/full_annotations.tsv")

annotations <- annotations %>% 
	mutate(locus_tag = case_when(
		gene == "Ctrl" ~ gene, !is.na(gene) & gene != "" & locus_tag != "" & locus_tag != "unknown" ~ locus_tag, 
		locus_tag == "" ~ gsub("PA14_", "", PAO1_locus_tag), 
		PA14_locus_tag != "" & locus_tag == "unknown" ~ PA14_locus_tag, 
		sequence == "CTGTGGGAGCTACGGGCATT" ~ "PA3145"), 
		gene = case_when(gene == "unnamed" | gene == "" ~ locus_tag, locus_tag == "PA3145" ~ "wbpL", 
										 TRUE ~ gene))

non_unique_genes <- annotations %>% select(locus_tag, gene) %>% unique %>% select(gene) %>% group_by(gene) %>% summarise(count = 1:n()) %>% filter(count > 1) %>% pull(gene)

annotations <- 
	annotations %>% 
	mutate(gene = case_when(
		gene %in% non_unique_genes ~ paste0(gene, "_" ,locus_tag), 
		TRUE ~ gene))

non_unique_locus_tags <- annotations %>% select(locus_tag, gene) %>% unique %>% select(locus_tag) %>% group_by(locus_tag) %>% summarise(count = 1:n()) %>% filter(count > 1) %>% pull(locus_tag)

annotations <- 
	annotations %>% 
	mutate(gene = case_when(
		locus_tag == "PA3145" ~ "wbpL",
		locus_tag == "RS09495" ~ "orfN",
		TRUE ~ gene))

##########################################################################################


exp_design <- fread("unified_counts/exp_design.tsv")

exp_design <- exp_design %>% group_by(media, gDNA_source, growth_condition) %>% mutate(rep = 1:n()) %>% data.table

exp_design[, sample_name := paste(media, gDNA_source, growth_condition, rep, sep = "_")]

exp_design[, sample_group := paste(media, gDNA_source, growth_condition, sep = "_")]


exp_design <- exp_design[!condition %in% c("dJMP3", "dJMP5")]

setorder(exp_design, condition)

##########################################################################################

all_counts1 <-
	fread("unified_counts/all_counts_seal1.tsv.gz", col.names = c("condition", "sequence", "promoter", "count"),
				na.strings = "NA")

all_counts2 <-
	fread("unified_counts/all_counts_seal2.tsv.gz", col.names = c("condition", "sequence", "promoter", "count"),
				na.strings = "NA")

all_counts <- rbind(all_counts1, all_counts2)

all_counts <- all_counts %>% filter(condition %like% "P1" | condition %like% "dJMP")

################################################################################
# raw count density plots
all_counts %>% 
	ggplot(aes(x = count)) + 
	geom_density(aes(fill = promoter), alpha = 0.25) + 	
	scale_x_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^(1:6)),
		labels = label_number_si()) + 
	facet_wrap(facets = "condition") +
	ggtitle("Raw Density Counts") -> p

print(p)
################################################################################
# let's pretend all the promoters are P1
all_counts %>% group_by(condition, sequence) %>% 
	summarise(count = sum(count)) %>% 
	ggplot(aes(x = count)) + 
	geom_density(aes(), alpha = 0.25) + 	
	scale_x_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^(1:6)),
		labels = label_number_si()) + 
	facet_wrap(facets = "condition") +
	ggtitle("Let's ignore promoters") -> p

print(p)
################################################################################
# let's pretend all the promoters are P1 & ignore stuff that wasn't in the mating strain
all_counts %>% 
	filter(sequence %in% (all_counts %>% filter(condition == "P1_mfdpir") %>% pull(sequence))) %>%
	group_by(condition, sequence) %>%
	summarise(count = sum(count)) %>%
	ggplot(aes(x = count)) + 
	geom_density(aes(), alpha = 0.25) + 	
	scale_x_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^(1:6)),
		labels = label_number_si()) + 
	facet_wrap(facets = "condition") +
	ggtitle("Let's ignore promoters, and extra guides found") -> p

print(p)
################################################################################

all_counts <- all_counts %>% 
	filter(sequence %in% (all_counts %>% filter(condition == "P1_mfdpir") %>% pull(sequence))) %>%
	group_by(condition, sequence) %>%
	summarise(count = sum(count))

#add zeroes and filter one last time to make sure only stuff in mating strain is present
all_counts <- all_counts %>% 
	pivot_wider(id_cols = sequence, names_from = condition, values_from = count, values_fill = 0) %>% 
	pivot_longer(!sequence, names_to = "condition", values_to = "count") %>%
	filter(sequence %in% (all_counts %>% filter(condition == "P1_mfdpir") %>% pull(sequence)))

#add a CPM value

all_counts <- 
	all_counts %>% group_by(condition) %>% mutate(CPM = 1e6 * count/sum(count))

################################################################################

all_counts %>%
ggplot(aes(x = CPM)) + 
	geom_density(aes(), alpha = 0.25) + 	
	scale_x_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^(1:6)),
		labels = label_number_si()) + 
	facet_wrap(facets = "condition") -> p

print(p)

################################################################################

all_counts %>% 
	inner_join(exp_design) %>%
	ggplot(aes(x = CPM)) + 
	geom_density(aes(fill = as.character(rep)), alpha = 0.25) + 	
	scale_x_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^(1:6)),
		labels = label_number_si()) + 
	facet_wrap(facets = "sample_group") -> p

print(p)

################################################################################

# Figure 2
all_counts %>% 
	inner_join(exp_design) %>%
	inner_join(annotations) %>% 
	mutate(type = case_when(gene %like% "Ctrl" ~ "control", TRUE ~ "knockdown")) %>% 
	filter(sample_group %like% "mating" | sample_group %like% "^inoculum_pellet") %>%
	mutate(sample_group = case_when(
		sample_group %like% "mating" ~ "Mating Strain",
		sample_group %like% "^inoculum_pellet" ~ "PA14 Library")) %>%
	arrange(sample_group) %>%
	mutate(sample_group = factor(sample_group, levels = unique(sample_group))) %>%
	ggplot(aes(x = CPM)) + 
	geom_density(aes(fill = as.character(rep)), alpha = 0.25) + 	
	scale_x_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^(1:6)),
		labels = label_number_si()) + 
	facet_grid(facets = c("type", "sample_group")) -> p

print(p)



################################################################################

setorder(all_counts, condition)

setorder(exp_design, condition)

all_counts <- all_counts %>% 
	data.table %>%
	`[`(condition %in% exp_design$condition)

all_counts <- all_counts %>% inner_join(annotations)

################################################################################

data_grid <- all_counts %>% pivot_wider(id_cols = sequence, names_from = condition, values_from = count)

data_grid_matrix <- data_grid %>% select(-sequence) %>% data.matrix

row.names(data_grid_matrix) <- data_grid$sequence

colnames(data_grid_matrix) <- data_grid %>% colnames %>% data.table(condition = .) %>% inner_join(exp_design) %>% pull(sample_name)

crossjoin_correlation_grid <- cor(data_grid_matrix)

# Create a square matrix from the list of pairwise condition correlations.

################################################################################

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
	clustering_distance_rows = "canberra",
	clustering_distance_cols = "canberra")

print(to_plot)

################################################################################

data_group <-
	factor(
		exp_design[, sample_group],
		levels = unique(exp_design$sample_group))

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
# 
# plot_grid <- cor(data_CPM)
# 
# plot_matrix <- data.matrix(plot_grid)
# 
# break_halves <- length(unique(as.vector(plot_matrix)))
# 
# breaks <- c(
# 	seq(min(plot_matrix),
# 			median(plot_matrix),
# 			length.out = break_halves)[-break_halves],
# 	seq(median(plot_matrix),
# 			max(plot_matrix),
# 			length.out = break_halves)
# )
# 
# breaks <- breaks[-length(breaks)]
# breaks <- c(breaks, 0.99999999)
# 
# plot_colors <-
# 	colorRampPalette(c("white", "#007ac1"))(break_halves * 2 - 1)[-c(1, break_halves)]
# 
# plot_colors <- c(plot_colors, "dark grey")
# 
# to_plot_title <-
# 	paste("Clustering of Conditional Correlations (CPM)")
# 
# to_plot <- pheatmap(
# 	plot_matrix,
# 	col = plot_colors,
# 	breaks = breaks,
# 	border_color = NA,
# 	cellwidth = 20,
# 	cellheight = 20,
# 	main = to_plot_title,
# 	angle_col = 315,
# 	# fontsize_col = 10,
# 	# fontsize_row = 10,
# 	# cluster_cols = FALSE,
# 	show_rownames = TRUE,
# 	show_colnames = TRUE,
# 	clustering_method = "ward.D2",
# 	clustering_distance_rows = "maximum",
# 	clustering_distance_cols = "maximum"
# )
# print(to_plot)

################################################################################
# label with what the sequencing experiment actually is
# 
# to_plot <- pheatmap(
# 	plot_matrix,
# 	col = plot_colors,
# 	breaks = breaks,
# 	border_color = NA,
# 	cellwidth = 20,
# 	cellheight = 20,
# 	main = to_plot_title,
# 	angle_col = 315,
# 	# fontsize_col = 10,
# 	# fontsize_row = 10,
# 	# cluster_cols = FALSE,
# 	show_rownames = TRUE,
# 	show_colnames = TRUE,
# 	clustering_method = "ward.D2",
# 	labels_row = factor(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")]),
# 	labels_col = factor(exp_design[,  paste(media, gDNA_source, growth_condition, sep = "_")]),
# 	clustering_distance_rows = "maximum",
# 	clustering_distance_cols = "maximum"
# )
# print(to_plot)

################################################################################

# data_CPM_by_group <- 
# 	copy(data_CPM)
# 
# colnames(data_CPM_by_group) <- 
# 	factor(exp_design[,  paste(
# 		media, gDNA_source, growth_condition, sep = "_")])
#
# print(to_plot)

################################################################################
# other diagnostic plots

# plotMDS(log2(data_CPM_by_group))
# plotQLDisp(data_fit)
# plotBCV(data_y)

################################################################################

contrast_levels <-
	c("LB_plated_6_generations - inoculum_pellet_t0",
		"mouse_plated_10x_inoculum_dilution - LB_plated_6_generations",
		"mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0")

data_contrast <- makeContrasts(
	contrasts = contrast_levels,
	levels = data_permut)

################################################################################

results_FDR <- all_counts[, .(genes = unique(sequence))]

results_LFC <- all_counts[, .(genes = unique(sequence))]

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
		on = .(genes, condition)] %>% 
	rename(sequence = genes) %>% inner_join(annotations)

melted_results <- melted_results[!is.na(FDR)]

melted_results_by_condition <- 
	melted_results[
		gene %like% "Ctrl",
		.(med_LFC = median(LFC)),
		keyby = .(condition)]

setkey(melted_results, condition)

melted_results[, LFC := melted_results_by_condition[
	melted_results, LFC - med_LFC, by = .EACHI]$V1]


################################################################################

median_melted_results <-
	melted_results[, .(
		medLFC = median(LFC),
		FDR = stouffer(FDR)$p),
		by = .(locus_tag, gene, condition)]

median_melted_results <- median_melted_results %>%
	mutate(type = case_when(
		gene == "Ctrl" ~ "control",
		gene != "Ctrl" ~ "knockdown"
	))

################################################################################

setorder(median_melted_results, FDR)

################################################################################
# add fancy names to median melted results

median_melted_results[gene != ".", gene_stylized := paste0("italic('", gene, "')")]
median_melted_results[gene == ".", gene_stylized := paste0("bold('", locus_tag, "')")]
median_melted_results[gene == "", gene_stylized := paste0("bold('", locus_tag, "')")]
median_melted_results[gene == "control", gene_stylized := paste0("bold('", locus_tag, "')")]

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

# CPM_melted <- melt(
# 	data.table(
# 		data_CPM, 
# 		keep.rownames = "spacer"), 
# 	id.vars = "spacer", 
# 	variable.name = "condition", 
# 	value.name = "CPM")
# 
# grouped_CPM <- copy(CPM_melted)
# 
# grouped_CPM[
# 	condition %in% c("Inoculum", "Mouse_P1_003"),  
# 	Condition := 'Inoculum (t0)']
# 
# grouped_CPM[
# 	condition %in% c("Mouse_P1_015", "Mouse_P1_016", "Mouse_P1_017"),  
# 	Condition := 'Lung']
# # plated 10x dilution 
# 
# grouped_CPM[
# 	condition %in% c("Mouse_P1_018", "Mouse_P1_019", "Mouse_P1_020", "Mouse_P1_021", "Mouse_P1_022"),  
# 	Condition := 'Plated ex-vivo 100× dilution']
# 
# grouped_CPM[
# 	condition %in% c("Mouse_P1_007", "Mouse_P1_008"),  
# 	Condition := 'Pelleted ex-vivo 10× dilution']
# 
# grouped_CPM[
# 	condition %in% c("dJMP1", "dJMP2", "dJMP3", "Mouse_P1_006"),  
# 	Condition := 'Inoculum grown on plates']
# 
# grouped_CPM[
# 	condition %in% c("dJMP4", "dJMP5", "Mouse_P1_004"),  
# 	Condition := 'Inoculum (6 gen in-vitro)']
# 
# grouped_CPM <- exp_design[grouped_CPM, on = .(condition)]
# 
# grouped_CPM <- annotated_key[grouped_CPM, on = .(name == spacer)]
# 
# grouped_CPM[!is.na(rep), verbose := paste(media, gDNA_source, growth_condition, rep, sep = "_")]
# 
# grouped_CPM[is.na(rep), verbose := paste(media, gDNA_source, growth_condition, sep = "_")]
# 
# setorder(median_melted_results, locus_tag)
# 
# results_summary <- melted_results[FDR < 0.05, .N, by = .(condition)]
# median_results_summary <- median_melted_results[FDR < 0.05, .N, by = .(condition)]

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
			filter(gene == this_gene | locus_tag == this_gene) %>%
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

################################################################################

median_melted_results %>% 
	mutate(
		Response = case_when(
			medLFC < -1 & FDR < 0.05 ~ "Vulnerable",
			medLFC > 1 & FDR < 0.05 ~ "Resistant",
			TRUE ~ "No Response")) %>%
	select(locus_tag, gene, condition, Response) %>% 
	mutate(condition = case_when(
		condition == "mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0" ~ "In vivo",
		condition == "mouse_plated_10x_inoculum_dilution - LB_plated_6_generations" ~ "In vivo v. in vitro",
		condition == "LB_plated_6_generations - inoculum_pellet_t0" ~ "In vitro")) %>% 
	filter(gene %like% "pur") %>% 
	pivot_wider(id_cols = c(gene), names_from = condition, values_from = Response) %>% 
	arrange(gene)

################################################################################

median_melted_results %>% 
	mutate(
		Response = case_when(
			medLFC < -1 & FDR < 0.05 ~ "Vulnerable",
			medLFC > 1 & FDR < 0.05 ~ "Resistant",
			TRUE ~ "No Response")) %>%
	select(locus_tag, gene, condition, Response) %>% 
	mutate(condition = case_when(
		condition == "mouse_plated_10x_inoculum_dilution - inoculum_pellet_t0" ~ "In vivo",
		condition == "mouse_plated_10x_inoculum_dilution - LB_plated_6_generations" ~ "In vivo v. in vitro",
		condition == "LB_plated_6_generations - inoculum_pellet_t0" ~ "In vitro")) %>% 
	filter(gene %like% "lpt") %>% 
	pivot_wider(id_cols = c(gene), names_from = condition, values_from = Response) %>% 
	arrange(gene)