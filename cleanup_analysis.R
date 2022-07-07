library(pacman)

p_load(data.table, tidyverse, edgeR, poolr, scales, pheatmap)

count_stats <- fread("stats.tsv", col.names = c("sequence", "count", "condition"))

all_guides <- fread("all_guides.tsv") %>% 
	group_by(type) %>% 
	mutate(
		index = 1:n(),
		gene = case_when(type == "control" ~ "control", gene == "-" ~ locus_tag, type == "knockdown" ~ gene),
		guide = case_when(type == "control" ~ paste(gene, index), type != "control" ~ paste(gene, offset)))

all_guides <- all_guides %>%
	filter(!sequence %in% c(
		"CCTTGATGCTGTTGAGGATC",
		"TGGTCTGGGTGCGCTCGGAG",
		"TAGTCAAAGATACCCCCGAA",
		"GCTCGCGGTTTACTTCGACC"
	))

map_stats <- fread("map_stats.tsv")

exp_design <- fread("unified_counts/exp_design.tsv")

exp_design <- exp_design %>% inner_join(map_stats)

exp_design <- exp_design %>% 
	#this one didn't get infected
	filter(!growth_condition %like% "100x_inoculum_dilution") %>%
	filter(!condition %in% c("dJMP3", "dJMP5")) %>%
	# filter(!condition %in% c("dJMP1")) %>%
	
	data.table

exp_design <- exp_design %>% mutate(sample_group = paste(gDNA_source, growth_condition, media, sep = "_"))
 
exp_design <- exp_design %>%
	mutate(match_prop = matched/reads)

quality_decision <- count_stats %>% 
	inner_join(all_guides) %>% 
	pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0) %>% 
	pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>% 
	group_by(condition) %>% 
	mutate(cpm = 1e6*(count/sum(count))) %>% 
	group_by(condition, type) %>% 
	summarise(cpm = sum(cpm)) %>% 
	pivot_wider(id_cols = condition, names_from = type, values_from = cpm) %>%
	arrange(desc(focused)) %>% inner_join(exp_design) %>% 
	select(condition, control, focused, knockdown, media, gDNA_source, match_prop, reads)

exp_design <- exp_design %>% 
	inner_join(quality_decision) %>% 
	filter(match_prop > 0.85 & reads > 5e6) %>%
	group_by(sample_group) %>% 
	mutate(rep = as.character(1:n())) 

exp_design <- exp_design %>% mutate(sample_name = paste(gDNA_source, growth_condition, media , rep, sep = "_"))

exp_design <- exp_design %>% filter(sample_group %in% (
	exp_design %>% group_by(sample_group) %>% summarise(reps = max(rep)) %>% filter(reps > 1) %>% pull(sample_group)))

exp_design <- exp_design %>% data.table

################################################################################

count_stats %>%
	inner_join(all_guides) %>%
	filter(type != "focused") %>%
	group_by(sequence, condition) %>% 
	inner_join(exp_design) %>% 
	pivot_wider(id_cols = sequence, names_from = sample_name, values_fill = 0, values_from = count) %>% 
	ungroup %>% 
	select(-sequence) %>% 
	data.matrix %>% 
	cor %>% 
	round(4) %>% 
	pheatmap(
		cutree_rows = 4,
		cutree_cols = 4,
		display_numbers = T, 
		number_format = "%.3f",
		breaks = seq(0.5, 1, length.out = 9999),
		color=c(colorRampPalette(c("navy","white", "red"))(9997), "gray"))

################################################################################

count_stats <- count_stats %>% filter(condition %in% exp_design$condition)

setorder(exp_design, condition)
setorder(count_stats, condition)


count_stats.mat <- count_stats %>% 
	inner_join(all_guides) %>% 
	pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0) 

count_stats.mat %>% 
	pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
	inner_join(all_guides)

count_stats <- count_stats.mat %>% 
	pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>%
	inner_join(all_guides) %>% 
	unite("guide", c(gene, offset), remove = F)

##########################################################################################
# which kind of gDNA is better? plates.
count_stats.mat %>% 
	select(-type, -sequence) %>% 
	data.matrix %>% 
	cpm %>% 
	as_tibble %>% 
	cbind(count_stats.mat %>% select(type, sequence)) %>% 
	pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "cpm") %>% 
	inner_join(exp_design) %>% 
	unite("growth_condition", c(growth_condition, media)) %>%
	ggplot(aes(x = cpm)) + 
	geom_density(aes(fill = rep), alpha = 0.25) + 	
	scale_x_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^(1:6)),
		labels = label_number_si()) + 
	facet_grid(facets = c("gDNA_source", "growth_condition")) +
	ggtitle("Counts per Million") -> plot_object

print(plot_object)

##########################################################################################
# is there a difference in the composition of the knockdown and control components
# only use plated non-focused

count_stats.mat.quality <- 
	count_stats.mat %>% 
	pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>% 
	inner_join(exp_design) %>% 
	# filter(gDNA_source == "plated") %>% 
	filter(type != "focused") %>%	
	pivot_wider(id_cols = c(type, sequence), names_from = condition, values_from = count, values_fill = 0) 

count_stats.mat.quality %>% 
	select(-type, -sequence) %>% 
	data.matrix %>% 
	cpm %>% 
	as_tibble %>% 
	cbind(count_stats.mat.quality %>% select(type, sequence)) %>% 
	pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "cpm") %>% 
	inner_join(exp_design) %>% 
	unite("growth_condition", c(growth_condition, media)) %>%
	ggplot(aes(x = cpm)) + 
	geom_density(aes(fill = rep), alpha = 0.25) + 	
	scale_x_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^(1:6)),
		labels = label_number_si()) + 
	facet_grid(facets = c("type", "sample_group")) +
	ggtitle("Counts per Million") -> plot_object

print(plot_object)

##########################################################################################

count_stats.quality <- 
	count_stats.mat.quality %>% 
	pivot_longer(cols = -c(type, sequence), names_to = "condition", values_to = "count") %>% 
	inner_join(all_guides) %>% unite("guide", c(gene, offset), remove = F) %>%
	mutate(guide = case_when(type == "control" ~ paste(gene, index), type == "knockdown" ~ paste(gene,offset)))
	
setorder(count_stats.quality, condition)
setorder(exp_design, condition)

conditions.quality <-
	count_stats.quality %>% 
	select(condition) %>% 
	unique

##########################################################################################


data_grid <- count_stats.quality %>% 
	pivot_wider(id_cols = sequence, names_from = condition, values_from = count)

data_grid_matrix <- data_grid %>% 
	select(-sequence) %>% data.matrix

row.names(data_grid_matrix) <- data_grid$sequence

colnames(data_grid_matrix) <- data_grid %>% 
	colnames %>% 
	data.table(condition = .) %>% 
	inner_join(exp_design) %>% 
	pull(sample_name)

data_group <-
	factor(
		exp_design[condition %in% conditions.quality$condition, sample_group],
		levels = unique(unique(exp_design[condition %in% conditions.quality$condition, sample_group])))

data_permut <- model.matrix(~ 0 + data_group)

colnames(data_permut) <- levels(data_group)

rownames(data_permut) <- colnames(data_grid_matrix)

data_permut_check <-
	melt(
		data.table(
			data_permut, 
			keep.rownames = "sample_name"), 
		id.vars = "sample_name") %>% 
	filter(value == 1) %>% 
	select(sample_name) %>% 
	inner_join(exp_design)

print(data_permut_check)

##########################################################################################

data_y <- DGEList(
	counts = data_grid_matrix,
	group = data_group,
	genes = row.names(data_grid_matrix))

data_keep <- filterByExpr(data_grid_matrix, data_group)

data_y <- data_y[data_keep, , keep.lib.sizes = FALSE]

data_y <- calcNormFactors(data_y)

data_y <- estimateDisp(data_y, data_permut)

data_fit <- glmFit(data_y, data_permut, robust = T)

data_CPM <- cpm(data_y, prior.count = 0)

################################################################################

contrast_levels <-
	c("plated_6_generations_LB - plated_t0_inoculum",
		"plated_10x_inoculum_dilution_mouse - plated_6_generations_LB",
		"plated_10x_inoculum_dilution_mouse - plated_t0_inoculum")

data_contrast <- makeContrasts(
	contrasts = contrast_levels,
	levels = data_permut)

################################################################################

results_FDR <- data.table(count_stats.quality)[, .(genes = unique(sequence))]

results_LFC <- data.table(count_stats.quality)[, .(genes = unique(sequence))]

################################################################################

for (i in 1:ncol(data_contrast)) {
	
	results <- glmLRT(data_fit, contrast = data_contrast[, i])
	
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
	rename(sequence = genes) %>% inner_join(all_guides)

melted_results <- melted_results[!is.na(FDR)]

melted_results_by_condition <- 
	melted_results[
		type == "control",
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

#################################################################################
