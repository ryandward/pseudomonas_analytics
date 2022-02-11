p_load("ggthemr")

problem_guides <- melted_results[FDR>0.9 & type != "control" & condition == "mouse_plated_10x_inoculum_dilution - inoculum_plated_t0"][, .(.N), by = .(gene_name)]
# to make a barplot
this_gene <- "orfN"

ggthemr("flat")

to_plot <- melted_results[gene_name == this_gene, .(name, gene_name, offset = as.character(offset), LFC, FDR, condition)]
ggplot(data = to_plot, aes(x = condition, y = -log10(FDR), fill = offset)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	ggtitle(paste("-log10(FDR) for guides in", this_gene)) +
	theme(axis.text.x = element_text(angle = 35, vjust = 1.0, hjust = 1)) +
	coord_flip()

	
ggthemr("dust")

to_plot <- melted_results[gene_name == this_gene, .(name, gene_name, offset = as.character(offset), LFC, FDR, condition)]
ggplot(data = to_plot, aes(x = offset, y = -log10(FDR), fill = condition)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	ggtitle(paste("-log10(FDR) for guides in", this_gene)) +
	theme(axis.text.x = element_text(angle = 35, vjust = 1.0, hjust = 1)) 

ggthemr_reset()


CPM_var <- melt(data.table(data_CPM, keep.rownames = "spacer"), id.vars = "spacer", variable.name = "condition", value.name = "CPM") 
CPM_var <- annotated_key[CPM_var, on = .(name == spacer)]

CPM_by_condition <- CPM_var[, .(condition_MED = var(CPM)), by = .(condition, locus_tag, gene_name)]
CPM_by_spacer <- CPM_var[, .(spacer_MED = var(CPM)), by = .(name, locus_tag, gene_name)]
bad_guides <- CPM_by_spacer[CPM_by_condition, on = .(locus_tag, gene_name), allow.cartesian = TRUE]

zzz <- bad_guides[, .(std_cond = std(condition_MED)), by = .(locus_tag, gene_name)]
zzz <- zzz[bad_guides[, .(std_spacer = std(spacer_MED)), by = .(locus_tag, gene_name)], on = .(locus_tag, gene_name)]

melt
