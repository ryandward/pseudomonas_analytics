source("counter_salvage.R")
# source("Nb_calculator.R")

CPM_melted <-
	melt(
		data.table(data_CPM, keep.rownames = "spacer"),
		id.vars = "spacer",
		variable.name = "condition",
		value.name = "CPM"
	)

annotated_key <- fread("annotated_key.tsv")

targeted_CPM <- CPM_melted[condition %in% c("Mouse_P1_003",
																						"Mouse_P1_006",
																						"Mouse_P1_015",
																						"Mouse_P1_016",
																						"Mouse_P1_017")]

targeted_CPM <- exp_design[targeted_CPM, on = .(condition)]
targeted_CPM <- annotated_key[targeted_CPM, on = .(name == spacer)]

targeted_CPM[!is.na(rep), verbose := paste(media, gDNA_source, growth_condition, rep, sep = "_")]

targeted_CPM[is.na(rep), verbose := paste(media, gDNA_source, growth_condition, sep = "_")]

plot_CPM <- function(this_gene) {
	to_plot <-
		targeted_CPM[gene_name == this_gene |
								 	locus_tag == this_gene, .(name, gene_name, locus_tag, condition, CPM, verbose, offset)]
	
	to_plot$offset <- factor(as.character(to_plot$offset),
													 levels = as.character(sort(unique(to_plot$offset))))
	
	this_plot <-
		ggplot(data = to_plot, aes(x = verbose, y = CPM, fill = offset)) +
		geom_bar(stat = "identity", position = position_dodge()) +
		ggtitle(paste("CPM for guides:", this_gene)) +
		theme(axis.text.x = element_text(
			angle = 55,
			vjust = 1.0,
			hjust = 1
		))
	
	ggthemr("flat")
	print(this_plot)
	rm(this_gene)
	ggthemr_reset()
	
}

# save as 1000 x 750

plot_CPM_error <- function(this_gene) {
	
	this_gene <- "orfN"
	to_plot <-
		targeted_CPM[gene_name == this_gene |
								 	locus_tag == this_gene, .(name, gene_name, locus_tag, condition, CPM, verbose, offset, rep)]
	
	to_plot[, group := gsub("_[0-9]+$", "", verbose) ]
	
	this_plot <-
		ggplot(data = to_plot, aes(x = verbose, y = CPM, fill = offset)) +
		geom_boxplot() + ggtitle("Box plot") +
		ggtitle(paste("CPM for guides:", this_gene)) +
		theme(axis.text.x = element_text(
			angle = 55,
			vjust = 1.0,
			hjust = 1
		))
	
	ggthemr("flat")
	print(this_plot)
	rm(this_gene)
	ggthemr_reset()
	
}


plot_CPM_controls <- function() {
	random_controls <-
		CPM_melted[spacer %like% "Ctrl", .(spacer = unique(spacer))][sample(.N, 8)]
	
	to_plot <-
		targeted_CPM[name %in% random_controls$spacer, .(name, gene_name, locus_tag, condition, CPM, verbose, offset)]
	
	to_plot$offset <- factor(as.character(to_plot$offset),
													 levels = as.character(sort(unique(to_plot$offset))))
	
	this_plot <-
		ggplot(data = to_plot, aes(x = verbose, y = CPM, fill = name)) +
		geom_bar(stat = "identity", position = position_dodge()) +
		ggtitle(paste("8 Random Control Guides. CPM")) +
		theme(axis.text.x = element_text(
			angle = 55,
			vjust = 1.0,
			hjust = 1
		))
	
	ggthemr("flat")
	print(this_plot)
	ggthemr_reset()
}

# save as 1000 x 750

