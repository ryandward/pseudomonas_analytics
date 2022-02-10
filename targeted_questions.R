source("presentation_analysis.R")
# source("Nb_calculator.R")

CPM_melted <-
	melt(
		data.table(data_CPM, keep.rownames = "spacer"),
		id.vars = "spacer",
		variable.name = "condition",
		value.name = "CPM"
	)

annotated_key <- fread("annotated_key.tsv")

targeted_CPM <- CPM_melted[condition %in% c("Inoculum",
																						"Mouse_P1_003",
																						# "Mouse_P1_006",
																						"Mouse_P1_015",
																						"Mouse_P1_016",
																						"Mouse_P1_017")]

targeted_CPM[condition %in% c("Inoculum", "Mouse_P1_003"),  
						 # fancy_condition := 'paste0(bold("pelleted "), "inoculum ", t[0])']
						 Condition := 'Pelleted inoculum t_0']
						 

targeted_CPM[condition == "Mouse_P1_006",  
						 # fancy_condition := 'paste0(bold("plated "), "inoculum ", t[0])']
						 Condition := 'Plated inoculum t_0']
						 
targeted_CPM[condition %in% c("Mouse_P1_015", "Mouse_P1_016", "Mouse_P1_017"),  
						 # fancy_condition := 'paste0(bold("Plated "), italic("ex-vivo "), "10× dilution")']
						 Condition := 'Plated ex-vivo 10× dilution']
						 
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

box_CPM <- function(this_gene) {
	
	to_plot <-
		targeted_CPM[gene_name == this_gene |
								 	locus_tag == this_gene, .(name, gene_name, locus_tag, condition, CPM, verbose, offset, rep, Condition)]
	
	to_plot[, group := gsub("_[0-9]+$", "", verbose) ]
	
	to_plot$offset <- factor(as.character(to_plot$offset),
													 levels = as.character(sort(unique(to_plot$offset))))
	this_plot <-
		ggplot(data = to_plot, aes(x = Condition, y = log2(CPM), fill = offset)) +
		geom_boxplot() + 
		scale_fill_brewer(palette = "Accent") +
		ylab("Log2 Counts per Million") +
		xlab("Condition") +
		ggtitle(paste("Guides recovered in pellet and mouse for", this_gene)) +
		theme(
			plot.title = element_text(hjust = 0.5, size = 20),
			axis.text = element_text(size = 14, color = "black"),
			axis.title = element_text(size = 14, color = "black"),
			legend.text = element_text(size = 8, color = "black"),
			legend.title = element_text(size = 14, color = "black"),
			legend.position = "bottom"
		)
	
	ggthemr("flat")
	print(this_plot)
	# rm(this_gene)
	ggthemr_reset()
	
}

box_CPM_controls <- function() {
	random_controls <-
		targeted_CPM[name %like% "Ctrl", .(name = unique(name))][sample(.N, 8)]
	
	to_plot <-
		targeted_CPM[name %in% random_controls$name, .(name, gene_name, locus_tag, condition, CPM, verbose, offset, rep, Condition)]
	
	to_plot[, group := gsub("_[0-9]+$", "", verbose) ]
	
	to_plot$offset <- factor(as.character(to_plot$offset),
													 levels = as.character(sort(unique(to_plot$offset))))
	this_plot <-
		ggplot(data = to_plot, aes(x = Condition, y = log2(CPM), fill = name)) +
		geom_boxplot() +
		scale_fill_brewer(palette = "Accent") +
		ggtitle(paste("Eight Random Control Guides"))

	
	ggthemr("flat")
	print(this_plot)
	# rm(this_gene)
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
			angle = 0,
			vjust = 1.0,
			hjust = 1
		))
	
	ggthemr("flat")
	print(this_plot)
	ggthemr_reset()
}

# save as 1000 x 750

