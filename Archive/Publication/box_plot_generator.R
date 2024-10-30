source("Publication/publication_counter.R")

p_

##################################################################################
# TAKE A SUB-SELECTION OF THE GROUPED CPM DATA THAT WE ARE INTERSETED IN EXAMINING
targeted_CPM <- grouped_CPM[
	condition %in% c(
		"Mouse_P1_015", 
		"Mouse_P1_016", 
		"Mouse_P1_017", 
		"Inoculum", 
		"Mouse_P1_003",
		"Mouse_P1_004",
		"dJMP4")]  

# save as 1000 x 750

box_CPM <- function(this_gene) {
	
	to_plot <-
		targeted_CPM[
			gene_name == this_gene |
				locus_tag == this_gene, .(name, gene_name, locus_tag, condition, CPM, verbose, offset, rep, Condition)]
	
	to_plot[
		, group := gsub("_[0-9]+$", "", verbose) ]
	
	to_plot$offset <- factor(
		as.character(to_plot$offset),
		levels = as.character(sort(unique(to_plot$offset))))
	
	this_plot <-
		ggplot(data = to_plot, aes(x = Condition, y = CPM, fill = offset)) +
		geom_boxplot() + 
		ylab("Counts per Million") +
		xlab("Condition") +
		ggtitle(paste("Guides recovered for", this_gene)) +
		theme(
			plot.title = element_text(hjust = 0.5, size = 20),
			axis.text = element_text(size = 14, color = "black"),
			axis.title = element_text(size = 14, color = "black"),
			legend.text = element_text(size = 8, color = "black"),
			legend.title = element_text(size = 14, color = "black"),
			legend.position = "bottom"
		) +
		theme_ipsum() +
		scale_fill_viridis(discrete = TRUE, alpha = 0.8) +
		scale_y_continuous(
			trans = scales::pseudo_log_trans(base = 10),
			breaks = c(0, 10^(1:6)),
			labels = label_number_si(),
			limits = c(0, max(to_plot$CPM)))
	
	print(this_plot)
	rm(this_gene)

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

	
	print(this_plot)
	# rm(this_gene)

}


bar_CPM_controls <- function() {
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
	
	print(this_plot)
}

# save as 1000 x 750

