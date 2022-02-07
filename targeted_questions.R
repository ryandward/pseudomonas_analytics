source("counter_salvage.R")
source("Nb_calculator.R")

annotated_key <- fread("annotated_key.tsv")

targeted_CPM <- CPM_melted[condition %in% c(
	"Mouse_P1_003",
	"Mouse_P1_006",
	"Mouse_P1_015", 
	"Mouse_P1_016", 
	"Mouse_P1_017")]

targeted_CPM <- exp_design[targeted_CPM, on = .(condition)]
targeted_CPM <- annotated_key[targeted_CPM, on = .(name == spacer)]

targeted_CPM[!is.na(rep), verbose := paste(
	media, gDNA_source, growth_condition, rep, sep = "_")]

targeted_CPM[is.na(rep), verbose := paste(
	media, gDNA_source, growth_condition, sep = "_")]

this_gene <- "rnhA"


to_plot <- targeted_CPM[gene_name == this_gene | locus_tag == this_gene, .(name, gene_name, locus_tag, condition, CPM, verbose)]
this_plot <- ggplot(data = to_plot, aes(x = verbose, y = CPM, fill = name)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	ggtitle(paste("CPM for guides:", this_gene)) +
	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1)) 

print(this_plot)
rm(this_gene)
# save as 1000 x 750

