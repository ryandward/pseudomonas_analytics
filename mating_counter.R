additional <- fread(
	"additional_sequences/additional.tsv.gz",
	col.names = c(
		"count",
		"condition",
		"promoter",
		"spacer"))

lost_and_found <- 
	annotated_key[
		additional[
			!spacer %in% melted_results$name, 
			.(spacer = unique(spacer))], 
		on = .( name == spacer )][
			, .(gene_name = unique(gene_name)), 
			by = .(locus_tag)][
				!gene_name %in% melted_results$gene_name]
