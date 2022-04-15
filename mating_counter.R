additional <- fread(
	"additional_sequences/additional.tsv.gz",
	col.names = c(
		"count",
		"condition",
		"promoter",
		"spacer"))

additional[, `:=`(c("claimed_promoter", "strain"), tstrsplit(condition, "_"))]


promoter_forensics <-
	dcast(additional, claimed_promoter + strain  ~ promoter, value.var = "count", fun.aggregate = sum)

promoter_forensics_matrix <- data.matrix(promoter_forensics[, -c(1:2)])

rownames(promoter_forensics_matrix) <- promoter_forensics[, paste(claimed_promoter, strain)]

barplot(
	t(promoter_forensics_matrix), 
	beside = T,
	legend.text = rownames(t(promoter_forensics_matrix)),
	args.legend = list(x = "top"))



additional <- additional[claimed_promoter == promoter]

lost_and_found <- 
	annotated_key[
		additional[condition %like% "mfdpir" &
			!spacer %in% melted_results$name, 
			.(spacer = unique(spacer))], 
		on = .( name == spacer )][
			, .(gene_name = unique(gene_name)), 
			by = .(locus_tag)][
				!gene_name %in% melted_results$gene_name]

additional_grid <- 
	dcast(additional, spacer ~ promoter + strain, 
				value.var = "count",
				fill = 0)

# remelt to get zeroes instead of NAs
additional <-
	melt(additional_grid, id.vars = "spacer", variable.name = "strain", value.name = "count")

additional[, count := count + 1]

additional_grid <- 
	dcast(additional, spacer ~ strain, 
				value.var = "count",
				fill = 0)

################################################################################

plot(density(log10(additional[strain == "P1_mfdpir" & spacer %like% "Ctrl"]$count)), main = "P1 mating", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional[strain == "P1_mfdpir" & !spacer %like% "Ctrl"]$count)), col = "red")

plot(density(log10(additional[strain == "P1_PA14" & spacer %like% "Ctrl"]$count)), main = "P1 PA14", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional[strain == "P1_PA14" & !spacer %like% "Ctrl"]$count)), col = "red")

################################################################################

plot(density(log10(additional[strain == "P2_mfdpir" & spacer %like% "Ctrl"]$count)), main = "P2 mating", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional[strain == "P2_mfdpir" & !spacer %like% "Ctrl"]$count)), col = "red")

plot(density(log10(additional[strain == "P2_PA14" & spacer %like% "Ctrl"]$count)), main = "P2 PA14", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional[strain == "P2_PA14" & !spacer %like% "Ctrl"]$count)), col = "red")

################################################################################

plot(density(log10(additional[strain == "P3_mfdpir" & spacer %like% "Ctrl"]$count)), main = "P3 mating", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional[strain == "P3_mfdpir" & !spacer %like% "Ctrl"]$count)), col = "red")

plot(density(log10(additional[strain == "P3_PA14" & spacer %like% "Ctrl"]$count)), main = "P3 PA14", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional[strain == "P3_PA14" & !spacer %like% "Ctrl"]$count)), col = "red")