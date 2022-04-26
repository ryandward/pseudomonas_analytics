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
	cpm(t(promoter_forensics_matrix)), 
	beside = T,
	legend.text = rownames(t(promoter_forensics_matrix)),
	args.legend = list(x = "top"),
	main = "CPM, promoters found in strains")



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
				fill = 0,
				fun.aggregate = sum)

# remelt to get 1s instead of NAs
additional <-
	melt(additional_grid, id.vars = "spacer", variable.name = "strain", value.name = "count")

additional[, count := count + 1]

additional_grid <- 
	dcast(additional, spacer ~ strain, 
				value.var = "count")

additional_grid_matrix <-
	data.matrix(additional_grid[, -1])

rownames(additional_grid_matrix) <- additional_grid$spacer

additional_grid_cpm <-
	cpm(additional_grid_matrix)

additional_cpm <-
	melt(data.table(additional_grid_cpm, keep.rownames = "spacer"), id.vars = "spacer", variable.name = "strain", value.name = "CPM")

################################################################################

plot(density(log10(additional_cpm[strain == "P1_mfdpir" & spacer %like% "Ctrl"]$CPM)), main = "P1 mating", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional_cpm[strain == "P1_mfdpir" & !spacer %like% "Ctrl"]$CPM)), col = "red")

plot(density(log10(additional_cpm[strain == "P1_PA14" & spacer %like% "Ctrl"]$CPM)), main = "P1 PA14", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional_cpm[strain == "P1_PA14" & !spacer %like% "Ctrl"]$CPM)), col = "red")

################################################################################

plot(density(log10(additional_cpm[strain == "P2_mfdpir" & spacer %like% "Ctrl"]$CPM)), main = "P2 mating", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional_cpm[strain == "P2_mfdpir" & !spacer %like% "Ctrl"]$CPM)), col = "red")

plot(density(log10(additional_cpm[strain == "P2_PA14" & spacer %like% "Ctrl"]$CPM)), main = "P2 PA14", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional_cpm[strain == "P2_PA14" & !spacer %like% "Ctrl"]$CPM)), col = "red")

################################################################################

plot(density(log10(additional_cpm[strain == "P3_mfdpir" & spacer %like% "Ctrl"]$CPM)), main = "P3 mating", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional_cpm[strain == "P3_mfdpir" & !spacer %like% "Ctrl"]$CPM)), col = "red")

plot(density(log10(additional_cpm[strain == "P3_PA14" & spacer %like% "Ctrl"]$CPM)), main = "P3 PA14", col = "blue", ylim = c(0,0.5))
lines(density(log10(additional_cpm[strain == "P3_PA14" & !spacer %like% "Ctrl"]$CPM)), col = "red")