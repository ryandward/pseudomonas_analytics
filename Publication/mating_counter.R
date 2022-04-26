library(pacman)
p_load(data.table, viridis, tidyverse)
p_load_current_gh("hrbrmstr/hrbrthemes")

 
# annotated_key <- fread("Publication/annotated_key.tsv")
# 
# additional <- fread(
# 	"additional_sequences/additional.tsv.gz",
# 	col.names = c(
# 		"count",
# 		"condition",
# 		"promoter",
# 		"spacer"))
# 
# 
# inoculum_counts <- 
# 	fread(
# 		"Publication/inoculum_counts.tsv.gz", 
# 		header = FALSE, 
# 		col.names = c(
# 			"promoter", 
# 			"spacer", 
# 			"count"))
# 
# 
# additional[, `:=`(c("claimed_promoter", "strain"), tstrsplit(condition, "_"))]
# 
# promoter_forensics <-
# 	dcast(additional, claimed_promoter + strain  ~ promoter, value.var = "count", fun.aggregate = sum)
# 
# promoter_forensics_matrix <- data.matrix(promoter_forensics[, -c(1:2)])
# 
# rownames(promoter_forensics_matrix) <- promoter_forensics[, paste(claimed_promoter, strain)]
# 
# barplot(
# 	cpm(t(promoter_forensics_matrix)), 
# 	beside = T,
# 	legend.text = rownames(t(promoter_forensics_matrix)),
# 	args.legend = list(x = "top"),
# 	main = "CPM, promoters found in strains")
# 
# 
# #There's a problem with promoter P3 contamination in the sample that supposedly
# #only contains P1-PA14.
# 
# #Upon further inspection, this contains the same level of contamination as 
# #the second "inoculum", which suggests it's identical.
# 
# #i.e. these two values are identical
# print(inoculum_counts[, .(sum(count)), by = .(promoter)])
# print(additional[condition == "P1_PA14", .(sum(count)), by = .(promoter)])
# 
# #Upon further inspection both of these correspond to files called:
# #"P1_PA14_S26_R1_001.fastq" and "P1_PA14_S26_R2_001.fastq.gz"
# #Furthermore, all files have the same md5sum signature
# 
# #It was suggested by collaborators to use Mouse_P1_001 from the original dataset
# #instead.
# 
# # BEGIN THE PROCESS FROM THE TOP AGAIN...
# ################################################################################

additional <- fread(
	"additional_sequences/additional.tsv.gz",
	col.names = c(
		"count",
		"condition",
		"promoter",
		"spacer"))

counts.batch1 <- 
	fread(
		"Publication/all_counts_seal.tsv.gz", 
		header = FALSE, 
		col.names = c(
			"promoter", 
			"spacer", 
			"count", 
			"condition"))

inoculum_counts <- 
	fread(
		"Publication/inoculum_counts.tsv.gz", 
		header = FALSE, 
		col.names = c(
			"promoter", 
			"spacer", 
			"count"))

inoculum_counts <- inoculum_counts %>% 
	mutate (condition = "Inoculum")

counts.batch1 <- counts.batch1 %>% 
	bind_rows(inoculum_counts)

# Remove P1_PA14 counts, replace with Mouse_P1_001 counts, and rename to be more 
# descriptive: PA_PA14 (round 1). Then separate into two columns 
# "claimed_promoter" and "strain", which are separated by "_".
additional <- additional %>% 
	filter(condition != "P1_PA14") %>%
	bind_rows(counts.batch1 %>% filter(condition == "Mouse_P1_001")) %>%
	mutate(
		condition = replace(
			condition, 
			condition == "Mouse_P1_001", 
			"P1_PA14 (round 1)")) %>%
	separate(
		condition, 
		into = c("claimed_promoter", "strain"), 
		sep = "_")

# Perform forensics to see the composition of the strain-promoter combos
promoter_forensics <- additional %>% 
	group_by(claimed_promoter, promoter, strain) %>% 
	summarise (count = sum(count)) %>%
	pivot_wider(id_cols = c( "strain", "claimed_promoter"), names_from = promoter, values_from = count) %>%
	unite ("claimed_strain", claimed_promoter, strain) 

# Counts per million by each (promoter, claimed_promoter+strain) combo, then 
# magrittr into ggplot and bar plot.
forensics.plot <- promoter_forensics %>% 
	select(-claimed_strain) %>% 
	data.matrix %>%
	t %>%
	cpm %>%
	t %>%
	as_tibble %>%
	bind_cols (promoter_forensics %>% select(claimed_strain)) %>%
	pivot_longer(!claimed_strain, names_to = "promoter", values_to = "CPM") %>%
	ggplot(aes(fill = promoter, x = claimed_strain, y = CPM)) +
	geom_bar(position = "dodge", stat = "identity") +
	theme_ipsum() + 
	scale_fill_viridis(discrete = TRUE, option = "cividis") + 
	ggtitle("Promoters (P1, P2, and P3) found in each strain (CPM)") +
	scale_y_continuous(labels = label_number_si())

print(forensics.plot)

# Now that we know the vast majority of each promoter encountered is in fact
# the claimed promoter, we can simply remove all instances where the promoter
# obtained from sequencing does not match the claimed promoter.
additional <- additional %>% filter(promoter == claimed_promoter)

################################################################################

# Use edgeR to determine log-fold changes
source("Publication/publication_counter.R")

lost_and_found <- 
	annotated_key[
		additional[strain == "mfdpir" &
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

additional_cpm[spacer %like% "Ctrl", type := "control"]

additional_cpm[!spacer %like% "Ctrl", type := "knockdown"]

additional_cpm[, Type := factor(type, levels = unique(type))]

additional_cpm[, strain := factor(strain, levels = unique(strain))]

################################################################################

for(i in unique(additional_cpm$strain)){
	
	cpm.plot <-
		additional_cpm %>% filter (strain == i) %>%
			ggplot(aes(x = CPM, fill = Type)) +
			geom_density() +
			theme_ipsum() +
			ggtitle(this_title) +
			scale_fill_viridis(discrete = TRUE, alpha = 0.5, option = "cividis", direction = -1) +
			scale_colour_manual("black") +
			scale_x_continuous(
				trans = scales::pseudo_log_trans(base = 10),
				breaks = c(0, 10^(1:6)),
				labels = label_number_si(),
				limits = c(0, max(additional_cpm$CPM))) +
			theme(legend.position = "bottom") + 
		ggtitle(i)
	
	print(cpm.plot)
}
