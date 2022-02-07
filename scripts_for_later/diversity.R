require('pacman')
source("annotate_key.R")
p_load(vegan)
################################################################################
all_counts <- fread("all_counts_bowtie.tsv")
all_counts <- all_counts[promoter == "P1"][, .(guide, condition, count)]

# accidentally called the spacers "guides", and changing the name here.
all_counts[, spacer := guide]
all_counts[, guide := NULL]

# get rid of all the other promoters

# accidentally called the spacers "guides", and changing the name here.
all_counts[, spacer := guide]
all_counts[, guide := NULL]

################################################################################
data_grid <- data.table::dcast(
	all_counts[condition %like% "Mouse" | condition %like% "Undetermined"], 
	spacer ~ condition,
	value.var = "count", 
	fill = 0)

data_grid_remelted <- melt(data_grid, variable.name = "condition", value.name = "count", id.vars = c('spacer'))
################################################################################

diversity_grid <- data_grid_remelted[spacer %like% "Ctrl", .(ctrl_diversity = diversity(index = "invsimpson", count)), by = .(condition)]
diversity_grid <- data_grid_remelted[!(spacer %like% "Ctrl"), .(non_ctrl_diversity = diversity(index = "invsimpson", count)), by = .(condition)][diversity_grid, on = .(condition)]
diversity_grid <- data_grid_remelted[, .(all_diversity = diversity(index = "invsimpson", count)), by = .(condition)][diversity_grid, on = .(condition)]
clipr::write_clip(diversity_grid)

diversity_list <- melt(diversity_grid, id.vars = "condition", value.name = "count")

ggplot(data = diversity_list, aes(x = condition, y = count, fill = variable)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	scale_fill_brewer(palette = "Paired") +
	ggtitle(paste("Shannon Diversity Indices")) +
	theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1))

