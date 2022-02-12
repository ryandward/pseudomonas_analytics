source("presentation_analysis.R")

p_load(pracma, data.table, ggplot2, ggthemr)

CPM_melted <- melt(
	data.table(
		data_CPM, 
		keep.rownames = "spacer"), 
	id.vars = "spacer", 
	variable.name = "condition", 
	value.name = "CPM")

########################################################################################

# !! assign a "Condition" to each "condition"
# conditions correspond to the unique value of each experimental sample
# Conditions correspond to the "group of conditions", i.e. no replicate indicator
# "Condition"s will be used when referring to the group of "conditions" in plots

# For future experiments, consider "sample" vs "condition".

grouped_CPM <- copy(CPM_melted)

grouped_CPM[
	condition %in% c("Inoculum", "Mouse_P1_003"),  
	Condition := 'Pelleted inoculum t_0']

# grouped_CPM[
# 	condition == "Mouse_P1_006",
# 	Condition := 'Plated inoculum t_0']

grouped_CPM[
	condition %in% c("Mouse_P1_015", "Mouse_P1_016", "Mouse_P1_017"),  
	Condition := 'Plated ex-vivo 10× dilution']

grouped_CPM[
	condition %in% c("Mouse_P1_018", "Mouse_P1_019", "Mouse_P1_020", "Mouse_P1_021", "Mouse_P1_022"),  
	Condition := 'Plated ex-vivo 100× dilution']

grouped_CPM[
	condition %in% c("Mouse_P1_007", "Mouse_P1_008"),  
	Condition := 'Pelleted ex-vivo 10× dilution']

########################################################################################
# extract information from other data.tables to give context to the Conditions

grouped_CPM <- exp_design[grouped_CPM, on = .(condition)]

grouped_CPM <- annotated_key[grouped_CPM, on = .(name == spacer)]

grouped_CPM[!is.na(rep), verbose := paste(media, gDNA_source, growth_condition, rep, sep = "_")]

grouped_CPM[is.na(rep), verbose := paste(media, gDNA_source, growth_condition, sep = "_")]

# cross-join between unique conditions and the unique groups of spacers
# i.e. conditions of grouped /x c(controls, knockdowns)

grouped_densities <- data.table(merge(grouped_CPM[, unique(condition)], grouped_CPM[, unique(type)]))

# the results of a crossjoin are called x and y
# change the names x and y to something more meaningful

setnames(grouped_densities, c("x", "y"), c("condition", "type"))

# !! coincidentally, crossjoins give x and y columns and so do density plots
# !! the two groups of x and y columns are COMPLETELY unrelated
# generate NAs in x and y column to create a prototype table

grouped_densities[, `:=`(x = NA, y = NA)]

########################################################################################
# quick, thorough, but complex method to generate density plots from all conditions 
# run through every "condition" to generate a points for a density plot

for (i in grouped_densities[, unique(condition)]) {
	
	for (j in grouped_densities[condition == i, unique(`type`)]) {
		
		density__ <- density(grouped_CPM[condition == i & type == j, log2(CPM)], from = -10, to = 15)
		
		density__ <- data.table(condition = i, type = j, x = density__$x, y = density__$y)
		
		grouped_densities <- merge(grouped_densities, density__, all.y = TRUE, all.x = TRUE)
		
		rm(density__)
	}
}

grouped_densities <- grouped_densities[complete.cases(grouped_densities)]

grouped_densities <- 
	grouped_densities[grouped_CPM[, .(
		Condition = unique(Condition)),
		by = .(condition)], 
		on = .(condition)]

grouped_densities_summary <- 
	grouped_densities[
		!is.na(Condition), 
		.(
			mean_Y = mean(y),
			sd_Y = pracma::std(y),
			.N,
			sem_Y = pracma::std(y)/sqrt(.N - 1)),
		by = .(type, Condition, x)]

# create a data.table summarizing the conditions across Conditions

grouped_densities_summary[is.nan(sd_Y), `:=`(sd_Y = 0, sem_Y = 0)]

grouped_densities_summary$CI_lower <- grouped_densities_summary$mean_Y + qt((1 - 0.85)/2, df = grouped_densities_summary$N - 1)*grouped_densities_summary$sem_Y

grouped_densities_summary$CI_upper <- grouped_densities_summary$mean_Y - qt((1 - 0.85)/2, df = grouped_densities_summary$N - 1)*grouped_densities_summary$sem_Y

grouped_densities_summary[is.na(CI_lower), CI_lower := mean_Y]

########################################################################################

for (i in grouped_densities_summary[, unique(Condition)]) {
	
	for (j in grouped_densities_summary[, unique(type)]) {

		if (j == "control") {plot_shade = "blue"}
		
		if (j != "control") {plot_shade = "red"}
		
		plot_object <- 
			ggplot(
				grouped_densities_summary[Condition == i & type == j], 
				aes(x = x, y = mean_Y)) +
			geom_line(
				data = grouped_densities[Condition == i & type == j], 
				aes(x = x, y = y, group = condition), 
				color = "dark grey") +
			geom_line(
				size = 1.5, 
				alpha = 0.8, 
				color = "black") +
			geom_ribbon(
				aes(ymin = CI_lower, ymax = CI_upper),
				fill = plot_shade, 
				alpha = 0.2) + 
			xlim(
				grouped_densities_summary[, min(x)], 
				grouped_densities_summary[, max(x)]) +
		ylim(
			grouped_densities_summary[, min(CI_lower, na.rm = T)], 
			grouped_densities_summary[, max(CI_upper, na.rm = T)]) +
			ggtitle(
				bquote(Log[2] ~ counts ~ per ~ million. ~ bold(.(i)) ~ italic(.(j)) ~ guides. )) +
			xlab(
				bquote(Log[2] ~ counts ~ per ~ million)) +
			ylab("Density")
		
		ggthemr("flat")
		print(plot_object)
		ggthemr_reset()
		
	}
}