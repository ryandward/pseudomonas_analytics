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

interest_CPM <- copy(CPM_melted)

interest_CPM[
	condition %in% c("Inoculum", "Mouse_P1_003"),  
	Condition := 'Pelleted inoculum t_0']

# interest_CPM[
# 	condition == "Mouse_P1_006",
# 	Condition := 'Plated inoculum t_0']

interest_CPM[
	condition %in% c("Mouse_P1_015", "Mouse_P1_016", "Mouse_P1_017"),  
	Condition := 'Plated ex-vivo 10× dilution']

interest_CPM[
	condition %in% c("Mouse_P1_018", "Mouse_P1_019", "Mouse_P1_020", "Mouse_P1_021", "Mouse_P1_022"),  
	Condition := 'Plated ex-vivo 100× dilution']

interest_CPM[
	condition %in% c("Mouse_P1_007", "Mouse_P1_008"),  
	Condition := 'Pelleted ex-vivo 10× dilution']

########################################################################################
# extract information from other data.tables to give context to the Conditions

interest_CPM <- exp_design[interest_CPM, on = .(condition)]

interest_CPM <- annotated_key[interest_CPM, on = .(name == spacer)]

interest_CPM[!is.na(rep), verbose := paste(media, gDNA_source, growth_condition, rep, sep = "_")]

interest_CPM[is.na(rep), verbose := paste(media, gDNA_source, growth_condition, sep = "_")]

# cross-join between unique conditions and the unique groups of spacers
# i.e. conditions of interest /x c(controls, knockdowns)

interest_densities <- data.table(merge(interest_CPM[, unique(condition)], interest_CPM[, unique(type)]))

# the results of a crossjoin are called x and y
# change the names x and y to something more meaningful

setnames(interest_densities, c("x", "y"), c("condition", "type"))

# !! coincidentally, crossjoins give x and y columns and so do density plots
# !! the two groups of x and y columns are COMPLETELY unrelated
# generate NAs in x and y column to create a prototype table

interest_densities[, `:=`(x = NA, y = NA)]

########################################################################################
# quick, thorough, but complex method to generate density plots from all conditions 
# run through every "condition" to generate a points for a density plot

for (i in interest_densities[, unique(condition)]) {
	
	for (j in interest_densities[condition == i, unique(`type`)]) {
		
		density__ <- density(interest_CPM[condition == i & type == j, log2(CPM)], from = -10, to = 15)
		
		density__ <- data.table(condition = i, type = j, x = density__$x, y = density__$y)
		
		interest_densities <- merge(interest_densities, density__, all.y = TRUE, all.x = TRUE)
		
		rm(density__)
	}
}

interest_densities <- interest_densities[complete.cases(interest_densities)]

interest_densities <- 
	interest_densities[interest_CPM[, .(
		Condition = unique(Condition)),
		by = .(condition)], 
		on = .(condition)]

interest_densities_summary <- 
	interest_densities[
		!is.na(Condition), 
		.(
			mean_Y = mean(y),
			sd_Y = pracma::std(y),
			.N,
			sem_Y = pracma::std(y)/sqrt(.N - 1)),
		by = .(type, Condition, x)]

# create a data.table summarizing the conditions across Conditions

interest_densities_summary[is.nan(sd_Y), `:=`(sd_Y = 0, sem_Y = 0)]

interest_densities_summary$CI_lower <- interest_densities_summary$mean_Y + qt((1 - 0.85)/2, df = interest_densities_summary$N - 1)*interest_densities_summary$sem_Y

interest_densities_summary$CI_upper <- interest_densities_summary$mean_Y - qt((1 - 0.85)/2, df = interest_densities_summary$N - 1)*interest_densities_summary$sem_Y

interest_densities_summary[is.na(CI_lower), CI_lower := mean_Y]

########################################################################################

for (i in interest_densities_summary[, unique(Condition)]) {
	
	for (j in interest_densities_summary[, unique(type)]) {

		if (j == "control") {plot_shade = "blue"}
		
		if (j != "control") {plot_shade = "red"}
		
		plot_object <- 
			ggplot(
				interest_densities_summary[Condition == i & type == j], 
				aes(x = x, y = mean_Y)) +
			geom_line(
				data = interest_densities[Condition == i & type == j], 
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
				interest_densities_summary[, min(x)], 
				interest_densities_summary[, max(x)]) +
		ylim(
			interest_densities_summary[, min(CI_lower, na.rm = T)], 
			interest_densities_summary[, max(CI_upper, na.rm = T)]) +
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