source("Nb_calculator.R")
source("interest_questions.R")
########################################################################################
##### Simple but complex way to get density plots from all conditions into one data.table




interest_CPM <- copy(CPM_melted)

interest_CPM[condition %in% c("Inoculum", "Mouse_P1_003"),  
						 Condition := 'Pelleted inoculum t_0']

interest_CPM[condition == "Mouse_P1_006",  
						 Condition := 'Plated inoculum t_0']

interest_CPM[condition %in% c("Mouse_P1_015", "Mouse_P1_016", "Mouse_P1_017"),  
						 Condition := 'Plated ex-vivo 10× dilution']

interest_CPM[condition %in% c("Mouse_P1_018", "Mouse_P1_019", "Mouse_P1_020", "Mouse_P1_021", "Mouse_P1_022"),  
						 Condition := 'Plated ex-vivo 100× dilution']

interest_CPM <- exp_design[interest_CPM, on = .(condition)]
interest_CPM <- annotated_key[interest_CPM, on = .(name == spacer)]

interest_CPM[!is.na(rep), verbose := paste(media, gDNA_source, growth_condition, rep, sep = "_")]

interest_CPM[is.na(rep), verbose := paste(media, gDNA_source, growth_condition, sep = "_")]



interest_densities <- data.table(merge(interest_CPM[, unique(condition)], interest_CPM[, unique(type)]))

setnames(interest_densities, c("x", "y"), c("condition", "type"))

interest_densities[, `:=`(x = NA, y = NA)]

for (i in interest_densities[, unique(condition)]) {
	for (j in interest_densities[condition == i, unique(`type`)]) {
		
		density__ <- density(interest_CPM[condition == i & type == j, log2(CPM)], from = -10, to = 15)
		
		density__ <- data.table(condition = i, type = j, x = density__$x, y = density__$y)
		
		interest_densities <- merge(interest_densities, density__, all.y = TRUE, all.x = TRUE)
		
	}
}

interest_densities <- interest_densities[complete.cases(interest_densities)]

interest_densities <- interest_densities[interest_CPM[, .(Condition = unique(Condition)), by = .(condition)], on = .(condition)]

interest_densities_summary <- interest_densities[
	, .(
		meany = mean(y),
		stdy = std(y),
		.N,
		semy = std(y)/sqrt(.N - 1)
	), 
	by = .(type, Condition, x)]

interest_densities_summary[is.nan(stdy), `:=`(stdy = 0, semy = 0)]

interest_densities_summary$CI_lower <- interest_densities_summary$meany + qt((1-0.85)/2, df=interest_densities_summary$N-1)*interest_densities_summary$semy
interest_densities_summary$CI_upper <- interest_densities_summary$meany - qt((1-0.85)/2, df=interest_densities_summary$N-1)*interest_densities_summary$semy

interest_densities_summary[is.na(CI_lower), CI_lower := meany]
########################################################################################

for (i in interest_densities_summary[, unique(Condition)]) {
	for (j in interest_densities_summary[, unique(type)]) {
		
		ggthemr("flat")
		
		if (j == "control") {plot_shade = "blue"}
		if (j != "control") {plot_shade = "red"}
		
		plot_object <- ggplot(
		interest_densities_summary[Condition == i & type == j], 
		aes(x = x, y = meany)) +
		geom_line(
			data = interest_densities[Condition == i & type == j], 
			aes(x = x, y = y, group = condition), 
			color = "dark grey") +
		geom_line(
			size = 1.5, 
			alpha = 0.8, 
			color = "black") +
		geom_ribbon(
			aes(
				ymin = CI_lower, ymax = CI_upper),
			fill = plot_shade, 
			alpha = 0.2) +
		xlim(interest_densities_summary[, min(x)], interest_densities_summary[, max(x)]) +
		ylim(interest_densities_summary[, min(CI_lower, na.rm = T)], interest_densities_summary[, max(CI_upper, na.rm = T)]) +
		ggtitle(bquote(Log[2] ~ counts ~ per ~ million. ~ bold(.(i)) ~ italic(.(j)) ~ guides. )) +
		xlab(bquote(Log[2] ~ counts ~ per ~ million)) +
		ylab("Density")
		
	print(plot_object)
	
	ggthemr_reset()
	}
}