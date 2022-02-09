source("Nb_calculator.R")
source("targeted_questions.R")
########################################################################################
##### Simple but complex way to get density plots from all conditions into one data.table

targeted_densities <- data.table(merge(targeted_CPM[, unique(condition)], targeted_CPM[, unique(type)]))

setnames(targeted_densities, c("x", "y"), c("condition", "type"))

targeted_densities[, `:=`(x = NA, y = NA)]

for (i in targeted_densities[, unique(condition)]) {
	for (j in targeted_densities[condition == i, unique(`type`)]) {
		
		density__ <- density(targeted_CPM[condition == i & type == j, log2(CPM)], from = -10, to = 15)
		
		density__ <- data.table(condition = i, type = j, x = density__$x, y = density__$y)
		
		targeted_densities <- merge(targeted_densities, density__, all.y = TRUE, all.x = TRUE)
		
	}
}

targeted_densities <- targeted_densities[complete.cases(targeted_densities)]

targeted_densities <- targeted_densities[targeted_CPM[, .(Condition = unique(Condition)), by = .(condition)], on = .(condition)]

targeted_densities_summary <- targeted_densities[
	, .(
		meany = mean(y),
		stdy = std(y),
		.N,
		semy = std(y)/sqrt(.N - 1)
	), 
	by = .(type, Condition, x)]

targeted_densities_summary[is.nan(stdy), `:=`(stdy = 0, semy = 0)]

targeted_densities_summary$CI_lower <- targeted_densities_summary$meany + qt((1-0.85)/2, df=targeted_densities_summary$N-1)*targeted_densities_summary$semy
targeted_densities_summary$CI_upper <- targeted_densities_summary$meany - qt((1-0.85)/2, df=targeted_densities_summary$N-1)*targeted_densities_summary$semy

targeted_densities_summary[is.na(CI_lower), CI_lower := meany]
########################################################################################

for (i in targeted_densities_summary[, unique(Condition)]) {
	for (j in targeted_densities_summary[, unique(type)]) {
		
		ggthemr("flat")
		
		if (j == "control") {plot_shade = "blue"}
		if (j != "control") {plot_shade = "red"}
		
		plot_object <- ggplot(
		targeted_densities_summary[Condition == i & type == j], 
		aes(x = x, y = meany)) +
		geom_line(
			data = targeted_densities[Condition == i & type == j], 
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
		xlim(targeted_densities_summary[, min(x)], targeted_densities_summary[, max(x)]) +
		ylim(targeted_densities_summary[, min(CI_lower, na.rm = T)], targeted_densities_summary[, max(CI_upper, na.rm = T)]) +
		ggtitle(bquote(Log[2] ~ counts ~ per ~ million. ~ bold(.(i)) ~ italic(.(j)) ~ guides. )) +
		xlab(bquote(Log[2] ~ counts ~ per ~ million)) +
		ylab("Density")
		
	print(plot_object)
	
	ggthemr_reset()
	}
}