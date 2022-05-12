source("Publication/publication_counter.R")

#################################################################
#################################################################


Nb_data <- 
	rbind(data_grid_remelted, 
				data_grid_remelted[
					condition %in% c("Inoculum", "Mouse_P1_003"), 
					.(condition = "Baseline", 
						count = median(count)),
					by = .(spacer)])

base_t0 = "Baseline"


#################################################################
#################################################################

botneck_summary <- Nb_data[count > 0, .N,  by = .(condition)]

botneck_summary <- botneck_summary[
	Nb_data[, .(med = median(count)),
										 by = .(condition)], on = .(condition)]

botneck_summary <- botneck_summary[
	Nb_data[spacer %like% "Ctrl",
										 .(control = median(count)),
										 by = .(condition)], on = .(condition)]

botneck_summary <- botneck_summary[
	Nb_data[!(spacer %like% "Ctrl"),
										 .(knockdown = median(count)),
										 by = .(condition)], on = .(condition)]

botneck_summary <- exp_design[
	botneck_summary,
	on = .(condition)]

botneck_summary[!is.na(rep), verbose := paste(
	media, gDNA_source, growth_condition, rep, sep = "_")]

botneck_summary[is.na(rep), verbose := paste(
	media, gDNA_source, growth_condition, sep = "_")]

#################################################################

for (i in botneck_summary[generations > 0, condition]) {

	botneck_calcs <- Nb_data[condition == base_t0, .(spacer, fio = count/sum(count))]
	botneck_calcs <- botneck_calcs[Nb_data[condition == i, .(spacer, fis = count/sum(count))], on = .(spacer)]
	botneck_calcs[, ratio := ((fis - fio)^2) / (fio * (1 - fio)^2)]
	botneck_calcs <- botneck_calcs[!is.na(ratio) & ratio != Inf]

	f_hat <- botneck_calcs[ratio != Inf & !(is.na(ratio)), sum(ratio) ] * ( 1 / length(Nb_data[, unique(spacer)]))
	Nb <- botneck_summary[condition == i, generations] / (f_hat - (1 / Nb_data[condition == base_t0, sum(count)]) - (1 / Nb_data[condition == base_t0, sum(count)]))

	botneck_summary[condition == i, Nb := ..Nb]}

# finding Nb for control_Nbs #######################################
#################################################################

for (i in botneck_summary[generations > 0, condition]) {

	botneck_calcs <- Nb_data[condition == base_t0 & spacer %like% "Ctrl", .(spacer, fio = count/sum(count))]
	botneck_calcs <- botneck_calcs[Nb_data[condition == i & spacer %like% "Ctrl", .(spacer, fis = count/sum(count))], on = .(spacer)]
	botneck_calcs[, ratio := ((fis - fio)^2) / (fio * (1 - fio)^2)]
	botneck_calcs <- botneck_calcs[!is.na(ratio) & ratio != Inf]

	f_hat <- botneck_calcs[ratio != Inf & !(is.na(ratio)), sum(ratio) ] * ( 1 / length(Nb_data[spacer %like% "Ctrl", unique(spacer)]))
	Nb <- botneck_summary[condition == i, generations] / (f_hat - (1 / Nb_data[condition == base_t0 & spacer %like% "Ctrl", sum(count)]) - (1 / Nb_data[condition == base_t0 & spacer %like% "Ctrl", sum(count)]))

	botneck_summary[condition == i, control_Nb := ..Nb]}

# # finding Nb for knockdown_Nb#######################################
# #################################################################

for (i in botneck_summary[generations > 0, condition]) {

	botneck_calcs <- Nb_data[condition == base_t0 & !(spacer %like% "Ctrl"), .(spacer, fio = count/sum(count))]
	botneck_calcs <- botneck_calcs[Nb_data[condition == i & !(spacer %like% "Ctrl"), .(spacer, fis = count/sum(count))], on = .(spacer)]
	botneck_calcs[, ratio := ((fis - fio)^2) / (fio * (1 - fio)^2)]
	botneck_calcs <- botneck_calcs[!is.na(ratio) & ratio != Inf]

	f_hat <- botneck_calcs[ratio != Inf & !(is.na(ratio)), sum(ratio) ] * ( 1 / length(Nb_data[!(spacer %like% "Ctrl"), unique(spacer)]))
	Nb <- botneck_summary[condition == i, generations] / (f_hat - (1 / Nb_data[condition == base_t0 & !(spacer %like% "Ctrl"), sum(count)]) - (1 / Nb_data[condition == base_t0 & !(spacer %like% "Ctrl"), sum(count)]))

	botneck_summary[condition == i, knockdown_Nb := ..Nb]}


to_barplot <- botneck_summary[, .(verbose, Nb, control_Nb, knockdown_Nb)]
to_barplot <- melt(to_barplot, id.vars = "verbose", variable = "group", value.name = "Nb")
to_barplot <- to_barplot[group != "Nb"]
to_barplot <- to_barplot[!is.na(Nb)]
to_barplot[, c("A", "B", "C", "D", "E", "F", "G") := tstrsplit(verbose, "_", type.convert = TRUE, fixed = TRUE)]
to_barplot[,index := 1:.N, by = .(group)]
to_barplot[, Label := paste0(A, "_", index)]

this_plot <- ggplot(
	data = to_barplot, aes(x = verbose, y = Nb, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	scale_fill_brewer(palette = "Paired") +
	ggtitle(paste("Bottleneck Number (Effective Population) by Condition: Controls, Knockdowns.")) +
	# ylim(c(0, 160000)) +
	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))

print(this_plot)
