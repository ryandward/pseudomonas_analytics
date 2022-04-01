source("deep_integration_analysis.R")

#################################################################
#################################################################

base_t0 = "Mouse_P1_003"

# #################################################################
# #################################################################
# 
# botneck_summary <- data_grid_remelted[count > 0, .N,  by = .(condition)]
# 
# botneck_summary <- botneck_summary[
# 	data_grid_remelted[, .(med = median(count)), 
# 										 by = .(condition)], on = .(condition)]
# 
# botneck_summary <- botneck_summary[
# 	data_grid_remelted[spacer %like% "Ctrl", 
# 										 .(control = median(count)), 
# 										 by = .(condition)], on = .(condition)]
# 
# botneck_summary <- botneck_summary[
# 	data_grid_remelted[!(spacer %like% "Ctrl"), 
# 										 .(knockdown = median(count)), 
# 										 by = .(condition)], on = .(condition)]
# 
# botneck_summary <- exp_design[
# 	botneck_summary, 
# 	on = .(condition)]
# 
# botneck_summary[!is.na(rep), verbose := paste(
# 	media, gDNA_source, growth_condition, rep, sep = "_")]
# 
# botneck_summary[is.na(rep), verbose := paste(
# 	media, gDNA_source, growth_condition, sep = "_")]
# 
# #################################################################
# 
# for (i in botneck_summary[generations > 0, condition]) {
# 	
# 	botneck_calcs <- data_grid_remelted[condition == base_t0, .(spacer, fio = count/sum(count))]
# 	botneck_calcs <- botneck_calcs[data_grid_remelted[condition == i, .(spacer, fis = count/sum(count))], on = .(spacer)]
# 	botneck_calcs[, ratio := ((fis - fio)^2) / (fio * (1 - fio))]
# 	botneck_calcs <- botneck_calcs[!is.na(ratio) & ratio != Inf]
# 	
# 	f_hat <- botneck_calcs[ratio != Inf & !(is.na(ratio)), sum(ratio) ] * ( 1 / length(data_grid_remelted[, unique(spacer)]))
# 	Nb <- botneck_summary[condition == i, generations] / (f_hat - (1 / data_grid_remelted[condition == base_t0, sum(count)]) - (1 / data_grid_remelted[condition == base_t0, sum(count)]))
# 	
# 	botneck_summary[condition == i, Nb := ..Nb]}
# 
# # finding Nb for control_Nbs #######################################
# #################################################################
# 
# for (i in botneck_summary[generations > 0, condition]) {
# 	
# 	botneck_calcs <- data_grid_remelted[condition == base_t0 & spacer %like% "Ctrl", .(spacer, fio = count/sum(count))]
# 	botneck_calcs <- botneck_calcs[data_grid_remelted[condition == i & spacer %like% "Ctrl", .(spacer, fis = count/sum(count))], on = .(spacer)]
# 	botneck_calcs[, ratio := ((fis - fio)^2) / (fio * (1 - fio))]
# 	botneck_calcs <- botneck_calcs[!is.na(ratio) & ratio != Inf]
# 	
# 	f_hat <- botneck_calcs[ratio != Inf & !(is.na(ratio)), sum(ratio) ] * ( 1 / length(data_grid_remelted[spacer %like% "Ctrl", unique(spacer)]))
# 	Nb <- botneck_summary[condition == i, generations] / (f_hat - (1 / data_grid_remelted[condition == base_t0 & spacer %like% "Ctrl", sum(count)]) - (1 / data_grid_remelted[condition == base_t0 & spacer %like% "Ctrl", sum(count)]))
# 	
# 	botneck_summary[condition == i, control_Nb := ..Nb]}
# 
# # finding Nb for knockdown_Nb#######################################
# #################################################################
# 
# for (i in botneck_summary[generations > 0, condition]) {
# 	
# 	botneck_calcs <- data_grid_remelted[condition == base_t0 & !(spacer %like% "Ctrl"), .(spacer, fio = count/sum(count))]
# 	botneck_calcs <- botneck_calcs[data_grid_remelted[condition == i & !(spacer %like% "Ctrl"), .(spacer, fis = count/sum(count))], on = .(spacer)]
# 	botneck_calcs[, ratio := ((fis - fio)^2) / (fio * (1 - fio))]
# 	botneck_calcs <- botneck_calcs[!is.na(ratio) & ratio != Inf]
# 	
# 	f_hat <- botneck_calcs[ratio != Inf & !(is.na(ratio)), sum(ratio) ] * ( 1 / length(data_grid_remelted[!(spacer %like% "Ctrl"), unique(spacer)]))
# 	Nb <- botneck_summary[condition == i, generations] / (f_hat - (1 / data_grid_remelted[condition == base_t0 & !(spacer %like% "Ctrl"), sum(count)]) - (1 / data_grid_remelted[condition == base_t0 & !(spacer %like% "Ctrl"), sum(count)]))
# 	
# 	botneck_summary[condition == i, knockdown_Nb := ..Nb]}
# 
# # Barplots                 ######################################
# #################################################################
# 
# to_barplot <- botneck_summary[, .(verbose, med, control, knockdown)]
# to_barplot <- melt(
# 	to_barplot, 
# 	id.vars = "verbose", 
# 	variable = "group", 
# 	value.name = "med")
# 
# to_barplot <- to_barplot[group != "med"]
# 
# this_plot <- ggplot(
# 	data = to_barplot, 
# 	aes(x = verbose, y = med, fill = group)) +
# 	geom_bar(stat = "identity", position = position_dodge()) +
# 	scale_fill_brewer(palette = "Paired") +
# 	ggtitle(paste("Median Count by Condition: Controls, Knockdowns.")) +
# 	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))
# 
# ggthemr("flat")
# print(this_plot)
# 
# #################################################################
# 
# to_barplot <- botneck_summary[, .(verbose, med, control, knockdown)]
# to_barplot <- melt(to_barplot, id.vars = "verbose", variable = "group", value.name = "med")
# to_barplot <- to_barplot[group != "med"]
# 
# to_barplot[
# 	, med := to_barplot[
# 		i  = group == "control",
# 		j  = .(ctrl = med),
# 		by = .(verbose)][
# 			i  = .SD, 
# 			on = .(verbose),
# 			j  = .(med = med / ctrl),
# 			by = .EACHI]$med]
# 
# this_plot <- ggplot(
# 	data = to_barplot, aes(x = verbose, y = med, fill = group)) +
# 	geom_bar(stat = "identity", position = position_dodge()) +
# 	scale_fill_brewer(palette = "Paired") +
# 	ggtitle(paste("Relative Raw Counts (Knockdown/Control)")) +
# 	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))
# 
# ggthemr("dust")
# print(this_plot)
# 
# #################################################################
# 
# to_barplot <- botneck_summary[, .(verbose, Nb, control_Nb, knockdown_Nb)]
# to_barplot <- melt(to_barplot, id.vars = "verbose", variable = "group", value.name = "Nb")
# to_barplot <- to_barplot[group != "Nb"]
# to_barplot <- to_barplot[!is.na(Nb)]
# to_barplot[, c("A", "B", "C", "D", "E", "F", "G") := tstrsplit(verbose, "_", type.convert = TRUE, fixed = TRUE)]
# to_barplot[,index := 1:.N, by = .(group)]
# to_barplot[, Label := paste0(A, "_", index)]
# 
# this_plot <- ggplot(
# 	data = to_barplot, aes(x = Label, y = Nb, fill = group)) +
# 	geom_bar(stat = "identity", position = position_dodge()) +
# 	scale_fill_brewer(palette = "Paired") +
# 	ggtitle(paste("Bottleneck Number (Effective Population) by Condition: Controls, Knockdowns.")) +
# 	ylim(c(0, 160000)) +
# 	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))
# 
# ggthemr("flat")
# print(this_plot)
# 
# #################################################################
# #################################################################
# #################################################################
# rm(botneck_summary)
# Same operations on Counts per Million (CPM)
#################################################################
#################################################################
#################################################################
CPM_melted <- melt(data.table(data_CPM, keep.rownames = "spacer"), id.vars = "spacer", variable.name = "condition", value.name = "CPM")

botneck_summary <- CPM_melted[CPM > 0, .N,  by = .(condition)]

botneck_summary <- botneck_summary[
	CPM_melted[, .(med = median(CPM)), 
						 by = .(condition)], on = .(condition)]

botneck_summary <- botneck_summary[
	CPM_melted[spacer %like% "Ctrl", 
						 .(control = median(CPM)), 
						 by = .(condition)], on = .(condition)]

botneck_summary <- botneck_summary[
	data_grid_remelted[!(spacer %like% "Ctrl"), 
										 .(knockdown = median(count)), 
										 by = .(condition)], on = .(condition)]

botneck_summary <- exp_design[botneck_summary, on = .(condition)]
botneck_summary[!is.na(rep), verbose := paste(media, gDNA_source, growth_condition, rep, sep = "_")]
botneck_summary[is.na(rep), verbose := paste(media, gDNA_source, growth_condition, sep = "_")]

#################################################################
for (i in botneck_summary[generations > 0, condition]) {
	
	botneck_calcs <- data_grid_remelted[condition == base_t0, .(spacer, fio = count/sum(count))]
	
	botneck_calcs <- botneck_calcs[data_grid_remelted[condition == i, .(spacer, fis = count/sum(count))], on = .(spacer)]
	botneck_calcs[, ratio := ((fis - fio)^2) / (fio * (1 - fio))]
	botneck_calcs <- botneck_calcs[!is.na(ratio) & ratio != Inf]
	
	f_hat <- botneck_calcs[ratio != Inf & !(is.na(ratio)), sum(ratio) ] * ( 1 / length(data_grid_remelted[, unique(spacer)]))
	Nb <- botneck_summary[condition == i, generations] / (f_hat - (1 / data_grid_remelted[condition == base_t0, sum(count)]) - (1 / data_grid_remelted[condition == base_t0, sum(count)]))
	
	botneck_summary[condition == i, Nb := ..Nb]}
# finding Nb for control_Nbs #######################################
#################################################################
for (i in botneck_summary[generations > 0, condition]) {
	
	botneck_calcs <- data_grid_remelted[condition == base_t0 & spacer %like% "Ctrl", .(spacer, fio = count/sum(count))]
	
	botneck_calcs <- botneck_calcs[data_grid_remelted[condition == i & spacer %like% "Ctrl", .(spacer, fis = count/sum(count))], on = .(spacer)]
	botneck_calcs[, ratio := ((fis - fio)^2) / (fio * (1 - fio))]
	botneck_calcs <- botneck_calcs[!is.na(ratio) & ratio != Inf]
	
	f_hat <- botneck_calcs[ratio != Inf & !(is.na(ratio)), sum(ratio) ] * ( 1 / length(data_grid_remelted[spacer %like% "Ctrl", unique(spacer)]))
	Nb <- botneck_summary[condition == i, generations] / (f_hat - (1 / data_grid_remelted[condition == base_t0 & spacer %like% "Ctrl", sum(count)]) - (1 / data_grid_remelted[condition == base_t0 & spacer %like% "Ctrl", sum(count)]))
	
	botneck_summary[condition == i, control_Nb := ..Nb]}
# finding Nb for knockdown_Nb#######################################
#################################################################
for (i in botneck_summary[generations > 0, condition]) {
	
	botneck_calcs <- data_grid_remelted[condition == base_t0 & !(spacer %like% "Ctrl"), .(spacer, fio = count/sum(count))]
	
	botneck_calcs <- botneck_calcs[data_grid_remelted[condition == i & !(spacer %like% "Ctrl"), .(spacer, fis = count/sum(count))], on = .(spacer)]
	botneck_calcs[, ratio := ((fis - fio)^2) / (fio * (1 - fio))]
	botneck_calcs <- botneck_calcs[!is.na(ratio) & ratio != Inf]
	
	f_hat <- botneck_calcs[ratio != Inf & !(is.na(ratio)), sum(ratio) ] * ( 1 / length(data_grid_remelted[!(spacer %like% "Ctrl"), unique(spacer)]))
	Nb <- botneck_summary[condition == i, generations] / (f_hat - (1 / data_grid_remelted[condition == base_t0 & !(spacer %like% "Ctrl"), sum(count)]) - (1 / data_grid_remelted[condition == base_t0 & !(spacer %like% "Ctrl"), sum(count)]))
	
	botneck_summary[condition == i, knockdown_Nb := ..Nb]}

# Barplots                 ######################################
#################################################################

to_barplot <- botneck_summary[, .(verbose, med, control, knockdown)]
to_barplot <- melt(to_barplot, id.vars = "verbose", variable = "group", value.name = "med")
to_barplot <- to_barplot[group != "med"]

this_plot <- ggplot(data = to_barplot, aes(x = verbose, y = med, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	ggtitle(paste("Median Count by Condition (CPM): Controls, Knockdowns.")) +
	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))

ggthemr("dust")
print(this_plot)

#################################################################

to_barplot <- botneck_summary[, .(verbose, med, control, knockdown)]
to_barplot <- melt(to_barplot, id.vars = "verbose", variable = "group", value.name = "med")
to_barplot <- to_barplot[group != "med"]

to_barplot[
	, med := to_barplot[
		i  = group == "control",
		j  = .(ctrl = med),
		by = .(verbose)][
			i  = .SD, 
			on = .(verbose),
			j  = .(med = med / ctrl),
			by = .EACHI]$med]

this_plot <- ggplot(data = to_barplot, aes(x = verbose, y = med, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	ggtitle(paste("Relative CPM (Knockdown/Control)")) +
	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))

ggthemr("dust")
print(this_plot)

#################################################################

to_barplot <- botneck_summary[, .(condition, med, control, knockdown)]
to_barplot <- melt(to_barplot, id.vars = "condition", variable = "group", value.name = "med")
to_barplot <- to_barplot[group != "med"]
to_barplot <- exp_design[to_barplot, on = .(condition)]
to_barplot[, level := paste(media, gDNA_source, growth_condition)]
to_barplot[, rep := as.character(rep)]
to_barplot$rep <- factor(to_barplot$rep, levels = c(1:10))

to_barplot[
	, med := to_barplot[
		i  = group == "control",
		j  = .(ctrl = med),
		by = .(verbose)][
			i  = .SD, 
			on = .(verbose),
			j  = .(med = med / ctrl),
			by = .EACHI]$med]

this_plot <- ggplot(data = to_barplot[group == "knockdown"], aes(x = level, y = med, fill = rep)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	ggtitle(paste("Relative CPM (Knockdown/Control)")) +
	scale_fill_brewer(palette = "Paired") +
	
	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))

ggthemr_reset()
print(this_plot)


#################################################################

to_barplot <- botneck_summary[, .(verbose, Nb, control_Nb, knockdown_Nb)]
to_barplot <- melt(to_barplot, id.vars = "verbose", variable = "group", value.name = "Nb")
to_barplot <- to_barplot[group != "Nb"]

this_plot <- ggplot(data = to_barplot, aes(x = verbose, y = Nb, fill = group)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	ggtitle(paste("Bottleneck Number (Effective Population, CPM) by Condition: Controls, Knockdowns.")) +
	theme(axis.text.x = element_text(angle = 55, vjust = 1.0, hjust = 1))

ggthemr("flat")
print(this_plot)

# save these as 1500 x 750

