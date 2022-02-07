# to make a barplot
this_spacer <- "141_PA1657_Ctrl_2"
promoter_spacer_combo <- dcast(all_counts[spacer == this_spacer][, .(promoter, condition, count)], promoter ~ condition, value.var = "count", fill = 0)
promoter_spacer_combo <- melt(promoter_spacer_combo, id.vars = "promoter", variable.name = "condition", value.name = "count")
ggplot(data = promoter_spacer_combo, aes(x = condition, y = count, fill = promoter)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	geom_text(aes(label = promoter), vjust = 2, color = "white",
						position = position_dodge(0.9), size = 3.5) +
	scale_fill_brewer(palette = "Paired") +
	ggtitle(paste("Promoter distribution for", this_spacer)) +
	theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1))
