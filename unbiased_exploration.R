require('pacman'); 
p_load(data.table, edgeR, reshape2, digest, forcats, pheatmap)

clipboard <- function(
	x , 
	sep = "\t" , 
	row.names = FALSE , 
	col.names = TRUE 
) {
	con <- pipe(
		"xclip -selection clipboard -i" , 
		open = "w"
	)
	
	write.table(
		x , 
		con , 
		sep = sep , 
		row.names = row.names , 
		col.names = col.names
	)
	close(con)}

RowVar <- function(x , ...) {
	rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

all_guides <- fread("awk 'NR%2==1' /home/ryandward/Mobile_CRISPRi_Neha_Ryan/GuidesOnly.fasta",
						col.names = c("guide"))

all_guides[, guide := gsub(">", "", guide)]

guide_matches <- fread("/home/ryandward/Mobile_CRISPRi_Neha_Ryan/UnbiasedP1lib_invivo_sequencing_reads/guide_stats/guide_matches.tsv",
											 col.names = c(
											 	"ID",
											 	"guide",
											 	"condition"))

guide_matches[, paste0("ID", 1:2) := tstrsplit(ID, " ", type.convert = TRUE, fixed = TRUE)]
guide_matches[, ID := NULL]
guide_matches[, ID2 := NULL]

setindex(guide_matches, condition, ID1)

counts <- guide_matches[
	, 
	.(count = .N),
	by = .(guide, 
				 condition)]

rm(guide_matches)

# Simplify names of groups
counts[, paste0("ID", 1:3) := tstrsplit(condition, "_", type.convert = TRUE, fixed = TRUE)]


counts[, run := as.numeric(run)]

###############
###############
###############
require('pacman'); 
p_load(data.table, edgeR, reshape2, digest, forcats, pheatmap)

#######################################################

setup <- fread("/home/ryandward/Mobile_CRISPRi_Neha_Ryan/setup.tsv")

counts <- setup[counts, on = .(condition)]

setnames(counts, c("guide", "condition", "count", "ID", "promoter", "rep"))

counts[is.na(promoter), promoter := "Unknown"]
counts[is.na(rep), rep := NA]


unbiased_counts <- data.table::dcast(counts[ID != "Undetermined" & guide %like% "Ctrl"] , guide ~ ID + promoter + rep, value.var = "count", fill = 0)

plot_matrix <- data.matrix(log2(unbiased_counts[, -1] + 1))
# plot_matrix <- data.matrix((unbiased_counts[, -1]))

rownames(plot_matrix) <- unbiased_counts$guide

break_halves <- length(unique(as.vector(plot_matrix)))

breaks <- c(
	seq(min(plot_matrix), 
			mean(plot_matrix),
			length.out = break_halves)[-break_halves], 
	seq(mean(plot_matrix), 
			max(plot_matrix), 
			length.out = break_halves))

plot_colors <- c(
	colorRampPalette(c("#ba000d", "white"))(break_halves)[-break_halves], 
	colorRampPalette(c("white", "#007ac1"))(break_halves)[-1])

to_plot_title <- paste("Log2 Unbiased Count + 1 (Only Controls)")

to_plot <- pheatmap(plot_matrix,
										col = plot_colors,
										breaks = breaks,
										border_color = NA,
										# cutree_rows = 4,
										# cutree_cols = 4,
										main = to_plot_title,
										angle_col = 315,
										# fontsize_col = 10,
										fontsize_row = 10,
										cluster_cols = FALSE,

										# cutree_cols = 10,
										clustering_method = "ward.D2",
										clustering_distance_rows = "maximum",
										clustering_distance_cols = "maximum"
)


