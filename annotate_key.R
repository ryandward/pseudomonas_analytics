require('pacman')
p_load(data.table, scales, edgeR, statmod, poolr, pheatmap, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, RColorBrewer)




exp_design <- fread("exp_design.tsv", 
										na.strings = c("NA"))

exp_design <- exp_design[condition != "Undetermined"]

################################################################################
# Fix the naming conventions

guide_key <- fread("guide_sequences.tsv", 
									 header = FALSE,
									 col.names = c("name", 
									 							"sequence"))
guide_key[, paste0("ID", 1:4) := tstrsplit(name, "_", type.convert = TRUE, fixed = TRUE)]

NC_002516.2 <- fread("NC_002516.2.bed",
										 header = FALSE,
										 col.names = c(
										 	"chromosome",
										 	"left",
										 	"right",
										 	"locus_tag",
										 	"gene_name",
										 	"strand",
										 	"feature",
										 	"completeness"
										 ))

NC_002516.2 <- NC_002516.2[feature %like% "gene"]
NC_002516.2[, strain := "PAO1"]

annotated_key <- merge(guide_key,
											 header = FALSE,
											 NC_002516.2, 
											 by.x = "ID2", 
											 by.y = "locus_tag", 
											 all.x = FALSE, 
											 all.y = FALSE)

NC_008463.1 <- fread("NC_008463.1.bed",
										 col.names = c(
										 	"chromosome",
										 	"left",
										 	"right",
										 	"locus_tag",
										 	"gene_name",
										 	"strand",
										 	"feature",
										 	"completeness"
										 ))
NC_008463.1 <- NC_008463.1[feature %like% "gene"]
NC_008463.1[, c("strain", "locus_tag") := tstrsplit(locus_tag, "_", type.convert = TRUE, fixed = TRUE)]
NC_008463.1[, strain := "PA14"]

annotated_key <- rbind(annotated_key, 
											 merge(guide_key, 
											 			NC_008463.1, 
											 			by.x = "ID2", 
											 			by.y = "locus_tag", 
											 			all.x = FALSE, 
											 			all.y = FALSE))

lost_guides <- guide_key[!(name %in% annotated_key$name)]

# fix a the "unk" PA3145 that seem to have been messed up... missing a column
lost_guides[ID2 == "unk", `:=` (ID2 = ID1, ID3 = NA, ID4 = ID3)]

# put it back into the annotated guides list
annotated_key <- rbind(annotated_key, 
											 merge(lost_guides, 
											 			NC_002516.2, 
											 			by.x = "ID2", 
											 			by.y = "locus_tag", 
											 			all.x = FALSE, 
											 			all.y = FALSE))

# one gene RS22570 that is still lost
lost_guides <- guide_key[!(name %in% annotated_key$name)]
lost_guides[ID3 == "-", `:=` (ID3 = "hypothetical")]

# add the appropriate columns to just shove it onto the list
lost_guides <- merge(lost_guides, 
										 NC_002516.2, 
										 by.x = "ID2", 
										 by.y = "locus_tag", 
										 all.x = TRUE, 
										 all.y = FALSE)

# just shove it onto the list
annotated_key <- rbind(annotated_key, lost_guides)

# fix some of the strain names
annotated_key[name %like% "PA14_", strain := "PA14" ]
annotated_key[!name %like% "PA14_", strain := "PAO1" ]
annotated_key[ID2 %like% "RS" & is.na(strain),  `:=` (strain = "PA14", chromosome = "NC_008463.1") ]

if(nrow(annotated_key[ID1 != strain])==0)
{
	annotated_key[, strain := ID1]
}

# guide_key[!(name %in% annotated_key$name)] should now be zero

unfound_guide_message <- paste("The number of guides with information not found:",nrow(guide_key[!(name %in% annotated_key$name)]), ".")
print(paste(unfound_guide_message, "If this number is not 0, something bad happened. Look carefully at the input."))

annotated_key[, type := "unknown"]

annotated_key[gene_name == "." & !(name %like% "Ctrl"), gene_name := ID3]
annotated_key[gene_name == "-", gene_name := "."]

annotated_key[, locus_tag := ID2]
annotated_key[, offset := ID4]

annotated_key[name %like% "Ctrl", `:=` (type = "control", locus_tag = "control", gene_name = "control")]

annotated_key <- annotated_key[, .(chromosome,
																	 left,
																	 right,
																	 locus_tag,
																	 gene_name,
																	 strand,
																	 feature,
																	 completeness,
																	 name,
																	 offset,
																	 type)]
#/fix
################################################################################