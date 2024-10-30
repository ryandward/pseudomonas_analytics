
################################################################################
# Some forensics...
# Use edgeR to determine log-fold changes
source("Publication/publication_counter.R")

spacers.not_in_mating_strains <- anti_join(
	annotated_key, 
	filter(additional, strain == "mfdpir"), 
	by = c("name" = "spacer"))

spacers.in_results <- semi_join(
	filter(additional, strain == "mfdpir"), 
	melted_results, 
	by = c("spacer" = "name")) %>% 
	select(spacer) %>% 
	unique

print(spacers.not_in_mating_strains)

spacers.not_in_results <- anti_join(
	filter(additional, strain == "mfdpir"), 
	melted_results, 
	by = c("spacer" = "name")) %>% 
	select(spacer) %>% 
	unique

genes.in_results <-
	spacers.in_results %>% 
	inner_join(
		annotated_key, 
		by = c("spacer" = "name")) %>% 
	select(gene_name) %>% unique

genes.missing_spacers <-
	spacers.not_in_results %>% 
	inner_join(
		annotated_key, 
		by = c("spacer" = "name")) %>% 
	select(gene_name) %>% unique

genes.not_in_results <-
	anti_join(
		(genes.missing_spacers), 
		genes.in_results)

print(genes.not_in_results)