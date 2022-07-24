#https://www.nature.com/articles/nmeth.3253

source("cleanup_analysis.R")

botneck_stats <- fread("stats.tsv", col.names = c("sequence", "count", "condition"))

botneck_stats <- botneck_stats %>% inner_join(all_guides) %>% filter(type != "focused")

exp_design_botneck <- fread("unified_counts/exp_design.tsv") %>% mutate(sample_group = paste(gDNA_source, growth_condition, media, sep = "_"))


botneck.t0 <- botneck_stats %>%
	inner_join(all_guides) %>% 
	inner_join(exp_design_botneck) %>%
	filter(condition %in% c("dJMP1", "dJMP2")) %>% 
	group_by(sequence, type) %>% 
	summarise(count = median(count)) %>%
	group_by(type) %>%
	mutate(
		fi0 = count/sum(count),
		count0 = count) %>%
	select(type, sequence, fi0, count0) %>%
	nest %>%
	rename(data0 = data) %>%
	mutate(s0 = map_dbl(data0, ~sum(.$count0)))

botneck <- botneck_stats %>%
	inner_join(all_guides) %>% 
	inner_join(exp_design_botneck) %>%
	filter(generations != 0) %>%
	group_by(type, condition, generations) %>%
	mutate(
		fis = count/sum(count)) %>%
	nest %>%
	mutate(
		ss = map_dbl(data, ~sum(.$count))) %>% 
	full_join(botneck.t0) %>% 
	mutate(data = map2(data, data0, inner_join))

botneck <- botneck %>% 
	mutate(
		data = map(
			data, ~.x %>% 
				mutate(
					ratio = ( (fis - fi0)^2) / ( fi0 * (1 - fi0)^2 ) )  ) )

botneck <- botneck %>% 
	mutate(
		f_hat = map_dbl(
			data, 
			~sum(.$ratio)) * (1 / map_dbl(data, ~n_distinct(.$sequence))),
		Nb = generations/(f_hat - 1/s0 - 1/ss))

botneck.stats <- botneck %>% inner_join(exp_design_botneck) %>%
	group_by(sample_group, type) %>% 
	summarise(
		Nb.med = median(Nb), 
		Nb.range = max(Nb) - min(Nb),
		Nb.mean = mean(Nb),
		Nb.sd = sd(Nb))

botneck.stats <- botneck.stats %>% mutate(Nb.sd = case_when(is.na(Nb.sd) ~ 0, TRUE ~ Nb.sd))

botneck.stats %>% inner_join(exp_design_botneck) %>% 
	filter(type != "focused") %>% 
	filter(!sample_group %like% "100x") %>%
	filter(!media %like% "Gent") %>%
	filter(gDNA_source == "plated") %>%
	mutate(sample_group = paste(media)) %>%
	ggplot(aes(
		x = sample_group, 
		y = Nb.mean, 
		fill = type)) +
	geom_bar(
		stat = "identity", 
		position = position_dodge(),
		colour = "black") +
	geom_errorbar(aes(
		x = sample_group,
		y = Nb.mean,
		ymin = Nb.mean - Nb.sd,
		ymax = Nb.mean + Nb.sd),
		width = 0.2,
		position = position_dodge(0.95)) +
	scale_fill_viridis(discrete = T, option = 'mako', alpha = 0.25) +
	theme_bw() +
	# theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 1)) +
	ggtitle("Population Complexity: Bottleneck Metric") +
	scale_y_continuous(
		trans = scales::pseudo_log_trans(base = 10),
		breaks = c(0, 10^seq(1,7)),
		labels = label_number(scale_cut = cut_short_scale())) +
	theme(axis.title.x = element_blank()) 

