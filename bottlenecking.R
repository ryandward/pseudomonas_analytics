
zzz <- data_grid[, .(
	spacer, 
	ratio = ((Mouse_P1_015/sum(Mouse_P1_015) - Mouse_P1_006/sum(Mouse_P1_006)))^2  /  Mouse_P1_006/sum(Mouse_P1_006) * ((1 - Mouse_P1_006/sum(Mouse_P1_006))^2)
)]

nrow(zzz)

f_hat <- zzz[ratio!=Inf, sum(ratio, na.rm = T)]*(1/nrow(zzz))

Nb <- 16/( f_hat - 1/data_grid[, sum(Mouse_P1_006)] - 1/data_grid[, sum(Mouse_P1_015)])

print(Nb)

####################

jan <- fread("jan.tsv")

bottlenecked_levels <- jan[, .(
	ratio = ((`Mice_noDox-9-Lung` / sum(`Mice_noDox-9-Lung`) - Pre1/sum(Pre1))^2) / 
		((Pre1/sum(Pre1)) * ( 1 - Pre1/sum(Pre1))))]

f_hat = bottlenecked_levels[ratio!=Inf & ! (is.na(ratio)), sum(ratio) ] * 
	( 1 / nrow(bottlenecked_levels))

Nb = 42 / ( f_hat - 1 / jan[, sum(Pre1)] - 1/jan[, sum(Pre1)])

#######################

bottlenecked_levels <- data_grid[, .(fis = (Mouse_P1_016 / sum(Mouse_P1_016)), fio = Mouse_P1_003/sum(Mouse_P1_003))]
bottlenecked_levels[, ratio := ((fis-fio)^2) / (fio * (1- fio))]
bottlenecked_levels <- bottlenecked_levels[!is.na(ratio) & ratio!=Inf]

f_hat = bottlenecked_levels[ratio!=Inf & ! (is.na(ratio)), sum(ratio) ] * 
	( 1 / nrow(bottlenecked_levels))

Nb = 10 / ( f_hat - 1 / data_grid[, sum(Mouse_P1_003)] - 1/data_grid[, sum(Mouse_P1_003)])

