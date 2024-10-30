dictionary_refseq <- fread("dictionary_refseq.tsv")
orthologs <- fread("orthologs.tsv")

dictionary <- annotated_key
dictionary[locus_tag %like% "RS", locus_tag := paste0("PA14_", locus_tag)]

zz <- dictionary[, .(locus_tag = unique(locus_tag))]

zz <- zz[dictionary_refseq, on = .(locus_tag == `locus tag`)]

zz <- orthologs[zz, on = .(`Locus Tag (Hit)` == PA14_ID)][!is.na(`Gene Name`)]

setnames(zz, c("Gene Name"), c("gene_name"))


pg_etc <- fread("CL930.tsv")

zzz <- zz[missing, on = .(gene_name == V1)]

zzzz <- table(zzz[,`Locus Tag (Query)` %in% pg_etc$identifier])
