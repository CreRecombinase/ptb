library(readr)
args <- commandArgs(trailingOnly=TRUE)
x <- as.numeric(args[1])

full_annot <- read_tsv("base2/full.annot.tsv.gz")

base_ix <- 2:50
nm <- names(full_annot)[50+x]
nm <- stringr::str_replace(nm, "_d$", "")
nm <- stringr::str_replace(nm, "_c$", "")
an <- full_annot[c(1:50, 50+x)]
write_tsv(an, path=paste0("base2/", nm, ".oneplus.annot.tsv.gz"))
