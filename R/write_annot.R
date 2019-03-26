library(readr)
library(purrr)
library(dplyr)

#args <- commandArgs(trailingOnly=TRUE)
                                        #out <- args[1]
fl <- snakemake@input[["dsc_beds"]]
chr <- snakemake@params[["chrom"]]
snpf <- snakemake@input[["bim"]]
out <- snakemake@output[["annot"]]



snps <- read_tsv(snpf,
                 col_names=c("CHR", "SNP", "cM", "BP", "A1", "A2"))

annot_df <- map(fl, function(x){read_tsv(x, col_names=c("chr", "start", "stop", "SNP"))})

for(i in seq_along(annots)){
    n <- annots[i]
    snps[[n]] <- with(snps, case_when(SNP %in% annot_df[[i]]$SNP ~ 1, TRUE ~ 0))

}
snps <- select(snps, -cM, -A1, -A2)
write_tsv(snps, path=out)


