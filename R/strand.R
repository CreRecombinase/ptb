library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

b <- args[1]
data_file <- args[2]
frq_file <- args[3]
out_dir <- args[4]

strand_flip <- function(b, data_file, frq_file, out_dir){
    #dat <- readRDS(file=paste0("data.", b, ".RDS"))
    dat <- readRDS(data_file)
    #frq <- read_table(paste0(b, ".frq"))
    frq <- read_table(frq_file)

    frq <- frq %>% rename(variant = SNP) 
    dat <- right_join(frq, dat, by="variant")
    #Strand switch if necessary
    strand <- which(!(dat$ref_allele == dat$A1 | dat$ref_allele == dat$A2))
    str_flp <- c("A" = "T", "G" = "C", "C" = "G", "T" = "A")
    dat$ref_allele[strand] <- str_flp[dat$ref_allele[strand]]
    dat$alt_allele[strand] <- str_flp[dat$alt_allele[strand]]
    #Remove SNPs that still don't match
    dat <- filter(dat, (ref_allele==A1 & alt_allele == A2) | (alt_allele == A1 & ref_allele ==A2))
    #Flip effects where necessary
    flip <- which(dat$alt_allele == dat$A1)
    dat_copy <- dat
    dat$ref_allele[flip] <- dat_copy$alt_allele[flip]
    dat$alt_allele[flip] <- dat_copy$ref_allele[flip]
    dat$effect[flip] <- -1*dat_copy$effect[flip]

    saveRDS(dat, file=paste0(out_dir, "data.filtered.", b, ".RDS"))
    write_tsv(select(dat, variant, tstat), 
              path=paste0(out_dir, "zscores.filtered.", b, ".tsv"), col_names=FALSE)
    write_lines(dat$variant, path=paste0(out_dir, "snps.filtered.", b, ".txt"))
}

strand_flip(b, data_file, frq_file, out_dir)
