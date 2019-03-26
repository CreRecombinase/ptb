library(tidyverse)

data_file <- snakemake@input[["data"]]
frq_file <- snakemake@input[["frq"]]
dfilter_f <- snakemake@output[["data"]]
sfilter_f <- snakemake@output[["snps"]]
zfilter_f <- snakemake@output[["zscores"]]



dat <- readRDS(data_file)

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

saveRDS(dat, file=dfilter_f)
write_tsv(select(dat, variant, tstat),
          path=zfilter_f, col_names=FALSE)
write_lines(dat$variant, path=sfilter_f)
