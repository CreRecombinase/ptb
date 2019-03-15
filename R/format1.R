library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
dir <- args[1]
args <- args[-1]
ff_format <- function(blocks, dir){
    zscores <- read_tsv("ga.zscores.tsv.gz")
    d <- readRDS("/project/mstephens/data/external_public_supp/gwas_summary_statistics/raw_summary_statistics/23_and_me_ptb/23andMe_ga.RDS")
    d <- d %>% rename(variant = snp)
    #vars_1kg <- read_tsv("/project2/xinhe/1kg/plink_format/EUR.1kg.bim",
    #                     col_names=c("CHR", "variant", "cM", "pos", "A1", "A2")) %>%

    for(b in blocks){
        prior_base <- read_delim(paste0("torus_base1/base/", b, ".prior"), delim="  ", col_names=FALSE) %>%    
                    rename(variant=X1, prior_base = X3) %>% select(-X2)
        prior_full <- read_delim(paste0("torus_base1/full/", b, ".prior"), delim="  ", col_names=FALSE) %>%    
                    rename(variant=X1, prior_full = X3) %>% select(-X2)

        dat <- filter(zscores, ldchunk==b) %>%
               full_join(., prior_base) %>%
               full_join(., prior_full) %>%
               left_join(., d)
        saveRDS(dat, file=paste0(dir, "data.", b, ".RDS"))
        write_lines(dat$variant, path=paste0(dir, "snps.", b, ".txt"))
    }
    
}
ff_format(args, dir)
