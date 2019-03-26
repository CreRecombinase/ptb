library(tidyverse)

dataif <- snakemake@input[["data"]]
zscoref <- snakemake@input[["zscores"]]
p_basef <- snakemake@input[["prior_base"]]
p_fullf <- snakemake@input[["prior_full"]]
data_f <- snakemake@output[["data_f"]]
snp_f <- snakemake@output[["snp_f"]]
zscores <- read_tsv(zscoref)
d <- readRDS(dataif)
d <- d %>% rename(variant = snp)

ff_format <- function(pbf,pff,df,sf){

    prior_base <- read_delim(pbf, delim="  ", col_names=FALSE) %>%
        rename(variant=X1, prior_base = X3) %>% select(-X2)
    prior_full <- read_delim(pff, delim="  ", col_names=FALSE) %>%
        rename(variant=X1, prior_full = X3) %>% select(-X2)

    dat <- filter(zscores, ldchunk==b) %>%
        full_join(., prior_base) %>%
        full_join(., prior_full) %>%
        left_join(., d)
    saveRDS(dat, file=df)
    write_lines(dat$variant, path=sf)
    }
    
}

pwalk(list(
    pbf=p_basef,
    pff=p_fullf,
    df=data_f,
    sf=snp_f),
    ff_format
    )
