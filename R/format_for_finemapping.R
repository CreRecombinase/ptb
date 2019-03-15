library(tidyverse)

#args = commandArgs(trailingOnly=TRUE)
#blocks <- as.numeric(args)
#b <- 512
ff_format <- function(blocks){
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
        saveRDS(dat, file=paste0("data.", b, ".RDS"))
        write_lines(dat$variant, path=paste0("snps.", b, ".txt"))
    }
    
}

# plink --bfile /project2/xinhe/1kg/plink_format/EUR.1kg \ 
# --extract snps.15.txt \
# --freq --out 15


# ~/dap/dap_src/dap-g -d_z zscores.filtered.15.tsv \
#                     -d_ld 15.filtered.ld \
#                     -o dapfinemap.15.result \
#                     --output_all -p torus_base1/full/15.prior

strand_flip <- function(b){
    dat <- readRDS(file=paste0("data.", b, ".RDS"))
    frq <- read_table(paste0(b, ".frq"))
    if(!nrow(dat) == nrow(frq)) cat("Warning: Some variants are being removed because they are not present in 1k genomes\n")
    frq <- frq %>% rename(variant = SNP) 
    dat <- right_join(frq, dat, by="variant")
    #Strand switch if necessary
    strand <- which(!(dat$ref_allele == dat$A1 | dat$ref_allele == dat$A2))
    str_flp <- c("A" = "T", "G" = "C", "C" = "G", "T" = "A")
    dat$ref_allele[strand] <- str_flp[dat$ref_allele[strand]]
    dat$alt_allele[strand] <- str_flp[dat$alt_allele[strand]]
    #Remove SNPs that still don't match
    remove <- which(!(dat$alt_allele == dat$A1 | dat$alt_allele == dat$A2))
    dat <- dat[-remove,]
    #Flip effects where necessary
    flip <- which(dat$alt_allele == dat$A1)
    dat_copy <- dat
    dat$ref_allele[flip] <- dat_copy$alt_allele[flip]
    dat$alt_allele[flip] <- dat_copy$ref_allele[flip]
    dat$effect[flip] <- -1*dat_copy$effect[flip]

    saveRDS(dat, file=paste0("data.filtered.", b, ".RDS"))
    write_tsv(select(dat, variant, tstat), path=paste0("zscores.filtered.", b, ".tsv"), col_names=FALSE)
    write_lines(dat$variant, path=paste0("snps.filtered.", b, ".txt"))
}

# plink --bfile /project2/xinhe/1kg/plink_format/EUR.1kg \ 
# --extract snps.filtered.15.txt --r square \ 
# --freq --out 15.filtered

run_susie <- function(b, L=10){
    dat <- readRDS(paste0("data.filtered.", b, ".RDS"))
    R <- read_tsv(paste0(b, ".filtered.ld"), col_names=FALSE)
    R <- as.matrix(R)
    n <- 43568
    fit_bhat_full <- with(dat,
                      susie_bhat(bhat = effect, shat = stderr, 
                      R = R, n = n, 
                      r_tol = 1e-4,
                      scaled_prior_variance=0.1, 
                      estimate_residual_variance=F,
                      prior_weights = prior_full, 
                      min_abs_corr = 0.1, 
                      L = L))

    fit_bhat_base <- with(dat,
                      susie_bhat(bhat = effect, shat = stderr, 
                      R = R, n = n, 
                      r_tol = 1e-4,
                      scaled_prior_variance=0.1, 
                      estimate_residual_variance=F,
                      prior_weights = prior_base, 
                      min_abs_corr = 0.1, 
                      L = L))

    fit_bhat_noprior <- with(dat,
                      susie_bhat(bhat = effect, shat = stderr, 
                      R = R, n = n, 
                      r_tol = 1e-4,
                      scaled_prior_variance=0.1, 
                      estimate_residual_variance=F,
                      min_abs_corr = 0.1, 
                      L = L)) 
    ret <- list(fit_bhat_full = fit_bhat_full, fit_bhat_base = fit_bhat_base, fit_bhat_noprior = fit_bhat_noprior)
    saveRDS(dat, file=paste0("susie.", b, ".RDS"))

}
