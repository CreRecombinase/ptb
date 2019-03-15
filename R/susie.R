library(tidyverse)
library(susieR)

args = commandArgs(trailingOnly=TRUE)
b <- args[1]
data_file <- args[2]
ld_file <- args[3]
dir <- args[4]

run_susie <- function(b, data_file, ld_file, dir,  L=10){
    #jdat <- readRDS(paste0("data.filtered.", b, ".RDS"))
    dat <- readRDS(data_file)
    #R <- read_tsv(paste0(b, ".filtered.ld"), col_names=FALSE)
    R <- read_tsv(ld_file, col_names=FALSE)
    R <- as.matrix(R)
    n <- 43568
    fit_bhat_full <- with(dat,
                      susie_bhat(bhat = effect, shat = stderr, 
                      R = R, n = n, 
                      r_tol = 1e-3,
                      scaled_prior_variance=0.1, 
                      estimate_residual_variance=F,
                      prior_weights = prior_full, 
                      min_abs_corr = 0.1, 
                      L = L))

    fit_bhat_base <- with(dat,
                      susie_bhat(bhat = effect, shat = stderr, 
                      R = R, n = n, 
                      r_tol = 1e-3,
                      scaled_prior_variance=0.1, 
                      estimate_residual_variance=F,
                      prior_weights = prior_base, 
                      min_abs_corr = 0.1, 
                      L = L))

    fit_bhat_noprior <- with(dat,
                      susie_bhat(bhat = effect, shat = stderr, 
                      R = R, n = n, 
                      r_tol = 1e-3,
                      scaled_prior_variance=0.1, 
                      estimate_residual_variance=F,
                      min_abs_corr = 0.1, 
                      L = L)) 
    ret <- list(fit_bhat_full = fit_bhat_full, fit_bhat_base = fit_bhat_base, fit_bhat_noprior = fit_bhat_noprior)
    saveRDS(ret, file=paste0(dir, "susie.", b, ".RDS"))

}

run_susie(b, data_file, ld_file, dir,  L=10)
