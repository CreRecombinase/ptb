library(tidyverse)
library(susieR)

args = commandArgs(trailingOnly=TRUE)
blocks <- as.numeric(args)

options(scipen=999)
#block <- 512
for(r in blocks){
    sres <- readRDS(paste0("susie_finemap/susie.", r, ".RDS"))
    x <- list(fit_bhat_full = summary(sres$fit_bhat_full),
           fit_bhat_base = summary(sres$fit_bhat_base),
           fit_bhat_noprior = summary(sres$fit_bhat_noprior))
    zscores <- readRDS(paste0("zscores2/data.filtered.", r, ".RDS"))
    dat <- zscores %>% select(CHR, variant, position, prior_full, prior_base)
    o <- order(x$fit_bhat_noprior$vars$variable)
    dat$np_pip <- x$fit_bhat_noprior$vars$variable_prob[o]
    dat$np_cs <- x$fit_bhat_noprior$vars$cs[o]
    o <- order(x$fit_bhat_base$vars$variable)
    dat$base_pip <- x$fit_bhat_base$vars$variable_prob[o]
    dat$base_cs <- x$fit_bhat_base$vars$cs[o]
    o <- order(x$fit_bhat_full$vars$variable)
    dat$full_pip <- x$fit_bhat_full$vars$variable_prob[o]
    dat$full_cs <- x$fit_bhat_full$vars$cs[o]
    saveRDS(dat, file=paste0("susie_finemap/susie.summary.", r, ".RDS"))

    for(prior in c("no prior", "base", "full")){
        header <- c(paste0("browser position chr", dat$CHR[1], ":", min(dat$position), "-", max(dat$position)), 
                   paste0("track type=wiggle_0 name=SuSIE description=\"SuSIE PIP ", prior, "\" visibility=full autoScale=off viewLimits=0:1 color=215,17,215 priority=10"), 
                   paste0("variableStep chrom=chr", dat$CHR[1], " span=1"))
        if(prior=="no prior") prf <- "np"
            else prf <- prior
        write_lines(header, path=paste0("susie_finemap/", r,".susie_", prf,  ".wig"))
        wig_dat <-  dat %>% 
               select(position, paste0(prf, "_pip"))  %>% 
               mutate(position =  as.integer(position))
        write.table(wig_dat, file=paste0("susie_finemap/", r,".susie_", prf,  ".wig"), 
                    append=TRUE, col.names=FALSE, row.names=FALSE, sep="\t")
    }
}


