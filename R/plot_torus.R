library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

#args = commandArgs(trailingOnly=TRUE)
#prefix <- args[1]
#file_ext <- args[2]
#ldscannot <- as.logical(args[3])


dsc <- readRDS("/project2/xinhe/jean/ptb/strat_ldscore/dsc_annotations/dsc_annot.RDS")
name_dict <- read_csv("/project2/mstephens/ptb/annot_name_dict.txt", col_names=c("annot", "name")) %>%
                filter(annot %in% dsc)
name_dict <- mutate(name_dict, file= paste0("ga.", annot, ".oneplus.rst"))
name_dict$name <- str_replace(name_dict$name, fixed(" (DSC)"), "")
### One plus enrichments


se_from_ci <- function(ci_width, max_se = 1000){
    if(is.na(ci_width)) return(NA)
    f <- function(s, ci_width){
        w <- qnorm(0.975, sd=s)-qnorm(0.025, sd=s)
        abs(w-ci_width)
    }
    opt <- optimize(f, interval=c(0, max_se), ci_width=ci_width)
    return(opt$minimum)
}
res <- apply(name_dict, 1, function(x){
              f <- x[3]
              nm <- x[1]
              rst <- read_lines(f, n_max=60)
              i <- grep(nm, rst)
              rst <- rst[i]
              rst <- unlist(str_split(rst, " "))
              rst <- rst[str_length(rst) > 0]
              est <- as.numeric(rst[2])
              ci_lower <- as.numeric(rst[3])
              ci_upper <- as.numeric(rst[4])
              se <- se_from_ci(ci_upper-ci_lower)
              return(c(est, ci_lower, ci_upper, se))
})
res <- data.frame(t(res))
names(res) <- c("Est", "CI_Lower", "CI_Upper", "SE")
res <- cbind(name_dict, res)

thresh <- 0.05/nrow(res)

res <- res %>% mutate(CI_Lower = Est -qnorm(1-thresh/2, sd=SE), 
              CI_Upper = Est +qnorm(1-thresh/2, sd=SE), 
              sig = abs(Est/SE) > qnorm(1-thresh/2)) 
saveRDS(res, file="results.oneplus.RDS")

res$name <- factor(res$name)
res$sig <- factor(res$sig, levels=c("FALSE", "TRUE"))
plt <- ggplot(res) + geom_point(aes(x=name, y=Est, color=sig)) + 
            geom_errorbar(aes(x=name, ymin=CI_Lower, ymax=CI_Upper, color=sig)) + 
            geom_hline(yintercept=1, linetype=2) + 
            ylab("Log Fold Enrichment") + 
            #geom_vline(xintercept=ix, color="#F8766D") + 
            theme_bw() + theme(axis.text.x = element_text(angle=90), 
                               axis.title.x=element_blank(), 
                               legend.position="none") 

ggsave(plt, file="ga.results.oneplus.png", height=4, width=6.5, units="in", dpi=300)

##### Full model enrichment
rst_full <- read_lines("ga.full.rst")
res_full <- apply(name_dict, 1, function(x){
                  line <- grep(x[1], rst_full)
                  rst <- rst_full[line]
                  rst <- unlist(str_split(rst, " "))
                  rst <- rst[str_length(rst) > 0]
                  est <- as.numeric(rst[2])
                  ci_lower <- as.numeric(rst[3])
                  ci_upper <- as.numeric(rst[4])
                  se <- se_from_ci(ci_upper-ci_lower)
                  return(c( est, ci_lower, ci_upper, se))
               })
res_full <- data.frame(t(res_full), stringsAsFactors=FALSE)
names(res_full) <- c("Est", "CI_Lower", "CI_Upper", "SE")

thresh <- 0.05/nrow(res)
res_full <- res_full %>% mutate(CI_Lower = Est -qnorm(1-thresh/2, sd=SE), 
              CI_Upper = Est +qnorm(1-thresh/2, sd=SE), 
              sig = abs(Est/SE) > qnorm(1-thresh/2)) 
res_full$name <- factor(name_dict$name)
res_full$sig <- factor(res_full$sig, levels=c("FALSE", "TRUE"))
saveRDS(res_full, file="results.full.RDS")


plt <- ggplot(res_full) + geom_point(aes(x=name, y=Est, color=sig)) + 
            geom_errorbar(aes(x=name, ymin=CI_Lower, ymax=CI_Upper, color=sig)) + 
            geom_hline(yintercept=1, linetype=2) + 
            ylab("Log Fold Enrichment") + 
            #geom_vline(xintercept=ix, color="#F8766D") + 
            theme_bw() + theme(axis.text.x = element_text(angle=90), 
                               axis.title.x=element_blank(), 
                               legend.position="none") 

ggsave(plt, file="ga.results.full_dsc.png", height=4, width=6.5, units="in", dpi=300)

res_full_temp <- map(rst_full[2:61], function(x){
                  x <- unlist(str_split(x, " "))
                  x <- x[str_length(x) > 0]
                  se <- se_from_ci(as.numeric(x[4])-as.numeric(x[3]))
                  x <- c(x, se)
                  #names(x) <- c("Name", "Est", "CI_Lower", "CI_Upper", "SE")
                  return(x)
               })
res_full <- data.frame(Name = map(res_full_temp, 1) %>% unlist, 
                       Est = map(res_full_temp, 2) %>% unlist, 
                       CI_Lower = map(res_full_temp, 3) %>% unlist, 
                       CI_Upper = map(res_full_temp, 4) %>% unlist, 
                       CI_SE = map(res_full_temp, 5) %>% unlist, stringsAsFactors=FALSE) 
res_full$Name <- str_replace(res_full$Name, fixed(".1"), "")
saveRDS(res_full, file="results.full_all.RDS")

#####FDR Results
fdr_results <- readRDS("results.RDS")
fdr_results <- fdr_results %>% mutate(color = case_when(FDR_full < 0.1 & FDR_base > 0.1 ~ 1,
                                     FDR_full < 0.1 & FDR_base < 0.1 ~ 4,
                                     FDR_full > 0.1 & FDR_base < 0.1 ~ 4,
                                     FDR_full > 0.1 & FDR_base > 0.1 ~ 4), 
                   color = factor(color))

plt <- ggplot(fdr_results) + geom_point(aes(x=FDR_base, y=FDR_full, col=color)) + 
        scale_color_manual(values=c("#58C14A", "black")) + 
        xlab("FDR with Baseline Annotations") + ylab("FDR with Baseline and Endometrial Annotations") + 
        geom_hline(yintercept=0.1, linetype=3) + geom_vline(xintercept=0.1, linetype=3) + 
        theme_bw() + theme(legend.position="none")


ggsave(plt, file="fdr.png", height=5, width=7, units="in", dpi=300)

#Region 614 is EBF1 on Chromosome 5
#Region 15 is WNT4 locus on Chrmosome 1
#Region 356 is EEFSEC on Chromosome 3
fdr_results <- fdr_results %>% mutate(color = case_when(
                                     region %in% c(614, 15, 356) ~ 5, 
                                     FDR_full < 0.1 & FDR_base > 0.1 ~ 1,
                                     FDR_full < 0.1 & FDR_base < 0.1 ~ 2,
                                     FDR_full > 0.1 & FDR_base < 0.1 ~ 3,
                                     FDR_full > 0.1 & FDR_base > 0.1 ~ 4),
                   color = factor(color))

plt <- ggplot(fdr_results) + geom_point(aes(x=FDR_base, y=FDR_full, col=color)) + 
        #scale_color_manual(values=c("#58C14A", "black")) + 
        xlab("FDR with Baseline Annotations") + ylab("FDR with Baseline and Endometrial Annotations") + 
        geom_hline(yintercept=0.1, linetype=3) + geom_vline(xintercept=0.1, linetype=3) + 
        theme_bw() + theme(legend.position="none")

ggsave(plt, file="fdr.png", height=5, width=7, units="in", dpi=300)

plts <- list()
plts[[1]] <- ggplot(fdr_results) + geom_point(aes(x=FDR_base, y=FDR_full ), alpha=0.5) + 
        #scale_color_manual(values=c("#E41A1C", "black")) + 
        xlab("FDR with Baseline(Gazal et al 2017) Annotations") + ylab("FDR with Baseline + Endometrial Annotations") + 
        geom_hline(yintercept=0.1, linetype=3) + geom_vline(xintercept=0.1, linetype=3) + 
        theme_bw() + theme(legend.position="none", axis.title=element_text(size=18))

plts[[2]] <- ggplot(fdr_results) + geom_point(aes(x=FDR_base, y=FDR_full, color=region %in% c(614, 15, 356) ), alpha=0.5) + 
        scale_color_manual(values=c("black", "#E41A1C")) + 
        xlab("FDR with Baseline(Gazal et al 2017) Annotations") + ylab("FDR with Baseline + Endometrial Annotations") + 
        geom_hline(yintercept=0.1, linetype=3) + geom_vline(xintercept=0.1, linetype=3) + 
        theme_bw() + theme(legend.position="none", axis.title=element_text(size=18))

plts[[3]] <- ggplot(fdr_results) + geom_point(aes(x=FDR_base, y=FDR_full, color=FDR_base < 0.1 ), alpha=0.5) + 
        scale_color_manual(values=c("black", "#E41A1C")) + 
        xlab("FDR with Baseline(Gazal et al 2017) Annotations") + ylab("FDR with Baseline + Endometrial Annotations") + 
        geom_hline(yintercept=0.1, linetype=3) + geom_vline(xintercept=0.1, linetype=3) + 
        theme_bw() + theme(legend.position="none", axis.title=element_text(size=18))


plts[[4]] <- ggplot(fdr_results) + geom_point(aes(x=FDR_base, y=FDR_full, color=FDR_full < 0.1 ), alpha=0.5) + 
        scale_color_manual(values=c("black", "#E41A1C")) + 
        xlab("FDR with Baseline(Gazal et al 2017) Annotations") + ylab("FDR with Baseline + Endometrial Annotations") + 
        geom_hline(yintercept=0.1, linetype=3) + geom_vline(xintercept=0.1, linetype=3) + 
        theme_bw() + theme(legend.position="none", axis.title=element_text(size=18))

for(i in 1:4) ggsave(plts[[i]], file=paste0("fdr", i, ".png"), height=6.2, width=8.7, units="in", dpi=300)


