library(tidyverse)

cbf_full <- read_delim("caviarbf_files/trim0.001_thresh0.05/onecausal_base1_lasso_lam0.gamma", 
                       delim=" ", col_names=FALSE)
nms <- read_lines("caviarbf_files/fullannot.names")
base1_ix <- c(1,2,4,6,8,10,11,13,15,17,19,21,23,24,26,27,29,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,59,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97)
nms <- nms[base1_ix]
nms <- str_replace(nms, "_c$", "")
nms <- str_replace(nms, "_d$", "")
nms <- c("Intercept", nms)
cbf_full$Name <- nms


torus_full <- readRDS("torus_base1/results.full_all.RDS")

torus_full$Est <- as.numeric(torus_full$Est)
res <- full_join(torus_full, cbf_full, by="Name") %>%
        select(-X2) %>% rename(cbf_full = X1)

cbf_base1m53 <- read_delim("caviarbf_files/trim0.001_thresh0.05/onecausal_base1m53_lasso_lam0.gamma", 
                       delim=" ", col_names=FALSE)
cbf_base1m53$Name <- nms[-54]

res <- full_join(res, cbf_base1m53) %>% select(-X2) %>% rename(cbf_m53=X1)

base1_only <- base1_ix[1:37]
saveRDS(base1_only, file="caviarbf_files/base1_only.RDS")
ix <- read_lines("caviarbf_files/trim0.001_thresh0.05/onecausal_base1only_lasso_topKcorr0.5.topk")
ix <- as.numeric(ix)
base1_only_keep <- base1_only[ix]
base1_keep <- c(base1_only[ix], base1_ix[38:59])
saveRDS(base1_keep, file="caviarbf_files/base1_keep.RDS")
z <- 38:59
z <- z[z!=53]
my_ix <- c(base1_only_keep, base1_ix[z])
cbf_base1keepm53 <- read_delim("caviarbf_files/trim0.001_thresh0.05/missing_dsc_53.gamma", 
                       delim=" ", col_names=FALSE)
nms2 <- read_lines("caviarbf_files/fullannot.names")
nms2 <- str_replace(nms2, "_c$", "")
nms2 <- str_replace(nms2, "_d$", "")
cbf_base1keepm53$Name <-c("Intercept", nms2[my_ix])
res <- full_join(res, cbf_base1keepm53) %>% select(-X2) %>% rename(cbf_base1keepm53=X1)


##Loooking up variants and comparing pips  
annots <- readRDS("../one_plus_annot_torus/full_annotations.torus.RDS")


cbf_null <- read_delim("caviarbf_files/trim0.001_thresh0.05/onecausal_noannot.marginalz",
                       col_names=c("ix", "pip", "locus", "snp", "zscore"), 
                       delim=" ") %>% 
            rename(cbf_null_pip = pip)
cbf_base <- read_delim("caviarbf_files/trim0.001_thresh0.05/onecausal_base1only_lasso_topKcorr0.5.marginalz",
                       col_names=c("ix", "pip", "locus", "snp", "zscore"), 
                       delim=" ") %>%
            rename(cbf_base_pip = pip)
cbf_full <- read_delim("caviarbf_files/trim0.001_thresh0.05/missing_dsc_53.marginalz",
                       col_names=c("ix", "pip", "locus", "snp", "zscore"), 
                       delim=" ") %>%
            rename(cbf_full_pip = pip)

cbf <- inner_join(cbf_null, cbf_base)
cbf <- inner_join(cbf, cbf_full) %>% rename(SNP = snp)

r <- 614
r <-15 
r <- 512
dat <- readRDS(paste0("zscores2/data.filtered.", r, ".RDS")) %>% 
            select(CHR, variant, position, tstat, effect, stderr, pvalue) %>%
            rename(SNP = variant)
dat <- left_join(dat, annots, by=c("SNP"))
susie <- readRDS(paste0("susie_finemap/susie.summary.", r, ".RDS")) %>% rename(SNP = variant)
dat <- left_join(dat, susie, by=c("CHR", "SNP", "position"))
dat <- left_join(dat, cbf, by=c("SNP"))

hic <- read_tsv("../Hi-C/DT1_dTL4_D_48h.ibed")

p <- 22470407
filter(hic, otherEnd_chr=="chr1" & otherEnd_start < p & otherEnd_end > p)


