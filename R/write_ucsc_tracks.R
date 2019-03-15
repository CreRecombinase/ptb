library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
blocks <- as.numeric(args)

#block <- 512

zscores <- read_tsv("ga.zscores.tsv.gz")


ld_files <- paste0("/project2/xinhe/LD/CEU/LD_DF/chr", 1:22, "AF0.05SNP0N0_CEU_CEU_F_CEU_T_0.01.RDS")
snp_info_files <- paste0("/project2/xinhe/LD/CEU/LD_DF_SNPList/chr", 1:22, "AF0.05SNP0N0_CEU_CEU_F_CEU_T_0.01.RDS")


snp_info <- map_df(snp_info_files, function(x){readRDS(x)})


#zz <- filter(zscores, ldchunk==block) %>% 
#      rename(SNP = variant) %>%
#      inner_join(., snp_info, by="SNP")
#

annot <- read_tsv("../../one_plus_annot_torus/base1/full.annot.tsv.gz")
#annot_dsc <- annot[, c(1, 38:60)]
#zz <- left_join(zz, annot_dsc, by="SNP")



#filter(zz, abs(tstat) > 2.5) %>% 
#       select(-SNP, -ldchunk, -tstat, -AF, -allele, -chr, -pos, -snp_id, 
#              -region_id, -map, -ld_snp_id) %>%
#      colSums(.)


#ggplot(zz) + geom_point(aes(x=pos, y=tstat, col=eQTL_0.05_FDR_d))


#annot1 <- read_tsv("../../one_plus_annot_torus/base1/eQTL_0.05_FDR.oneplus.annot.tsv.gz")
#yy <- filter(zscores, ldchunk==block) %>% 
#      rename(SNP = variant) %>%
#      inner_join(., snp_info, by="SNP") %>%
#      left_join(., annot1, by="SNP")
#ggplot(yy) + geom_point(aes(x=pos, y=tstat, col=eQTL_0.05_FDR_d))

#Write bed files for UCSC
assoc <- read_tsv("/project/mstephens/data/external_public_supp/gwas_summary_statistics/raw_summary_statistics/23_and_me_ptb/23andMe_ga.tsv.gz")
res <- readRDS("torus_base1/results.RDS")

for(block in blocks){
cat(block, " ")
dat <- filter(zscores, ldchunk==block) %>% 
       rename(snp = variant) %>%
       left_join(., assoc) %>%
       rename(SNP = snp) %>%
       left_join(., annot)

#ggplot(dat) + geom_point(aes(x=position, y=tstat, col=`h3k27ac-final-D-48h_d`))
fdr1 <- res$FDR_base[res$region==block]
fdr2 <- res$FDR_full[res$region==block]


zscore_header <- c(paste0("browser position ", dat$chr[1], ":", min(dat$position), "-", max(dat$position)), 
                   paste0("track type=wiggle_0 name=gestational_age description=\"Negative log10 p-value, FDR_base=", fdr1, " FDR_full=", fdr2, "\" visibility=full autoScale=off viewLimits=0:7 color=50,150,255 priority=10"), 
                   paste0("variableStep chrom=", dat$chr[1], " span=1"))
write_lines(zscore_header, path=paste0("zscores.", block, ".wig"))
zscore_dat <-  dat %>% mutate(neglogp = -log10(2*pnorm(-abs(tstat)))) %>% 
               select(position, neglogp)  %>% 
               mutate(position =  as.integer(position))
write_tsv(zscore_dat, path=paste0("zscores.", block, ".wig"), append=TRUE, col_names=FALSE)


dat <- mutate(dat, p_value = 2*pnorm(-abs(tstat))) 

dsc_annot <- names(dat[49:70])
dsc_annot <- str_replace(dsc_annot, "_d$", "")
dsc_annot <- str_replace(dsc_annot, "_c$", "")
bed_files <- paste0("/project2/mstephens/ptb/dsc_annotation/", dsc_annot, ".bed")

for(i in seq_along(dsc_annot)){
    header = c(paste0("browser position ", dat$chr[1], ":", min(dat$position), "-", max(dat$position)), 
               paste0("track name=", dsc_annot[i], " description=\"", dsc_annot[i], "\" visibility=1"))
    write_lines(header, path=paste0("dsc_annot.", block, ".bed"), append=i>1)
    bed <- read_tsv(bed_files[i], col_names=FALSE)
    bed <- filter(bed, X1==dat$chr[1] & X2 >= min(dat$position) & X3 <= max(dat$position)) %>%
            select(X1, X2, X3)
    write_tsv(bed, path=paste0("dsc_annot.", block, ".bed"), append=TRUE, col_names=FALSE)

}


}

