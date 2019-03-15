library(tidyverse)
library(spaa)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
blocks <- as.numeric(args)
#b <- 512

zscores <- read_tsv("ga.zscores.tsv.gz")

snp_info <- readRDS("/project/compbio/LD/CEU/SNP_chunks.RDS")

block2chunk <- readRDS("/project/compbio/LD/CEU/block2chunk.RDS")

#annot <- read_tsv("../../one_plus_annot_torus/base1/full.annot.tsv.gz")

for(b in blocks){
    myz <- filter(zscores, ldchunk==b) %>%
       rename(rs = variant) %>%
       inner_join(., snp_info)
    chunks <- unique(myz$chunk)
    blocks <- filter(block2chunk,  row_chunk %in% chunks & col_chunk %in% chunks)

    chr <- myz$chr[1]
    files <- paste0("/project/compbio/LD/CEU/CEU/LD_DF/chr", chr, "/AF0.01chr", chr, "_CEU_F_omni_T_0.01_", blocks$block_ind, "_190.RDS")
    ld <- map_df(files, function(x){
             readRDS(x)})
    ld <- filter(ld, rowsnp %in% myz$rs & colsnp %in% myz$rs)

#myz.save <- myz
#ld.save <- ld
#myz <- myz[1:10,]
#ld <- filter(ld, rowsnp %in% myz$rs & colsnp %in% myz$rs)

    ld$ix1 <- match(ld$rowsnp, myz$rs)
    ld$ix2 <- match(ld$colsnp, myz$rs)

    cor_mat <- sparseMatrix(i=ld$ix1,j=ld$ix2,x=ld$r,symmetric=T)
    cor_mat <- data.frame(as.matrix(cor_mat))
    write_tsv(cor_mat, path=paste0("ld.", b, ".txt"), col_names=FALSE)

    myz <- myz %>% rename(variant=rs) %>% select(variant, ldchunk, tstat)
    write_tsv(myz, path=paste0("ga.zscores.", b, ".tsv"))
}

