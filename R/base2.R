library(tidyverse)


base_files <- paste0("../1000G_ldscores/1000G_Phase3_baselineLD_ldscores/baselineLD.", 1:22, ".annot.tsv.gz")
ldsc_files <- paste0("../1000G_ldscores/1000G_Phase3_baselineLD_ldscores_full/baselineLD.1000G.EUR.", 1:22, ".l2.ldscore.gz")
dsc_files <- paste0("../dsc_annotation/dsc.", 1:22, ".annot.gz")


base_annot <- map_df(base_files, function(x){read_tsv(x)})
dsc_annot <- map_df(dsc_files, function(x){read_tsv(x)})
ldsc <- map_df(ldsc_files, function(x){read_tsv(x)})
ldsc <- select(ldsc,  SNP, baseL2)


#data <- read_tsv("/project/mstephens/data/external_public_supp/gwas_summary_statistics/raw_summary_statistics/23_and_me_ptb/23andMe_ga.tsv.gz")
zscores <- read_tsv("/project2/xinhe/jean/ptb/torus/23andMe_ga/23andMe_ga.zscores.tsv.gz")


dsc_annotations <- names(dsc_annot)[-c(1:3)]
base_annotations <- c("baseL2", names(base_annot)[-c(1:5)])


full_annot <- inner_join(base_annot, dsc_annot)
full_annot <- full_annot %>% select(-CHR, -BP, -CM, -base) %>% 
                filter(SNP %in% zscores$variant)
full_annot <- inner_join(ldsc, full_annot)

dim(full_annot)
#[1] 8862220      98

cont_ix <- apply(full_annot[1:100,], 2, function(x){any(!x %in% c(0, 1))})
cont_ix <- as.numeric(which(cont_ix))[-1]
disc_ix <- apply(full_annot[1:100,], 2, function(x){all(x %in% c(0, 1))})
disc_ix <- as.numeric(which(disc_ix))

nms <- names(full_annot)
nms[cont_ix] <- paste0(nms[cont_ix], "_c")
nms[disc_ix] <- paste0(nms[disc_ix], "_d")
names(full_annot) <- nms


saveRDS(full_annot, file="full_annotations.torus.RDS")

#write_tsv(full_annot, "full_annotation.torus.tsv.gz")

zscores <- filter(zscores, variant %in% full_annot$SNP)
o <- match(full_annot$SNP, zscores$variant)
zscores <- zscores[o,]

write_tsv(zscores, "../torus/gestational_age/ga.zscores.tsv.gz")

base_annotations <- names(full_annot)[2:76]
dsc_annotations <- names(full_annot)[77:98]
#We need to reduce the background model so that torus can run in a reasonable amount of time
#Base1
#base_annotations <- base_annotations[grep("extend", base_annotations, invert=TRUE)]
#base_annotations <- base_annotations[grep("MAF", base_annotations, invert=TRUE)]

#Base2
extend <- base_annotations[grep("extend", base_annotations)]
ext_name <- stringr::str_replace(extend, ".extend.500", "")
base_annotations <- base_annotations[!base_annotations %in% ext_name]

base_ix <- match(base_annotations, names(full_annot))
dsc_ix <- match(dsc_annotations, names(full_annot))


full_annot_keep <- full_annot[c(1, base_ix, dsc_ix)]
write_tsv(full_annot_keep, path="base2/full.annot.tsv.gz")

base <- full_annot[, c(1, base_ix)]
write_tsv(base, path="base2/base.annot.tsv.gz")

#######End#######
