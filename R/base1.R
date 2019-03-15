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
#write_tsv(full_annot, "full_annotation.torus.tsv.gz")

zscores <- filter(zscores, variant %in% full_annot$SNP)
o <- match(full_annot$SNP, zscores$variant)
zscores <- zscores[o,]

write_tsv(zscores, "../torus/gestational_age/ga.zscores.tsv.gz")


#We need to reduce the background model so that torus can run in a reasonable amount of time
base_annotations <- base_annotations[grep("extend", base_annotations, invert=TRUE)]
base_annotations <- base_annotations[grep("MAF", base_annotations, invert=TRUE)]

#dsc_annot_ix <- match(dsc_annotations, names(full_annot))
#base_ix <- match(base_annotations, names(full_annot))
#for(i in seq_along(dsc_annotations)){
#    an <- full_annot[c(1, base_ix, dsc_annot_ix[i])]
#    nm <- dsc_annotations[i]
#    cat(nm, " ")
#    write_tsv(an, path=paste0(nm, ".oneplus.annot.tsv.gz"))
#}

full_annot_keep <- full_annot[c(1, base_ix, dsc_annot_ix)]
write_tsv(full_annot_keep, path="full.annot.tsv.gz")


full_annot <- read_tsv("full.annot.tsv.gz")
names(full_annot) <- paste0(names(full_annot), 
                            rep(c("", "_c", "_d", "_c", "_d"), c(1, 1, 32, 4, 22)))
write_tsv(full_annot, path="full.annot.tsv.gz")

base <- full_annot[, 1:38]
write_tsv(base, path="base.annot.tsv.gz")


##########
full_annot <- read_tsv("full.annot.tsv.gz")
keep_col <- c(1:38, 39:40, 43:44, 53:60)
partial <- full_annot[, keep_col]
write_tsv(partial, path="base_and_dsc_noTF.annot.tsv.gz")
