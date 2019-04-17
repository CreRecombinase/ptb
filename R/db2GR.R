library(dbplyr)
library(tidyverse)
library(ldshrink)

input_db <- snakemake@input[["input_f"]]
input_f <- snakemake@input[["ldf"]]
min_size <- as.integer(snakemake@params[["ld_reg"]])
outf <- snakemake@output[["outf"]]

beta_v <- snakemake@params[["beta_v"]]
se_v <- snakemake@params[["se_v"]]
stopifnot(!is.null(beta_v),!is.null(se_v))
save.image(paste0("ts.",beta_v,"RData"))

ld_df <- read_tsv(input_f, col_types = cols(
                              chr = col_character(),
                              start = col_integer(),
                              stop = col_integer()
                          )) %>%
    rename(chrom = chr) %>%
    group_by(chrom) %>%
    mutate(start = if_else(start == min(start),
                           0L, start)) %>%
    mutate(stop = if_else(stop == max(stop),
                          .Machine$integer.max, stop)) %>%
    ungroup() %>%
    mutate(region_id = 1:n())


stopifnot(file.exists(input_db))
output_txtf <- snakemake@output[["txtf"]]
output_rds <- snakemake@output[["rds"]]
gwas_df <- dplyr::tbl(dplyr::src_sqlite(path = input_db, create = F),
                      "gwas") %>%
    collect()
gwas_z <- dplyr::select(gwas_df,
                        SNP = id,
                        chrom, pos,
                        beta = !! beta_v,
                        se = !! se_v) %>%
    dplyr::mutate(`z-val` = beta / se)


m_gwas_df <- select(gwas_df,
                    chrom,
                    pos,
                    MAJOR = A1,
                    MINOR = A2) %>%
    mutate(chrom = as.integer(chrom),
           pos = as.integer(pos),
           SNP = paste0(chrom, ":", pos))

snp_df <- gwas_z %>%
    select(SNP, chr = chrom, pos, `z-val`) %>%
    mutate(chr = as.integer(gsub("chr", "", chr)), pos = as.integer(pos))  %>%
    distinct(chr, pos, .keep_all = T) %>%
    arrange(chr, pos)
snp_df <- mutate(snp_df,
                 region_id = assign_snp_block(break_chr = ld_df$ch,
                                              break_start = ld_df$start,
                                              break_stop = ld_df$stop,
                                              break_id = ld_df$region_id,
                                              snp_chr = chr,
                                              pos = pos,
                                              assign_all = T))

snp_df <- dplyr::group_by(snp_df, chr, region_id) %>%
    dplyr::summarise(ct = n()) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(snp_df) %>%
    dplyr::mutate(region_id = ifelse(ct < min_size,
                                     region_id - 1,
                                     region_id)) %>%
    dplyr::select(-ct) %>%
    dplyr::ungroup()

snp_df %>%
    select(SNP, locus = region_id, `z-val`) %>%
    write_tsv(output_txtf)
gr_df <- select(gwas_z, chrom, start = pos) %>%
    mutate(chrom = paste0("chr", chrom),
           start = as.integer(start),
           end = start + 1L) %>%
    GenomicRanges::makeGRangesFromDataFrame()

saveRDS(gr_df, output_rds)
