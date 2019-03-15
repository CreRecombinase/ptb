library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

write_caviar_wig <- function(marginalz_file, data_rds_list, 
                             locus_name_list, track_name){
    out_prefix <-  str_replace(marginalz_file, fixed(".marginalz"), "")


    pips <- read_delim(marginalz_file, 
                       col_names=c("ix", "pip", "locus", "snp", "zscore"), 
                       delim=" ")
    stopifnot(length(data_rds_list) == length(locus_name_list))
    for( i in seq_along(data_rds_list)){
        fn <- paste0(out_prefix, ".", locus_name_list[i],  ".wig")
        dat <- readRDS(data_rds_list[i]) %>% rename(snp = variant)
        header <- c(paste0("browser position chr", dat$CHR[1], ":", 
                           min(dat$position), "-", max(dat$position)), 
                   paste0("track type=wiggle_0 name=", track_name, " description=\"caviarbf PIP\" visibility=full autoScale=off viewLimits=0:1 color=23,128,48 priority=10"), 
                   paste0("variableStep chrom=chr", dat$CHR[1], " span=1"))
        
        write_lines(header, path=fn)
        wig_dat <-  inner_join(pips, dat) %>% 
               select(position, pip)  %>% 
               mutate(position =  as.integer(position)) %>%
               arrange(position)
        write.table(wig_dat, file=fn,
                    append=TRUE, col.names=FALSE, row.names=FALSE, sep="\t")
    }
}

data_rds_list <- list.files("zscores2", "data.filtered", full.names=TRUE)
locus_name_list <- str_replace(data_rds_list, "zscores2/data.filtered.", "")
locus_name_list <- str_replace(locus_name_list, ".RDS", "")
