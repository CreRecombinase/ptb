library(readr)
library(dplyr)

annot_for_caviarbf <- function(zscore_file_list, annot_rds, output_file_list, header_file){
    stopifnot(length(zscore_file_list) == length(output_file_list))
    annot <- readRDS(annot_rds) 
    for(i in seq_along(zscore_file_list)){
        zscores <- read_tsv(zscore_file_list[i], col_names=c("SNP", "zscore"))
        my_annot <- left_join(zscores, annot, by="SNP") %>% select(-zscore)
        stopifnot(all(my_annot$SNP == zscores$SNP))
        my_annot <- select(my_annot, -SNP)
        write_tsv(my_annot, output_file_list[i], col_names=FALSE)
    }
    write_lines(names(my_annot), header_file)

}

library(stringr)
zscore_file_list <- list.files("zscores2/", "zscores.filtered", full.names=TRUE)
output_file_list <- str_replace(zscore_file_list, fixed("zscores."), "fullannot.")
header_file <- "zscores2/fullannot.names"
annot_rds <- "../one_plus_annot_torus/full_annotations.torus.RDS"
annot_for_caviarbf(zscore_file_list, annot_rds, output_file_list, header_file)

if(FALSE){
#mkdir caviarbf_files
#cp zscores2/zscores.filtered.* caviarbf_files/
#cp plink_output/*filtered*ld caviarbf_files/
#cp zscores2/fullannot.* caviarbf_files/
fn1 <- list.files(".", "zscores")
fn <- str_replace(fn1, fixed("zscores.filtered."), "")
fn <- str_replace(fn, fixed(".tsv"), ".filtered")
cmds <- paste0("mv ", fn1, " ", fn)
sapply(cmds, function(x){system(x)})

fn1 <- list.files(".", "ld")
fn <- str_replace(fn1, "ld", "LD")
cmds <- paste0("mv ", fn1, " ", fn)
sapply(cmds, function(x){system(x)})

fn1 <- list.files(".", "fullannot")
fn <- str_replace(fn1, fixed("fullannot.filtered."), "")
fn <- str_replace(fn, fixed(".tsv"), ".filtered.annotations")
cmds <- paste0("mv ", fn1, " ", fn)
sapply(cmds, function(x){system(x)})
}

