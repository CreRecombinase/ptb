library(readr)
library(caviarbf)

args <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])


dir <- "trim0.001_thresh0.05"

base1_ix <- c(1,2,4,6,8,10,11,13,15,17,19,21,23,24,26,27,29,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,59,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97)
base1_only <- base1_ix[1:37]

ix <- read_lines("trim0.001_thresh0.05/onecausal_base1only_lasso_topKcorr0.5.topk")
ix <- as.numeric(ix)
base1_only_keep <- base1_only[ix]


cat(i, "\n")
if(i ==0){ my_ix <- base1_only_keep
} else{my_ix <- c(base1_only_keep, base1_ix[i])}

res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/adding_dsc_", i), 
                                lambda=0,
                                eps = 0.01, 
                                annotationIndices = my_ix)

