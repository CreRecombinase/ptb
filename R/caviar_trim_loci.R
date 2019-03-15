library(readr)
caviar_trim_loci <- function(prefix, out_prefix, trim_pvalue = 0.001, thresh_pvalue){
    zscores <- read_tsv(prefix, col_names=c("SNP", "zscore"))
    trim <- abs(qnorm(trim_pvalue/2))
    if(all(abs(zscores$zscore) < trim)){
        cat("All zscores below threshold\n")
        return(0)
    }
    ix1 <- match(TRUE, abs(zscores$zscore) > trim)
    ix2 <- nrow(zscores) - match(TRUE, rev(abs(zscores$zscore) > trim)) + 1
    thresh <- abs(qnorm(thresh_pvalue/2))
    ix <- which(abs(zscores$zscore) > thresh)
    ix_keep <- ix[ix >= ix1 & ix <= ix2]
    cat("Retaining ", length(ix_keep), " SNPs of ", nrow(zscores), "\n")
    zscores <- zscores[ix_keep,]
    write_tsv(zscores, path=out_prefix, col_names=FALSE)

    ld <- read_tsv(paste0(prefix, ".LD"), col_names=FALSE)
    ld <- as.matrix(ld)
    ld <- ld[ix_keep, ix_keep]
    ld <- data.frame(ld)
    write_tsv(ld, path=paste0(out_prefix, ".LD"), col_names=FALSE)

    annot <- read_tsv(paste0(prefix, ".annotations"), col_names=FALSE)
    annot <- annot[ix_keep,]
    write_tsv(annot, path=paste0(out_prefix, ".annotations"), col_names=FALSE)
}



if(FALSE){
    library(stringr)
fl <- list.files(".", fixed("filtered.LD"))
fl <- str_replace(fl, fixed(".LD"), "")
fl <- str_replace(fl, fixed(".trim0.0001"), "")
fl <- unique(fl)
trim_pval <- 0.001
thresh_pval <- 0.05
dir <- paste0("trim", trim_pval, "_thresh", thresh_pval)
system(paste0("mkdir ", dir))
new_prefix <- paste0(dir, "/", fl, ".", dir)
for(i in seq_along(fl)){
    cat(fl[i], "\n")
    caviar_trim_loci(fl[i], new_prefix[i], trim_pvalue=trim_pval, thresh_pvalue=thresh_pval)
}
x <- paste0(fl, ".", dir)
write_lines(x, path=paste0(dir, "/ga_loci.txt"))
for(i in seq_along(x)){
    write_lines(x[i], path=paste0(dir, "/ga_locus.", i, ".txt"))
}
#test caviarbf
library(caviarbf)
times <- list()
for(i in 4:18){
    times[[i]] <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_locus.", i, ".txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "lasso",  hyperParamSelection = "topK", K = 10, 
                                useParallel = F,
                                outputPrefix =paste0(dir, "/test_", i), 
                                eps = 0.01, 
                                annotationIndices = -1)
    )
}
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  hyperParamSelection = "topK", K = 50, 
                                useParallel = F,
                                outputPrefix =paste0(dir, "/onecausal_topK"), 
                                eps = 0.01, 
                                annotationIndices = (1:97)[-70])
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  hyperParamSelection = "cv",
                                useParallel = F,
                                outputPrefix =paste0(dir, "/onecausal_lasso_cv"), 
                                eps = 0.01, 
                                annotationIndices = (1:97)[-70])
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 3, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  hyperParamSelection = "cv",
                                useParallel = F,
                                outputPrefix =paste0(dir, "/threecausal_lasso_cv"), 
                                eps = 0.01, 
                                annotationIndices = (1:97)[-70])
    )

dir <- "trim0.001"
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  hyperParamSelection = "cv",
                                useParallel = F,
                                outputPrefix =paste0(dir, "/onecausal_lasso_cv"), 
                                eps = 0.01, 
                                annotationIndices = (1:97)[-70])
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  hyperParamSelection = "topK", K = 50,
                                useParallel = F,
                                outputPrefix =paste0(dir, "/onecausal_topK"), 
                                eps = 0.01, 
                                annotationIndices = (1:97)[-70])
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetENETMin",  hyperParamSelection = "cv",
                                useParallel = F,
                                outputPrefix =paste0(dir, "/onecausal_enet_cv"), 
                                eps = 0.01, 
                                annotationIndices = (1:97)[-70])
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetENETMin",  hyperParamSelection = "cv",
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/onecausal_base1only_enet_cv"), 
                                eps = 0.01, 
                                annotationIndices = base1_only)
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  hyperParamSelection = "topK", K = 100,
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/onecausal_lasso_topK_corr0.5"), 
                                pvalueThreshold=1, rThreshold = 0.5,
                                eps = 0.01, 
                                annotationIndices = (1:97)[-70])
    )

base1_ix <- c(1,2,4,6,8,10,11,13,15,17,19,21,23,24,26,27,29,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,59,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97)
saveRDS(base1_ix, file="caviarbf_files/base1_ix.RDS")
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/onecausal_base1_lasso_lam0"), 
                                lambda = 0,
                                eps = 0.01, 
                                annotationIndices = base1_ix)
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/onecausal_base1m53_lasso_lam0"), 
                                lambda = 0,
                                eps = 0.01, 
                                annotationIndices = base1_ix[-53])
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/onecausal_base1_lasso_lam1e-303"), 
                                lambda = 1e-303,
                                eps = 0.01, 
                                annotationIndices = base1_ix)
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                hyperParamSelection = "topK", K = 100,
                                outputPrefix =paste0(dir, "/onecausal_base1_lasso_topKcorr0.5"), 
                                pvalueThreshold=1, rThreshold = 0.5,
                                eps = 0.01, 
                                annotationIndices = base1_ix)
    )

t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetENETMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/onecausal_base1m53_enet_cv"), 
                                eps = 0.01, 
                                annotationIndices = base1_ix[-53])
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/onecausal_base1only_lasso_cv"), 
                                eps = 0.01, 
                                annotationIndices = base1_ix[1:37])
    )

t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/onecausal_noannot"), 
                                lambda = 1,
                                eps = 0.01, 
                                annotationIndices = c(1))
    )


base1_only <- base1_ix[1:37]
saveRDS(base1_only, file="base1_only.RDS")
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                hyperParamSelection = "topK", K = 100,
                                outputPrefix =paste0(dir, "/onecausal_base1_lasso_topKcorr0.5pval0.5"), 
                                pvalueThreshold=0.5, rThreshold = 0.5,
                                eps = 0.01, 
                                annotationIndices = base1_ix)
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/onecausal_base1only_lasso_lam0"), 
                                lambda = 0,
                                eps = 0.01, 
                                annotationIndices = base1_only)
    )
t <- system.time(
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/test_lam0"), 
                                lambda = 0,
                                eps = 0.01, 
                                annotationIndices = c(1, 2, 4, 97))
    )

ix <- read_lines("caviarbf_files/trim0.001_thresh0.05/onecausal_base1only_lasso_topKcorr0.5.topk")
ix <- as.numeric(ix)
base1_only_keep <- base1_only[ix]
base1_keep <- c(base1_only[ix], base1_ix[38:59])
saveRDS(base1_keep, file="base1_keep.RDS")
for(i in 38:59){
    cat(i, "\n")
    my_ix <- c(base1_only_keep, base1_ix[i])
    t <- system.time(
        res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), 
                                maxCausal = 1, nSample = 43568,
                                fittingMethod = "glmnetLASSOMin",  
                                useParallel = T, ncores=8,
                                outputPrefix =paste0(dir, "/adding_dsc_", i), 
                                lambda=0,
                                eps = 0.01, 
                                annotationIndices = my_ix)
    )
}

screen sinteractive --partition=broadwl --account=pi-xinhe --cpus-per-task=8 --ntasks=1 --mem-per-cpu=7G --time=10:00:00
}
