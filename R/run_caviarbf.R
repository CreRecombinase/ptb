library(caviarbf)

args <- commandArgs(trailingOnly=TRUE)
fit <- args[1]
cores <- as.numeric(args[2])
eps <- as.numeric(args[3])
selection <- args[4]
dir <- args[5]
lambda <- as.numeric(args[6])

if(fit == "lasso"){
    fit_method <- "glmnetLASSOMin"
}else if(fit == "enet"){
    fit_method <- "glmnetENETMin"
}else{
    stop("fit method: ", fit, " not valid. Use lasso or enet\n")         
}
out_prefix <- paste0(dir, "/", fit, "_eps", eps, "_", selection)
cat(out_prefix, "\n")

#res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), maxCausal = 5, nSample = 43568,
#                            fittingMethod = fit_method,  hyperParamSelection = selection, K = 50, 
#                            useParallel = T, ncores = cores,
#                            outputPrefix = out_prefix, eps = eps)
if(lambda >=0){
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), maxCausal = 1, nSample = 43568,
                            fittingMethod = fit_method,  hyperParamSelection = selection, K = 100, 
                            pvalueThreshold = 1, rThreshold = 0.8,
                            #useParallel = T, ncores = cores,
                            lambda = lambda,
                            outputPrefix = out_prefix, eps = eps, 
                            annotationIndices = (1:97)[-70])
}else{
    res1 <- caviarbfFineMapping(lociListFile = paste0(dir, "/ga_loci.txt"), maxCausal = 1, nSample = 43568,
                            fittingMethod = fit_method,  hyperParamSelection = selection, K = 100, 
                            pvalueThreshold = 1, rThreshold = 0.8,
                            #useParallel = T, ncores = cores,
                            #lambda = lambda,
                            outputPrefix = out_prefix, eps = eps, 
                            annotationIndices = (1:97)[-70])

}
