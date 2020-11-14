########
# scDD #
########

###########################
# Load Essential Packages #
###########################

# pacman, devtools, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load('reshape2')

if(! require("SummarizedExperiment")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
}

if(! require("scran")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("scran")
}

if(! require("scDD")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("scDD")
}

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scDD))


#########################
# Fit scDD To A Dataset #
#########################

fit.scDD = function(features, 
                    metadata, 
                    libSize, 
                    ID, 
                    transformation,
                    multiple_qvalues) {
  
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default scDD model. Use NONE.')
  
  ##########################
  # Standard scDD pipeline #
  ##########################
  
  scDatList <- list()
  name_metadata <- names(metadata)
  groups <- unique(metadata[,name_metadata])
  for (i in 1:length(groups)) {
    scDatList[[paste0("G", i)]] <- as.matrix(metadata[metadata$celltype == groups[i],])
  }
  sce <- SingleCellExperiment(assays=list(normcounts=t(features)), colData=metadata)
  datNorm.scran <- scDD::preprocess(sce, condition = name_metadata,zero.thresh = 0.75, median_norm = F)
  condition <- metadata
  condition <- as.numeric(as.factor(as.matrix(condition)))
  names(condition) <- colnames(datNorm.scran)
  SDSumExp <- SingleCellExperiment(assays = list(normcounts = datNorm.scran@assays@data@listData$normcounts), colData = data.frame(condition))
  prior_param <- list(alpha = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01)
  scd <- scDD(SDSumExp, prior_param = prior_param, testZeroes = FALSE, param = BiocParallel::MulticoreParam(workers = 1), condition = "condition", min.size = 3, min.nonzero = NULL)
  fit <- results(scd)
  
  ###################
  # Combine results #
  ###################
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit[,3,1], 'coef')
    pvalMatrix<-get_pval_scDD(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    paras<-cbind.data.frame(coef.vector, pvalue.vector)
    paras<-paras[, !duplicated(colnames(paras))]
  }
  else{
    coef<- rep(NA, dim(fit)[1])
    pval<-fit$nonzero.pvalue
    paras<-cbind.data.frame(coef,pval)
    paras$feature<-fit$gene
    paras$metadata<- names(metadata)
  }
  
  ###############################################
  # Calculate multiple qvalues only if prompted #
  ###############################################
  
  if(multiple_qvalues){
    paras<-append_qvalues(features, metadata, paras)
  } else{
    paras$qval_BH<-as.numeric(p.adjust(paras$pval, method = 'BH'))
  }
  
  #################
  # Return output #
  #################
  
  paras<-paras[order(paras$qval_BH, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  rownames(paras)<-NULL
  return(paras)   
}

##################################
# Fit scDD To A List of Datasets #
##################################

list.scDD<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.scDD","rename.features", "get_pval_scDD"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "SummarizedExperiment", "scran","scDD", "reshape2"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.scDD(features, metadata, libSize, ID, transformation, multiple_qvalues)
      DD$pairwiseAssociation<-paste('pairwiseAssociation', 1:nrow(DD), sep='')
      wh.TP<-intersect(grep("[[:print:]]+\\_TP$", DD$metadata), grep("[[:print:]]+\\_TP$", DD$feature))
      newname<-paste0(DD$pairwiseAssociation[wh.TP], "_TP")
      DD$pairwiseAssociation[wh.TP]<-newname
      DD<-dplyr::select(DD, c('pairwiseAssociation', 'feature', 'metadata'), everything())
      stop.time<-Sys.time()
      time<-as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
      DD$time<-time
      return(DD)
    }
}
