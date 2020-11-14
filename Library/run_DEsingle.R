############
# DEsingle #
############

###########################
# Load Essential Packages #
###########################

# pacman, devtools, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load('reshape2')

if(! require("edgeR")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("edgeR")
}

if(! require("DEsingle")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DEsingle")
}
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DEsingle))


##########################
# Fit MAST To A Dataset #
##########################

fit.DEsingle = function(features, 
                        metadata, 
                        libSize, 
                        ID, 
                        transformation,
                        multiple_qvalues) {
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default DEsingle model. Use NONE.')
  
  ##############################
  # Standard DEsingle pipeline #
  ##############################
  
  name_metadata <- names(metadata)
  results <- DEsingle(counts = round(t(features)), group = factor(metadata[,name_metadata]))
  results <- data.frame(results, stringsAsFactors = FALSE)
  fit <- results
  
  ###################
  # Combine results #
  ###################
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit[,3,1], 'coef')
    pvalMatrix<-get_pval_DEsingle(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    paras<-cbind.data.frame(coef.vector, pvalue.vector)
    paras<-paras[, !duplicated(colnames(paras))]
  }
  else{
    coef<- fit$foldChange
    pval<-fit$pvalue
    paras<-cbind.data.frame(coef,pval)
    paras$feature<-rownames(fit)
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


######################################
# Fit DEsingle To A List of Datasets #
######################################

list.DEsingle<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.DEsingle","rename.features", "get_pval_DEsingle"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "DEsingle","edgeR","reshape2"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.DEsingle(features, metadata, libSize, ID, transformation, multiple_qvalues)
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


