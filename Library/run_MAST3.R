#########
# MAST3 #
#########

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
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

if(! require("MAST")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("MAST")
}
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(MAST))

##########################
# Fit MAST3 To A Dataset #
##########################

fit.MAST3 <-function(features, 
                    metadata, 
                    libSize, 
                    ID, 
                    transformation,
                    multiple_qvalues) {
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default CPLM model. Use NONE.')
  
  ###############################################
  # Standard MAST pipeline (CPM + No Covariate) #
  ###############################################
  
  name_metadata <- names(metadata)
  grp <- metadata[, name_metadata]
  counts = t(features)
  tpms <- counts*1e6/colSums(counts)
  sca <- FromMatrix(exprsArray = log2(tpms + 1), cData = data.frame(wellKey=rownames(features), grp = grp), fData = data.frame(primerid= colnames(features)))
  
  # Differential expression
  zlmdata <- zlm(~ grp, sca = sca)
  fit <- lrTest(zlmdata, "grp")
  
  ###################
  # Combine results #
  ###################
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit[,3,1], 'coef')
    pvalMatrix<-get_pval_MAST(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    paras<-cbind.data.frame(coef.vector, pvalue.vector)
    paras<-paras[, !duplicated(colnames(paras))]  
    }
  else{
    coef<-fit[,3,1]
    pval<-fit[,3,3]
    # coef<-fit[,1,1]
    # pval<-fit[,1,3]
    paras<-cbind.data.frame(coef, pval)
    paras$feature<-rownames(paras)
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
  
  paras<-paras[order(paras$qval_BH, decreasing = FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  rownames(paras)<-NULL
  return(paras)   
}

###################################
# Fit MAST3 To A List of Datasets #
###################################

list.MAST3<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.MAST3","rename.features", "get_pval_MAST"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "MAST","edgeR","reshape2"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.MAST3(features, metadata, libSize, ID, transformation, multiple_qvalues)
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
