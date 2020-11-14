########
# BPSC #
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

if(! require("edgeR")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("edgeR")
}

if(! require("BPSC")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("BPSC")
}
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(BPSC))

##########################
# Fit BPSC To A Dataset #
##########################

fit.BPSC = function(features, 
                    metadata, 
                    libSize, 
                    ID, 
                    transformation,
                    multiple_qvalues) {
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default BPSC model. Use NONE.')
  
  ##########################
  # Standard MAST pipeline #
  ##########################
  
  name_metadata <- names(metadata)
  cpms <- edgeR::cpm(t(features), lib.size = colSums(t(features)) * edgeR::calcNormFactors(t(features)))
  controlIds <- which(metadata[,name_metadata] == levels(factor(metadata[,name_metadata]))[1])
  design <- model.matrix(~ metadata[,name_metadata])
  coef <- 2
  resbp <- BPglm(data = cpms, controlIds = controlIds, design = design, coef = coef, estIntPar = FALSE) 
  fit = resbp
  
  ###################
  # Combine results #
  ###################
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit[,3,1], 'coef')
    pvalMatrix<-get_pval_BPSC(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    paras<-cbind.data.frame(coef.vector, pvalue.vector)
    paras<-paras[, !duplicated(colnames(paras))]
  }
  else{
    coef<- fit$TVAL
    pval<-fit$PVAL
    paras<-cbind.data.frame(coef, pval)
    paras$feature<-names(fit$PVAL)
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
# Fit BPSC To A List of Datasets #
##################################

list.BPSC<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.BPSC","rename.features", "get_pval_BPSC"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "BPSC","edgeR","reshape2"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.BPSC(features, metadata, libSize, ID, transformation, multiple_qvalues)
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
