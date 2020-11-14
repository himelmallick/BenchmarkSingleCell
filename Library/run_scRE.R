#####################
# scREhurdle (scRE) #
#####################

###########################
# Load Essential Packages #
###########################

# pacman, devtools, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

if(! require("scREhurdle")) {
  devtools::install_github("mnsekula/scREhurdle")
}

suppressPackageStartupMessages(library(scREhurdle))

#########################
# Fit scRE To A Dataset #
#########################

fit.scRE <-function(features, 
                    metadata, 
                    libSize, 
                    ID, 
                    transformation,
                    multiple_qvalues) {
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default CPLM model. Use NONE.')
  
  ##########################
  # Standard scRE pipeline #
  ##########################
  
  name_metadata <- names(metadata)
  grp <- metadata[,name_metadata]
  scDat <- as.data.frame(t(features))
  mod.IRE <- scREhurdle(Y=scDat, treatGroup=grp, typeRE = "ind", stan_seed=523)
  res <- mod.IRE$deTab
  
  ###################
  # Combine results #
  ###################
  
  paras<-cbind.data.frame(coef = res$chisq, pval = res$chisq.pval)
  # paras<-cbind.data.frame(coef = res$C.Z, pval = res$C.pval)
  paras$feature<-rownames(res)
  paras$metadata<- names(metadata)

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

##################################
# Fit scRE To A List of Datasets #
##################################

list.scRE<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.scRE"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "scREhurdle"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.scRE(features, metadata, libSize, ID, transformation, multiple_qvalues)
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
