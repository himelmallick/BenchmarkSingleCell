###########################################################################################
# Zero-inflated Compound Poisson with libray size covariate in the count submodel (ZICP2) #
###########################################################################################

###########################
# Load Essential Packages #
###########################

# pacman, dplyr, tidyverse, fdrtool, ashr, GMPR, swfdr, genefilter, IHW
load_essential_packages()

###########################################
# Load Dedicated Method-specific Packages #
###########################################

pacman::p_load("pbapply", "glmmTMB", "cplm")

if (packageVersion("statmod")!="1.4.33"){
  devtools::install_version("statmod", version = "1.4.33", repos = "http://cran.us.r-project.org")
}

if (packageVersion("cplm")!="0.7-8"){
  devtools::install_version("cplm", version = "0.7-8", repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages(library(statmod))
suppressPackageStartupMessages(library(cplm))


##########################
# Fit ZICP2 To A Dataset #
##########################

fit.ZICP2 <- function(features, 
                     metadata, 
                     libSize, 
                     ID, 
                     transformation,
                     multiple_qvalues){
  
  #########################
  # Transformation if any #
  #########################
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default ZICP2 model. Use NONE.')
  
  #####################
  # Per-feature model #
  #####################
  
  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
    
    ###############################
    # Extract features one by one #
    ###############################
    
    featuresVector <- features[, x]
    
    #################################
    # Create per-feature input data #
    #################################
    
    dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata, libSize, ID)
    dat_sub2<- dat_sub[, !colnames(dat_sub) %in% c('expr', 'ID')]
    formula<-as.formula(paste("expr ~ ", paste(colnames(dat_sub2), collapse= "+")))
    
    #######################
    # Random effect model #
    #######################
    
    if(!length(ID) == length(unique(ID))){ 
      formula<-update(formula, . ~ . + (1|ID))
      fit <- tryCatch({
        fit1 <- glmmTMB::glmmTMB(formula = formula,  
                                 data = dat_sub, 
                                 family = glmmTMB::tweedie(link = "log"), 
                                 ziformula = ~1) # Intercept-only ZI model
      }, error=function(err){
        fit1 <- try({glmmTMB::glmmTMB(formula = formula,  
                                      data = dat_sub, 
                                      family = glmmTMB::tweedie(link = "log"), 
                                      ziformula = ~1)}) # Intercept-only ZI model
        return(fit1)
      })
      
      ###################################
      # Summarize Coefficient Estimates #
      ###################################
      
      if (class(fit) != "try-error"){
        para<-as.data.frame(coef(summary(fit))$cond)[-1,-3]
        para<-para[-nrow(para),] # Remove library size coefficients
      } else{
        print(paste("Fitting problem for feature", x, "returning NA"))
        para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=3))
      }
      colnames(para)<-c('coef', 'stderr', 'pval')
      para$metadata<-colnames(metadata)
      para$feature<-colnames(features)[x]
    }
    
    #######################
    # Fixed effects model #
    #######################
    
    else{ 
      fit <- tryCatch({
        fit1 <- cplm::zcpglm(formula, 
                            data = dat_sub) # Intercept-only ZI model
      }, error=function(err){
        fit1 <- try({cplm::zcpglm(formula, 
                                 data = dat_sub)}) # Intercept-only ZI model
        return(fit1)
      })
      
      ###################################
      # Summarize Coefficient Estimates #
      ###################################
      
      if (class(fit) != "try-error"){
        para<-as.data.frame(summary(fit)$coefficients$tweedie)[-1,-3]
        para<-para[-nrow(para),] # Remove library size coefficients
      } else{
        print(paste("Fitting problem for feature", x, "returning NA"))
        para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=3))  
      }
      colnames(para)<-c('coef', 'stderr', 'pval')
      para$metadata<-colnames(metadata)
      para$feature<-colnames(features)[x]
    }
    return(para)
  })
  
  ###################
  # Combine results #
  ###################
  
  paras<-do.call(rbind, paras)
  
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
# Fit ZICP2 To A List of Datasets #
##################################

list.ZICP2<-function(physeq, transformation = 'NONE', multiple_qvalues = TRUE){
  foreach(physeq = physeq, 
          .export = c("pvalueAdjustment_HM", "append_qvalues",
                      "fit.ZICP2"), 
          .packages = c("tidyverse", "fdrtool", "ashr", "GMPR", "swfdr", "genefilter", "IHW",
                        "pbapply", "cplm", "glmmTMB"),
          .errorhandling = "remove") %dopar% 
    {
      start.time<-Sys.time()
      features<-physeq$features
      metadata<-physeq$metadata
      libSize<-physeq$libSize
      ID<-physeq$ID
      DD<-fit.ZICP2(features, metadata, libSize, ID, transformation, multiple_qvalues)
      DD$pairwiseAssociation<-paste('pairwiseAssociation', 1:nrow(DD), sep='')
      wh.TP<-intersect(grep("[[:print:]]+\\_TP$", DD$metadata), grep("[[:print:]]+\\_TP$", DD$feature))
      newname<-paste0(DD$pairwiseAssociation[wh.TP], "_TP")
      DD$pairwiseAssociation[wh.TP]<-newname
      DD<-dplyr::select(DD, c('pairwiseAssociation', 'feature', 'metadata'), everything())
      stop.time<-Sys.time()
      time<-as.numeric(round(difftime(stop.time, start.time, units="min"),3), units = "mins")
      DD$time<-time
      return(DD)
    }
}

