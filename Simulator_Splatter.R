
#############################################
##### Loading Source data and Libraries #####
#############################################

SC.sim<-function(NUM)
{
start.time <- Sys.time()
reqpkg<-c("data.table","splatter","pacman",
          "stringr","pkgmaker","optparse","parallel","foreach","MASS",
          "doParallel","future","ROCR","plyr")
for (i in reqpkg){
print(i)
print(packageVersion(i))
suppressPackageStartupMessages(library(i,character.only=TRUE)) 
}
workingDirectory<-'/gpfs/gsfs11/users/chatterjees7/TweedieVerse/SingleCell/Splatter'
pkgmaker::source_files(paste(workingDirectory,'/Library',sep=''),'*.R')
  
##################################
###### Assigning Parameters ###### 
##################################
  
ZeroInflate<-c("ZI","NoZI") 
RandomEffect<-FALSE 
metadataType<-'UVB' 
nCells<-c(100,1000) 
nPerSubject<-1
nGenes<-c(1000,5000) 
nMetadata<-1 
effectSize<-c(0.5,1,2) 
pDE<-0.1
sourceData<-c("Islam","PBMC")
nIterations<-100
seed<-1234

#########################
### Organizing Inputs ###
#########################

if (ZeroInflate=="ZI"){
  dropout.type<-"experiment" 
  ZeroInflated<-'/ZeroInflated'}
if (ZeroInflate=="NoZI"){
  dropout.type<-"none" 
  ZeroInflated<-'/noZeroInflated'}
if (RandomEffect==TRUE){RandomEffect<-'RandomEffect'}
if (RandomEffect==FALSE){RandomEffect<-'noRandomEffect'}
inputString<-apply(expand.grid(sourceData,ZeroInflated,RandomEffect,metadataType, 
                               nCells,nPerSubject,nGenes,nMetadata,
                               effectSize,pDE,nIterations),1,paste0,collapse = "_")
inputStringlabels<-c('sourceData','ZeroInflated','RandomEffect','metadataType',
                       'nCells','nPerSubject','nGenes','nMetadata',
                       'effectSize',"pDE","nIterations")
inputDirectory <- file.path(workingDirectory,'Input')
if (!dir.exists(inputDirectory))
  {
  dir.create(inputDirectory)
}

if(sourceData=="Islam"){
  load(paste(workingDirectory,"Data","Islam.RData",sep="/"))
} else if(sourceData=="PBMC"){
  load(paste(workingDirectory,"Data","PBMC.RData",sep="/"))
} else {
  print(paste("Data Mismatch"))
  stop()
}
  
###########################################
##### Setting Simulation Environment ######
###########################################
  
message("Generating Data for Scenario",NUM)
input<-str_split(inputString[NUM],pattern="_")[[1]]
names(input)<-inputStringlabels
inputStringR<-paste(inputDirectory,paste(inputString[NUM],'.RData',sep=''),sep='')
inputStringC<-paste(inputDirectory,paste(inputString[NUM],'.csv',sep=''),sep='')
sourceData<-as.numeric(input["sourceData"])
nGenes<-as.numeric(input["nGenes"])
nCells<-as.numeric(input["nCells"])
effectSize<-as.numeric(input["effectSize"])
pDE<-as.numeric(input["pDE"])

#########################################
### Estimating parameters Source Data ###
#########################################

params<-splatEstimate(sce)
dropout<-assay(sim[[j]],"")

##################################
### Generating Hyperparameters ###
##################################

sim<-list()
for(i in 1:nIterations){
params<-setParams(params,update=list(nGenes=nGenes,
                                     seed=seed+i,
                                     de.facLoc=effectSize,
                                     batchCells=nCells,
                                     group.prob=c(0.5,0.5),
                                     de.prob=pDE,
                                     .type=.type))
sim[[i]]<-splatSimulate(params,method="groups",verbose=F)
}
stop.time<- Sys.time()
minutes<-round(difftime(stop.time,start.time,units="min"),1)
message("Data generation time: ",minutes," minutes")

######################
### Gathering Data ###
######################

simlist0<-list()
for(j in 1:length(sim))
{
colInfo<-data.frame(colData(sim[[j]]))
rowInfo<-data.frame(rowData(sim[[j]]))
features<-assay(sim[[j]],"counts")
features<-data.frame(t(features))
trueCounts<-assay(sim[[j]],"TrueCounts")

############################
### Summarizing Metadata ###
############################

colInfo$Meta<-ifelse(colInfo$Group=="Group1",0,1)
metadata<-data.frame(Metadata_TP=colInfo$Meta,row.names=rownames(colInfo))

############################
### Summarizing features ###
############################

rowInfo$TP<-ifelse(rowInfo$DEFacGroup1==1 & rowInfo$DEFacGroup2==1,FALSE,TRUE)
rownames(rowInfo)<-ifelse(rowInfo$TP,paste(rownames(rowInfo),
                                           "TP",sep="_"),rownames(rowInfo))
names(features)<-rownames(rowInfo)
rownames(features)<-rownames(metadata)

###############################
### Summarizing LibSize, ID ###
###############################

libSize<-rowSums(features)
ID<-rownames(features)
simlist0[[j]]<-list(features=features,metadata=metadata,libSize=libSize,
                    ID=ID,trueCounts=trueCounts)
}

###############################
### Labeling Simulated Data ###
###############################

newname<-paste(inputString[NUM],1:nIterations,sep="_")
names(simlist0)<-newname

#############################
### Saving Simulated Data ###
#############################

save(simlist0,file=inputStringR)

###########################################
### Checking Per Feature Zero Inflation ###
###########################################

inputplotR<-paste(inputDirectory,paste(inputString[NUM],'.pdf',sep=''),sep='')
pdf(inputplotR)
par(mfrow=c(2,3))
simulation<-simlist0[[1]]
metadata<-simulation$metadata
metadata$Cells<-rownames(metadata)

########################################
### Zero-Inflation in simulated Data ###
########################################

featureRS<-data.frame(t(simulation$features))
obsZeroes<-unname(rowSums(featureRS==0))
obsZeroes<-data.frame(Gene=rownames(featureRS),Obs.Zero=obsZeroes/dim(featureRS)[2])

##################################
### Zero-Inflation in DE Genes ###
##################################

DEfeatureRS<-featureRS[grepl("_TP",rownames(featureRS)),]
DEZeroes<-unname(rowSums(DEfeatureRS==0))
DEZeroes<-data.frame(Gene=rownames(DEfeatureRS),DE.Zero=DEZeroes/dim(DEfeatureRS)[2])

###################################
### True Zero in simulated Data ###
###################################

trueZero<-unname(rowSums(simulation$trueCounts==0))
trueZero<-data.frame(Gene=rownames(simulation$trueCounts),
                     True.Zero=trueZero/dim(simulation$trueCounts)[2])

####################################
### False Zero in simulated Data ###
####################################

falseZero<-data.frame(Gene=rownames(featureRS),
                      False.Zero=abs(obsZeroes$Obs.Zero-trueZero$True.Zero))

######################################
### Mean Expression in DE genes G1 ###
######################################

cells<-metadata[metadata$Metadata_TP==1,]
LmeanExprs1<-featureRS[,names(featureRS)%in%cells$Cells]
DE1<-LmeanExprs1[grepl("_TP",rownames(LmeanExprs1)),]
DE1<-log(as.numeric(apply(DE1,1,mean))+0.0001)

######################################
### Mean Expression in DE genes G2 ###
######################################

cells<-metadata[metadata$Metadata_TP==0,]
LmeanExprs2<-featureRS[,names(featureRS)%in%cells$Cells]
DE2<-LmeanExprs2[grepl("_TP",rownames(LmeanExprs2)),]
DE2<-log(as.numeric(apply(DE2,1,mean))+0.00001)

###############################
### Plotting Zero-Inflation ###
###############################

plot.info<-paste("(Cells = ",nCells,", Genes = ",nGenes,",  = ",ZeroInflate,
                 ", Effect Size = ",effectSize,")",sep="")

a<-hist(trueZero$True.Zero,breaks=,plot=F)
plot(a,ylab="Number of Genes",main="Biological Zero",
     xlab="Fraction Zero",col="steelblue")

b<-hist(falseZero$False.Zero,plot=F)
plot(b,ylab="Number of Genes",main=" Zero",
     xlab="Fraction Zero",col="steelblue")

c<-hist(obsZeroes$Obs.Zero,plot=F)
plot(c,ylab="Number of Genes",main="Zero-inflation (All Genes)",
     xlab="Fraction Zero",col="steelblue")

d<-hist(DEZeroes$DE.Zero,plot=F)
plot(d,ylab="Number of Genes",main="Zero-inflation (DE Genes)",
     xlab="Fraction Zero",col="steelblue")

e<-hist(DE1,plot=F)
plot(e,ylab="Number of Genes",main="Average Expression in DE Genes (Cell Type 1)",
     xlab="Log Mean Expression",col="steelblue")

f<-hist(DE2,plot=F)
plot(f,ylab="Number of Genes",main="Average Expression in DE Genes (Cell Type 2)",
     xlab="Log Mean Expression",col="steelblue")
mtext(plot.info,side=3,cex=0.8,line=-1,font=2.5,outer=TRUE) 
dev.off()

####################################
### Simulation Descriptive table ###
####################################

SimTab<-data.frame(Gene=rownames(featureRS),
                   EffectSize=rep(effectSize,nrow(featureRS)),
                   nCells=rep(nCells,nrow(featureRS)),
                   Obs.Zer=obsZeroes$Obs.Zero,
                   True.Zero=trueZero$True.Zero,
                   False.Zero=falseZero$False.Zero,
                   DE=rep(nrow(DEZeroes),nrow(featureRS)))
fwrite(SimTab,inputStringC)                                 
}





