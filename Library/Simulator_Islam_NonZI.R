#####################
### Clear Console ###
#####################

rm(list=ls())

###############################################
##### Initializing Function and Libraries #####
###############################################

Sim_NonZI<-function(NUM)
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
workDirectory<-'/gpfs/gsfs11/users/chatterjees7/TweedieVerse/SingleCell/Splatter'
  
##################################
###### Assigning Parameters ###### 
##################################
  
ZeroInflate<-"FALSE" 
RandomEffect<-"FALSE"
metadataType<-'UVB' 
nCells<-c(100,1000) 
nPerSubject<-1
nGenes<-c(1000,10000) 
nMetadata<-1 
effectSize<-c(0.5,1,2)
dropout.mid<-"None"
SourceData<-"Islam"
ScTech<-"Smart-Seq"
pDE<-0.1
nIterations<-100
seed<-1234

#############################
### Organizing Parameters ###
#############################

if (ZeroInflate=="TRUE"){
  dropout.type<-"experiment" 
  ZeroInflated<-'ZeroInflated'
  Source<-paste('/',SourceData,sep='')
  }
if (ZeroInflate=="FALSE"){
  dropout.type<-"none" 
  ZeroInflated<-'NonZeroInflated'
  Source<-paste('/',SourceData,sep='')
  }

if (RandomEffect=="TRUE"){
  RandomEffect<-'RandomEffect'
  }
if (RandomEffect=="FALSE"){
  RandomEffect<-'noRandomEffect'
  }

######################################
### Getting Parameter Combinations ###
######################################

inputString<-apply(expand.grid(Source,ZeroInflated,RandomEffect,metadataType, 
                               nCells,nPerSubject,nGenes,nMetadata,
                               effectSize,dropout.mid,pDE,nIterations),1,paste0,collapse = "_")
inputStringlabels<-c('Source','ZeroInflated','RandomEffect','metadataType',
                      'nCells','nPerSubject','nGenes','nMetadata',
                      'effectSize',"dropout.mid","pDE","nIterations")

###########################################
##### Setting Simulation Environment ######
###########################################
  
print(paste("Generating Data for Scenario:",NUM,sep=' '))
input<-str_split(inputString[NUM],pattern="_")[[1]]
names(input)<-inputStringlabels
nGenes<-as.numeric(input["nGenes"])
nCells<-as.numeric(input["nCells"])
effectSize<-as.numeric(input["effectSize"])
dropout.mid<-as.character(input["dropout.mid"])
pDE<-as.numeric(input["pDE"])

#########################
### Input Source Data ###
#########################

load(paste(workDirectory,"Data","Islam.RData",sep="/"))

#########################################
### Estimating parameters Source Data ###
#########################################

params<-splatEstimate(sce)

##################################
### Generating Hyperparameters ###
##################################

sim<-list()
for(i in 1:nIterations){
  params<-setParams(params,
                    nGenes=nGenes,
                    seed=seed+i,
                    de.facLoc=effectSize,
                    batchCells=nCells,
                    group.prob=c(0.5,0.5),
                    dropout.type=dropout.type,
                    de.prob=pDE)
  sim[[i]]<-splatSimulate(params,method="groups",verbose=F)
  }
stop.time<- Sys.time()
minutes<-round(difftime(stop.time,start.time,units="min"),1)
print(paste("Data generation time: ",minutes," minutes",sep=' '))

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
simlist0[[j]]<-list(features=features,metadata=metadata,libSize=libSize,ID=ID)
}

###############################
### Labeling Simulated Data ###
###############################

newname<-paste(inputString[NUM],1:nIterations,sep="_")
names(simlist0)<-newname

################################
### Setting Output Directory ###
################################

OutDirectory<-file.path(workDirectory,paste('Input',SourceData,'NonZI',sep='_'))
if (!dir.exists(OutDirectory))
  {
  dir.create(OutDirectory)
  }
OutStringR<-paste(OutDirectory,paste(inputString[NUM],'.RData',sep=''),sep='')
OutStringC<-paste(OutDirectory,paste(inputString[NUM],'.csv',sep=''),sep='')

#############################
### Saving Simulated Data ###
#############################

save(simlist0,file=OutStringR)

###########################################
### Checking Per Feature Zero Inflation ###
###########################################

inputplotR<-paste(OutDirectory,paste(inputString[NUM],'.pdf',sep=''),sep='')
pdf(inputplotR,width=7,height=5)
par(mfrow=c(2,2))
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

plot.info<-paste("(Cells = ",nCells,", Genes = ",nGenes,", Added Zero = ","No",
                 ", Effect Size = ",effectSize,", Source = ",ScTech,")",sep="")

c<-hist(obsZeroes$Obs.Zero,plot=F)
plot(c,ylab="Number of Genes",main="Zero-inflation (All Genes)",
     xlab="Fraction Zero",col="steelblue")

d<-hist(DEZeroes$DE.Zero,plot=F)
plot(d,ylab="Number of Genes",main="Zero-inflation (DE Genes)",
     xlab="Fraction Zero",col="steelblue")

e<-hist(DE1,plot=F)
plot(e,ylab="Number of Genes",main="Expression Cell Type 1 (DE Genes)",
     xlab="Log Mean Expression",col="steelblue")

f<-hist(DE2,plot=F)
plot(f,ylab="Number of Genes",main="Expression Cell Type 2 (DE Genes)",
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
                   dropout.Rate=dropout.mid,
                   DE=rep(nrow(DEZeroes),nrow(featureRS)))
fwrite(SimTab,OutStringC)                                 
}





