
#############################################
##### Loading Source data and Libraries #####
#############################################

SC.sim<-function(NUM)
{
reqpkg<-c("data.table","SPsimSeq","pacman",
          "stringr","pkgmaker","optparse","parallel","foreach","MASS",
          "doParallel","future","ROCR","plyr")
for (i in reqpkg){
  print(i)
  print(packageVersion(i))
  suppressPackageStartupMessages(library(i,character.only=TRUE)) 
}
workingDirectory<-'/gpfs/gsfs11/users/chatterjees7/TweedieVerse/SingleCell'
pkgmaker::source_files(paste(workingDirectory,'/Library',sep=''),'*.R')
data<-data.frame(read.delim(file.path(workingDirectory,"Source_Data","Expression.data.Islam.csv"),
                            sep=",",row.names=1),check.names=F)
data<-data[rowSums(data>0)>=5,]

##################################
###### Assigning Parameters ###### 
##################################

ZeroInflate<-TRUE # stands for zero-inflation
RandomEffect<-FALSE # stands for longitudinal design
metadataType<-'UVB' # stands for Univariate Binary design matrix
nSamples<-c(20,50,100,200,500) # stands for number of cells
nPerSubject<-1
nGenes<-c(500,1000,2000) # stands for number of Genes
nMetadata<-1 # stands for number of metadata
effectSize<-c(1,2,2.5,3,3.5) # stands for effect size
minFracZeroes<-0.25
pDE<-.1
nIterations<-100 # stands for how many times an experiment is repeated
rSeed<-568910 # stands for reproducibility index (random seed)
treatment<-ifelse(str_detect(names(data),"ESC"),1,0)
if (ZeroInflate==TRUE){ZeroInflated<-'/ZeroInflated'}
if (ZeroInflate==FALSE){ZeroInflated<-'/noZeroInflated'}
if (RandomEffect==TRUE){RandomEffect<-'RandomEffect'}
if (RandomEffect==FALSE){RandomEffect<-'noRandomEffect'}
inputString<-apply(expand.grid(ZeroInflated,RandomEffect,metadataType, 
                               nSamples,nPerSubject,nGenes,nMetadata,
                               effectSize,minFracZeroes,pDE,
                               nIterations),1,paste0,collapse = "_")
inputStringlabels<-c('ZeroInflated','RandomEffect','metadataType',
                     'nSamples','nPerSubject','nGenes','nMetadata',
                     'effectSize','minFracZeroes',"pDE","nIterations")
inputDirectory <- file.path(workingDirectory,'Input')
if (!dir.exists(inputDirectory))
{
  dir.create(inputDirectory)
}

######################################
##### Simulating ScRNA-Seq data ######
######################################

set.seed(rSeed)
start.time <- Sys.time()
message("Generating Data for Scenario ",NUM)
input<-str_split(inputString[NUM],pattern="_")[[1]]
names(input)<-inputStringlabels
inputStringR<-paste(inputDirectory,paste(inputString[NUM],'.RData',sep=''),sep='')
n.genes<-as.numeric(input["nGenes"])
model.zero.prob<-ZeroInflate
pDE<-as.numeric(input["pDE"])
lfc.thrld<-as.numeric(input["effectSize"])
nSamples<-as.numeric(input["nSamples"])
minFracZeroes<-as.numeric(input["minFracZeroes"])
sub.data<-as.matrix(data[sample(nrow(data),n.genes),])
simlist<-SPsimSeq(n.sim=nIterations,s.data=sub.data,
                  group=treatment,n.genes=n.genes, 
                  batch.config=1,group.config=c(0.5,0.5),
                  pDE=pDE,cand.DE.genes=NULL,lfc.thrld=lfc.thrld,
                  t.thrld=2.5,llStat.thrld=5,tot.samples=nSamples,
                  model.zero.prob=model.zero.prob,genewiseCor=TRUE,
                  log.CPM.transform=TRUE,lib.size.params=NULL,
                  variable.lib.size=FALSE,w=NULL,result.format="list",
                  verbose=FALSE,prior.count=1,return.details=FALSE,
                  n.mean.class=0.2,const.mult=1e6,
                  minFracZeroes=minFracZeroes)
stop.time<- Sys.time()
minutes<-round(difftime(stop.time,start.time,units="min"),1)
message("Data generation time: ",minutes," minutes")
simlist0<-list()
for(j in 1:length(simlist))
{
features<-data.frame(t(simlist[[j]][[1]]),check.names=F)
metadata<-data.frame(Metadata_TP=simlist[[j]][[2]]$Group,check.names=F)
row.names(metadata)<-row.names(features)
libSize<-simlist[[j]][[2]]$sim.Lib.Size
ID<-row.names(features)
simlist[[j]][[3]]$genes<-names(features)
simlist[[j]][[3]]$Truth<-ifelse(simlist[[j]][[3]]$DE.ind=="TRUE",
                                paste(simlist[[j]][[3]]$genes,"TP",sep="_"),
                                simlist[[j]][[3]]$genes)
names(features)<-simlist[[j]][[3]]$Truth
simlist0[[j]]<-list(features=features,metadata=metadata,libSize=libSize,ID=ID)
}
newname<-paste(inputString[NUM],1:nIterations,sep="_")
names(simlist0)<-newname
save(simlist0,file=inputStringR)

##########################
### Filtering Counts=0 ###
##########################

Threshold_Abundance = 0
Threshold_Prevalence = 0.1
simlist.filtered<-list.SimpleFilter(simlist0,Threshold_Abundance,Threshold_Prevalence)
fitlist<-simlist.filtered
plotlist<-fitlist[[sample(1:length(fitlist),1)]]$features

###########################################
#### Checking Retained DE, post filter ####
###########################################

frac.retained<-c()
for(m in 1:length(fitlist))
{
tmp<-fitlist[[m]]$features
frac.retained<-c(frac.retained,ncol(tmp[,grepl("_TP",names(tmp))])/(pDE*n.genes))
}
avg.DE<-round(mean(frac.retained)*100,2)
if(avg.DE<=(pDE*0.5)*n.genes)
{
message("Average Retained DE is critically Low ")
}

###############################
### Checking Zero Inflation ###
###############################

inputplotR<-paste(inputDirectory,paste(inputString[NUM],'.pdf',sep=''),sep='')
fraction.zeros.full<-round(colSums(plotlist==0)/nSamples,1)*100
fraction.zeros.DE<-plotlist[,grepl("_TP",names(plotlist))]
fraction.zeros.DE<-round(colSums(fraction.zeros.DE==0)/nSamples,1)*100
plot.info<-paste("(S=",nSamples,",G=",n.genes,",E=",lfc.thrld,
                 ",DE=",length(fraction.zeros.DE),")",sep="")
pdf(inputplotR)
par(mfrow=c(1,2))
g<-hist(fraction.zeros.full,plot=F)
g$density = g$counts/sum(g$counts)*100
plot(g,freq=FALSE,xlab="% Zero",main="Full Inflation",
     ylab="% Genes",col="steelblue")
h<-hist(fraction.zeros.DE,plot=F)
h$density=h$counts/sum(h$counts)*100
plot(h,freq=FALSE,xlab="% Zero",main="DE Inflation",
     ylab="% Genes",col="steelblue")
mtext(plot.info,side=3,cex=0.8,line=-1,font=2.5,outer=TRUE) 
dev.off()
}




  
  