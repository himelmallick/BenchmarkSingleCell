
#########################
### Loading Libraries ###
#########################

reqpkg<-c("data.table","pacman","stringr",
          "pkgmaker","optparse","sjmisc",
          "parallel","foreach","MASS",
          "doParallel","future","ROCR","plyr")
for (i in reqpkg){
print(i)
print(packageVersion(i))
suppressPackageStartupMessages(library(i,character.only=TRUE)) 
}
workingDirectory<-'/gpfs/gsfs11/users/chatterjees7/TweedieVerse/SingleCell'
pkgmaker::source_files(paste(workingDirectory,'/Library',sep=''),'*.R')

#######################
### Calling Methods ###
#######################

SC.method<-function(NUM,methods)
{
methods.all<-c("BPSC","CPLM","DESeq2","DEsingle","edgeR","limmaVOOM",
               "MAST","metagenomeSeq","monocle","negbin",
               "scDD","Wilcoxon","ZIB","ZICP","ZINB","zingeR")

if(str_contains(methods,"all"))
{
methods.selected<-methods.all
}else{
methods.selected<-intersect(methods,methods.all)
}
inputDirectory<-'/gpfs/gsfs11/users/chatterjees7/TweedieVerse/SingleCell/Input'
inputlist<-list.files(inputDirectory,pattern=".RData")
inputString<-apply(expand.grid(methods.selected,inputlist),1,paste0,collapse = "_")
select.string<-inputString[NUM]
inputString.dat<-substring(select.string,regexpr("_",select.string)+1,
                           nchar(select.string))
load(file.path(inputDirectory,inputString.dat))
run.method<-sub("_.*", "",select.string)

###################################
### Gathering data & Parameters ###
###################################

ncpus<-48
cl<-makeCluster(ncpus)
registerDoParallel(cl)
transfMethod<-"NONE"
print(paste(run.method,"scenario",NUM,"running",sep=" "))
inputname<-tools::file_path_sans_ext(basename(inputString.dat))
inputStringlabels<-c('ZeroInflated','RandomEffect','metadataType',
                     'nSamples','nPerSubject','nGenes','nMetadata',
                     'effectSize','minFracZeroes',"pDE","reps")
input<-str_split(inputname,pattern="_")[[1]]
names(input)<-inputStringlabels
ZeroInflated<-as.character(input["ZeroInflated"])
RandomEffect<-as.character(input["RandomEffect"])

##########################
### Filtering Counts=0 ###
##########################

Threshold_Abundance = 0
Threshold_Prevalence = 0.1
simlist.filtered<-list.SimpleFilter(simlist0,Threshold_Abundance,Threshold_Prevalence)
fitlist<-simlist.filtered

#########################################
### Checking Retained DE, post filter ###
#########################################

n.genes<-as.numeric(input["nGenes"])
pDE<-as.numeric(input["pDE"])
frac.retained<-c()
for(m in 1:length(fitlist))
{
tmp<-fitlist[[m]]$features
frac.retained<-c(frac.retained,ncol(tmp[,grepl("_TP",names(tmp))])/(pDE*n.genes))
}
avg.DE<-round(mean(frac.retained)*100,2)
if(avg.DE<=(pDE*0.5)*n.genes)
{
message("Average Retained DE is critically Low")
}

######################
### Running Method ###
######################

if(run.method=="BPSC"){  
output<-list.BPSC(fitlist,transformation=transfMethod)
}else if(run.method=="CPLM"){  
output<-list.CPLM(fitlist,transformation=transfMethod)
}else if(run.method=="DESeq2"){  
output<-list.DESeq2(fitlist,transformation=transfMethod)
}else if(run.method=="DEsingle"){  
output<-list.DEsingle(fitlist,transformation=transfMethod)
}else if(run.method=="edgeR"){
output<-list.edgeR(fitlist,transformation=transfMethod)
}else if(run.method=="limmaVOOM"){
output<-list.limmaVOOM(fitlist,transformation=transfMethod)
}else if(run.method=="MAST"){
output<-list.MAST(fitlist,transformation=transfMethod)
}else if(run.method=="metagenomeSeq"){
output<-list.metagenomeSeq(fitlist,transformation=transfMethod)
}else if(run.method=="monocle"){
output<-list.monocle(fitlist,transformation=transfMethod)
}else if(run.method=="negbin"){
output<-list.negbin(fitlist,transformation=transfMethod)
}else if(run.method=="scDD"){
output<-list.scDD(fitlist,transformation=transfMethod)
}else if(run.method=="Wilcoxon"){
output<-list.Wilcoxon(fitlist,transformation=transfMethod)
}else if(run.method=="ZIB"){
output<-list.ZIB(fitlist,transformation=transfMethod)
}else if(run.method=="ZICP"){
output<-list.ZICP(fitlist,transformation=transfMethod)
}else if(run.method=="ZINB"){
output<-list.ZINB(fitlist,transformation=transfMethod)
}else if(run.method=="zingeR"){
output<-list.zingeR(fitlist,transformation=transfMethod)
}
if (length(output)<1){
print(paste('Consistent error in the model fitting. No output returned for',run.method,sep=""))
}
outputString<-paste(run.method,ZeroInflated,RandomEffect,
                    unname(input["metadataType"]),unname(input["nSamples"]), 
                    unname(input["nPerSubject"]),unname(input["nGenes"]), 
                    unname(input["nMetadata"]),unname(input["effectSize"]), 
                    unname(input["minFracZeroes"]),unname(input["pDE"]),sep='_')
                   
names(output)<-paste(outputString,1:length(output),sep='_')
dflist<-lapply(output,eval_res_list_repeated)
names(dflist)<-names(output)
simparamslabels<-c("methodName","ZeroInflate","RandomEffect",
                   "metadataType","nSamples","nPerSubject", 
                   "nGenes","nMetadata","effectSize",'minFracZeroes',
                   "pDE","rep")
df<-make_power_df(dflist,simparamslabels)

stopCluster(cl)
outputDirectory<-file.path(workingDirectory,'Outputs')
if (!dir.exists(outputDirectory))
{
dir.create(outputDirectory)
}
outputStringR<-paste(outputDirectory,paste(outputString,'.csv',sep=''),sep='/')
write.table(df.averaged,file=outputStringR,sep=",",row.names=F)
}

