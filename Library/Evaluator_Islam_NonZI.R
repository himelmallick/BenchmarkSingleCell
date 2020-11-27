#####################
### Clear Console ###
#####################

rm(list=ls())

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
workDirectory<-'/gpfs/gsfs11/users/chatterjees7/TweedieVerse/SingleCell/Splatter'
pkgmaker::source_files(paste(workDirectory,'/Library',sep=''),'*.R')

#######################
### Calling Methods ###
#######################

Eval_NonZI<-function(NUM,methods)
{
methods.all<-c("CPLM",'CPLM.SCRAN',"DESeq2",
              "DEsingle","edgeR","limmaVOOM",
              "MAST","metagenomeSeq","monocle",
              "scRE","Wilcoxon","ZICP",
              'ZICP.SCRAN',"zingeR")
  
if(str_contains(methods,"all"))
{
methods.selected<-methods.all
}else{
methods.selected<-intersect(methods,methods.all)
}

###################################
### Extracting Input file names ###
###################################

inputDirectory<-paste(workDirectory,'Input_Islam_NonZI',sep='/')
inputlist<-list.files(inputDirectory,pattern=".RData")
inputString<-apply(expand.grid(methods.selected,inputlist),1,paste0,collapse = "_")
  
####################
### Loading Data ###
####################
  
select.string<-inputString[NUM]
inputString.dat<-substring(select.string,regexpr("_",select.string)+1,
                           nchar(select.string))
load(file.path(inputDirectory,inputString.dat))
run.method<-sub("_.*", "",select.string)
print(paste("Scenario",inputString.dat,sep=" "))
print(paste("Method:",run.method,sep=" "))
  
############################
### Parallel Environment ###
############################
  
ncpus<-16
cl<-makeCluster(ncpus)
registerDoParallel(cl)
transfMethod<-"NONE"
  
###########################
### Input Configuration ###
###########################

inputname<-tools::file_path_sans_ext(basename(inputString.dat))
inputStringlabels<-c('Source','ZeroInflated','RandomEffect','metadataType',
                     'nCells','nPerSubject','nGenes','nMetadata',
                     'effectSize',"dropout.mid","pDE","nIterations")
input<-str_split(inputname,pattern="_")[[1]]
names(input)<-inputStringlabels
ZeroInflated<-as.character(input["ZeroInflated"])
RandomEffect<-as.character(input["RandomEffect"])
SourceData<-as.character(input["Source"])
  
############################
### Filtering Low Counts ###
############################
  
Threshold_Abundance = 0
Threshold_Prevalence = 0.1
simlist.filtered<-list.SimpleFilter(simlist0,Threshold_Abundance,Threshold_Prevalence)
fitlist<-simlist.filtered
  
##############################
### Applying Normalization ###
##############################
    
if(run.method=="CPLM.SCRAN"){
normlist<-list.SCRAN(fitlist)
}else if(run.method=="ZICP.SCRAN"){
normlist<-list.SCRAN(fitlist)
}

######################
### Running Method ###
######################
    
if(run.method=="CPLM"){  
output<-list.CPLM(fitlist,transformation=transfMethod)
}else if(run.method=="CPLM.SCRAN"){  
output<-list.CPLM(normlist,transformation=transfMethod)
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
output<-list.metagenomeSeq(fitlist,transformation=transfMethod)
}else if(run.method=="scRE"){
output<-list.scRE(fitlist,transformation=transfMethod)
}else if(run.method=="Wilcoxon"){
output<-list.Wilcoxon(fitlist,transformation=transfMethod)
}else if(run.method=="ZICP"){
output<-list.ZICP(fitlist,transformation=transfMethod)
}else if(run.method=="ZICP.SCRAN"){
output<-list.ZICP(normlist,transformation=transfMethod)
}else if(run.method=="zingeR"){
output<-list.zingeR(fitlist,transformation=transfMethod)
}
if (length(output)<1){
print(paste('No Output for:',run.method,sep=" "))
}
    
###############################
### Shaping Output's format ###
###############################
    
outputString<-paste(run.method,ZeroInflated,RandomEffect,
                    unname(input["metadataType"]),unname(input["nCells"]), 
                    unname(input["nPerSubject"]),unname(input["nGenes"]), 
                    unname(input["nMetadata"]),unname(input["effectSize"]),
                    unname(input["dropout.mid"]),unname(input["pDE"]),
                    unname(input["nIterations"]),sep='_')
    
names(output)<-paste(outputString,1:length(output),sep='_')
dflist<-lapply(output,eval_res_list_repeated)
names(dflist)<-names(output)
simparamslabels<-c("methodName","ZeroInflate","RandomEffect",
                    "metadataType","nCells","nPerSubject", 
                    "nGenes","nMetadata","effectSize",
                    "dropout.mid","pDE","nIterations","rep")
                       
##########################
### Summarizing Output ###
##########################
    
df<-make_power_df(dflist,simparamslabels)
df.summarized<-combine_results(df)
    
#######################
### Printing Output ###
#######################

stopCluster(cl)
outputDirectory<-paste(workDirectory,paste('Outputs',SourceData,'NonZI',sep='_'),sep="/")
if (!dir.exists(outputDirectory))
{
dir.create(outputDirectory)
}
outputStringR<-paste(outputDirectory,paste(SourceData,"_",outputString,'.csv',sep=''),sep='/')
write.table(df.summarized,file=outputStringR,sep=",",row.names=F)
}



