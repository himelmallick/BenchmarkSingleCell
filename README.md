Large-scale Benchmarking of Single-cell Differential Expression Analysis Methods
================

![Overview](https://github.com/himelmallick/BenchmarkSingleCell/raw/master/cover.png)


## Contents
- [Introduction](#Introduction)
- [Simulation Designs](#SimulationDesigns)
- [Simulation Methods](#SimulationMethods)
- [Examples](#Examples)
- [Output](#Output)
- [Folder Descriptions](#FolderDescriptions)
- [References](#References)
- [Citation](#Citation)



## Introduction
In translational single-cell RNA-seq (scRNA-seq) studies, effective downstream analysis and mechanistic interpretation represents key statistical and computational challenges. While normalization and differential expression analysis have been extensively studied in the bulk RNA-seq experiments that measure average transcript abundances summarized over thousands of cells, the specific characteristics of single-cell RNA-seq data have necessitated the development of dedicated analytical tools (Berge2018, Sekula2019, Finak2015, Wang2019, Korthauer2016, Jaakkola2017). Specifically, several zero-corrected methods in the context of differential expression of scRNA-seq count data have been developed. Nearly all these methods rely on some form of normalization (to account for differences in library size or sequencing depth), and they differ primarily in the choice of the data transformation, regression model, and statistical inference. For example, MAST \citep{Finak2015} employs a two‐part hurdle model consisting of logistic and Gaussian sub-components that analyzes continuous scRNA‐seq expression levels (following a transformation of the raw counts). In a similar vein, scREHurdle (Sekula2019) employs a variational two‐part discrete hurdle model consisting of independent logistic and zero‐truncated negative binomial models. Among other published methods, scDE (Korthauer2016) uses a mixture of a negative binomial distribution and a low‐level Poisson distribution to simultaneously account for overdispersed expression counts from detected transcripts as well as to accommodate genes that are undetected across some of the cells.

Although the commonly used methods described above have been useful to identify differentially expressed (DE) genes between cell populations, they do not accurately incorporate the wide diversity of expression profiles arising from extensively replicated scRNA-seq experiments. We develop a statistical method that addresses the aforementioned shortcomings of previous methods for scRNA-seq DE analysis. We argue that in order to capture the large dynamic range of scRNA-Seq expression profiles, as induced by both heavy tails and zero-inflation, a three-parameter Tweedie family of distributions (Zhang2013, Tweedie1984,  Jorgensen1987) is a more natural choice. Rather than separately modeling the zeroes and non-zeroes, we borrow information across cells by considering all observations together in a single zero-adjusted model. Consequently, by utilizing the cellular sequencing depth as an offset, we bypass the need for an adhoc normalization and/or transformation prior to DE testing. We evaluate the performance of our proposed models through simulations and also performed a large-scale benchmarking by comparing our models to 10 differential expression analysis methods with respect to a variety of performance metrics.



## SimulationDesigns
We simulate data varying a combination of parameters such as number of cells (nSamples), number of genes (nGenes), effect size (effectSize). In the table below we give the list of parameters and their values used in our simulation. 

### Parameters used for simulating datasets using SC.sim() in SingleCell_Simulator
  
| Parameter                 | Values used in the simulation | Description   |	
| :------------------------ |:-------------:| :-------------|
|ZeroInflate                |True           | logical parameter for zero-inflation in the data|
|RandomEffect               |False          | only when longitudinal design is used|
|metadataType               |'UVB'          | Univariate Binary design matrix
|nSamples                   |20,50,100,200,500| number of cells|
|nPerSubject                | 1             | |
|nGenes                     |500,1000,2000| number of Genes|
|nMetadata                  |1               |number of metadata|
|effectSize                 |1,2,2.5,3,3.5|effect size|
|minFracZeroes              |    0.25        | minimum fraction of zeros|
|pDE                        | 0.1            | proportion of differentially expressed genes|
|nIterations                | 100            |number of times an experiment is repeated|
|rSeed                      |568910          |reproducibility index (random seed)|


## Simulation Methods
We performed a large-scale benchmarking by comparing our models to 10 existing differential expression analysis methods. A complete list of tested methods are given below:

```
CPLM
ZICP
BPSC
MAST 
metagenomeSeq
monocle
scDD
scRE
Wilcoxon 
zingeR 
DESeq2 
edgeR
limmaVOOM
```


## Examples

## Output

### Output for SingleCell_Simulaor
```
List of simulated datasets
```


###  Output for SingleCell_Evaluator

```
Table of marginal p-values, adjusted p-values (q-values after FDR correction), fold change for each feature for every scenario foe every method
```



## Folder Descriptions

**BenchmarkSingleCell** repository contains the R scripts for simulating datasets and evaluating methods of differential expression analysis for the simulating datasets. Specifically the folder Library contains the codes for conducting simulations and evaluations: 

(1) **SingleCell_Simulator**: This script contains R function for simulating single cell RNAseq and bulk RNAseq data using simulators such as Splatter.

(2) **SingleCell_Evaluator**: This script contains R function for applying different methods of differential analysis on the simulated datasets and evaluate the performance of the methods based on metrics such as sensitivity, false discovery rate (FDR) and false poitive rate (FPR). 

(3) **allUtilityFunctions_SC**: This script contains R functions for different methods of normalization of scRNAseq and bulk RNAseq data.

(4) **run_BPSC**, ..., **run_zingeR**: These scripts contain R function for different methods of differential expression analysis for single cell data as well as bulk RNAseq data.  


References
----------

Zappia et al. (2017). [Splatter: Simulation of Single-cell RNA Sequencing Data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0). Genome Biology, 18(1): 1-15.

Citation
--------

Mallick H, Chatterjee S, Chowdhury S, Chatterjee S, Hicks S (2021+). Differential Expression Analysis of Single-cell RNA-Seq Data using Self-adaptive Tweedie Generalized Linear Models (In Submission).


