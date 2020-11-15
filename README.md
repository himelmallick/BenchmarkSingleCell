Large-scale Benchmarking of Single-cell Differential Expression Analysis Methods
================

![Overview](https://github.com/himelmallick/BenchmarkSingleCell/raw/master/cover.png)


## Cotents
- [Overview](#Overview)
- [Description](#Description)
- [Usage](#Usage)
- [Output](#Output)
- [Examples](#Examples)
- [Contributions](#contributions)


## Overview

**BenchmarkSingleCell** repository contains the R scripts for simulating datasets and evaluating methods of differential expression analysis for the simulating datasets. Specifically the folder Library contains the codes for conducting simulations and evaluations: 

(1) **SingleCell_Simulator**: This script contains R function for simulating single cell RNAseq and bulk RNAseq data using simulators such as Splatter.

(2) **SingleCell_Evaluator**: This script contains R function for applying different methods of differential analysis on the simulated datasets and evaluate the performance of the methods based on metrics such as sensitivity, false discovery rate (FDR) and false poitive rate (FPR). 

(3) **allUtilityFunctions_SC**: This script contains R functions for different methods of normalization of scRNAseq and bulk RNAseq data.

(4) **run_BPSC**, ..., **run_zingeR**: These scripts contain R function for different methods of differential expression analysis for single cell data as well as bulk RNAseq data.  


## Description


## Usage

**SingleCell_Simulator**
```

```


**SingleCell_Evaluator**
```

```


## Parameters

### Parameters used for simulating datasets using SC.sim() in SingleCell_Simulator
  
| Parameter                 | Values used in the simulation | Description   |	
| :------------------------ |:-------------:| :-------------|
|ZeroInflate                |True           | logical parameter for zero-inflation in the data|
|RandomEffect               |False          | only when longitudinal design is used|
|metadataType               |'UVB'          | Univariate Binary design matrix
|nSamples                   |c(20,50,100,200,500)| number of cells|
|nPerSubject                | 1             | |
|nGenes                     |c(500,1000,2000)| number of Genes|
|nMetadata                  |1               |number of metadata|
|effectSize                 |c(1,2,2.5,3,3.5)|effect size|
|minFracZeroes              |    0.25        | minimum fraction of zeros|
|pDE                        | 0.1            | proportion of differentially expressed genes|
|nIterations                | 100            |number of times an experiment is repeated|
|rSeed                      |568910          |reproducibility index (random seed)|



### Parameters used for evaluating benchmarking methods using SC.method in SingleCell_Evaluator
  
| Parameter                 | Values       | Description   |	
| :------------------------ |:-------------:| :-------------|
|Methods                | "CPLM","DESeq2","edgeR","MAST","ZICP","zingeR" etc. |all benchmarking methods|
|fitlist (argument for running run_"")   |          |simulated datasets after applying filtering based on pre-specified abundance threshld ad prevalence|
|transfMethod (argument for running run_"")|   default "none"      |  normalzation method   



## Output

### Output for SingleCell_Simulaor
```
List of simulated datasets
```


###  Output for SingleCell_Evaluator

```
Table of marginal p-values, adjusted p-values (q-values after FDR correction), fold change for each feature for every scenario foe every method
```



## Examples
```

```

References
----------

Zappia et al. (2017). [Splatter: Simulation of Single-cell RNA Sequencing Data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0). Genome Biology, 18(1): 1-15.

Citation
--------

Mallick H, Chatterjee S, Chowdhury S, Chatterjee S, Hicks S (2021+). Differential Expression Analysis of Single-cell RNA-Seq Data using Self-adaptive Tweedie Generalized Linear Models (In Submission).


