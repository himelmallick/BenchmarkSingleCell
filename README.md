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


## Arguments

### Arguments for SingleCell_Simulator
  
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|




### Arguments for SingleCell_Evaluator
  
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|



## Output

### Output for SingleCell_Simulaor

a list of three components

| Object       | Description   |
| :------------------------ | :-------------|


###  Output for SingleCell_Evaluator

a list of three components

| Object       | Description   |
| :------------------------ | :-------------|



## Examples
```

```

References
----------

Zappia et al. (2017). [Splatter: Simulation of Single-cell RNA Sequencing Data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0). Genome Biology, 18(1): 1-15.

Citation
--------

Mallick H, Chatterjee S, Chowdhury S, Chatterjee S, Hicks S (2021+). Differential Expression Analysis of Single-cell RNA-Seq Data using Self-adaptive Tweedie Generalized Linear Models (In Submission).


