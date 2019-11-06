# Consistency for Stochastic Block Model

This package works on consistency for stochastic block model.

# Installation

```
devtools::install_github("tbilab/bisbmSC")

library(tidyverse)
library(furrr)
library(future)
library(Matrix)
library(dplyr)
library(doParallel)
library(foreach)
library(pheSBMR)
library(aws.s3)
library(igraph)
library(RSpectra)
library(expm)
library(cluster)
```

# Example

1.Iteratively run stochastic block model multiple times and cluster the subjects (nodes) to generate an ensemble of partitions; Will return a list of dataframes

```
#5 is the number of cores
#dat is the long format dataframe with two columns of a nodes and b nodes. Every instance of b code occuring gets a row.
#1000 is the times of running sbm

sbm_sim <- cleanSBM(1000,5,dat)
```

2.Create the affinity matrix(similarity matrix) for each partition, whose element aij represents whether the subjects (nodes) i and j are assigned to the same community for each partition

```
start <- Sys.time()
for (i in 1:2000) {
  s3saveRDS(create_affinity(i),bucket = "/affinitymatrix",object = paste("sbmmatrix_",i,".rds",sep = ""))
}
end <- Sys.time()
```

3.Randomly choose a number of affinity matrices and computed a consensus graph represented by an aggregated probability affinity matrix P, whose element pij counts probability the subjects (nodes) i and j are assigned to the same community across the full ensemble. Then cluster based on different approaches.(This step is used for comparing different approaches. If you are confident on one approach and have an aggregated probability affinity matrix P, you can directly skip this step and use `build_tree` function to get a top-down recursive partition to learn the hierarchical structure) 

```
plan(multiprocess,workers=10)
aug_result <- future.apply::future_lapply(seq(20,1000,20), function(x){
  library(aws.s3)
  library(igraph)
  library(RSpectra)
  library(expm)
  library(cluster)
  augmentation_matrix(x,2000,339,c(1:10),10,"all","/affinitymatrix",0.1,0.95,level4.R)}
)
```

4.Input the aggregated probability affinity matrix P and take a top-down recursive partition approach to learn the hierarchical structure of the subjects(nodes).

```
cluster_result <- build_tree(A, xi.loc.labels=list(), ncl=0, cl.labels=1:339, n.min=25, D=NULL)
```

5.To get the community memberships matrix and visualize the cluster result. Will return two dataframes of plot data and membership dataframe.

```
membership(cluster_result)
```