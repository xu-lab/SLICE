# SLICE

SLICE is an algorithm that utilizes single-cell RNA-seq (scRNA-seq) data to quantitatively measure cellular differentiation states based on single cell entropy and predict cell differentiation lineages via the construction of entropy directed cell trajectories.

## Installation

* In R or RStudio, type the following command to install devtools
  
  ```
  install.packages("devtools")
  library(devtools)
  ```
  
* Then, use devtools to install SINCERA from github
  
  ```
  devtools::install_github("xu-lab/SLICE")
  ```

* Use library() to activate SINCERA

  ```
  library(SLICE)
  ```

## Demonstration

* A demonstration of using SLICE to reconstruct a two-branched lung fibroblast differentiation lineage from E16.5 mouse lung single cell data can be found at https://github.com/minzheguo/SLICE/blob/master/demo/FB.R. 


## Citation


* Minzhe Guo, Erik L. Bao, Michael Wagner, Jeffrey A. Whitsett, Yan Xu. 2016. SLICE: determing cell differentiation and lineage based on single cell entropy. Nucleic Acids Research. doi:10.1093/nar/gkw1278. (MG and ELB are co-first authors) 
