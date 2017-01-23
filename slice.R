# SLICE: Single Cell Lineage Inference Using Cell Expression Similarity and Entropy

# ==============

# Copyright (C) 2016  Minzhe Guo (minzhe.guo@cchmc.org)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# ================

# Author: Minzhe Guo (minzhe.guo@cchmc.org)
# Version: a12242016

# This is the R script of the SLICE algorithm for determing cell differentiation and lineage based on single cell entropy. 
# To use this script, you will need the R statistical computing environment (version 3.2.5 or later) and several packages freely available through Bioconductor and CRAN, including
# * Bioconductor::Biobase, R::ggplot2, R::igraph, R::reshape2, R::entropy, R::cluster, Bioconductor::graph, Bioconductor::BioNet, R::princurve, R: lmtest, R::mgcv
# SLICE will try to resolve dependencies automatically. If the dependencies cannot be resolved, please refer to the website of each package for more information.

# SLICE is under active development. Core features have been implemented. We are improving the documentation and the visualization functions, and refining the user interfaces.
# Updates of SLICE will be distributed primarily through the SLICE website at: http://research.cchmc.org/pbge/slice.html.

# If you publish results obtained using SLICE, please cite
#   Guo M, Bao EL, Wagner M, Whitsett JA, Xu Y. 2016. SLICE: determing cell differentiation and lineage based on single cell entropy. Nucleic Acids Research. doi:10.1093/nar/gkw1278.

# ================

#######################################################
#                 Resolving Dependencies              #
#######################################################

if(!require(Biobase)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biobase")
}
if(!require(ggplot2)) {
    install.packages('ggplot2', dep=T)
}
if (!require(igraph)) {
    install.packages('igraph', dep=T)
}
if (!require(reshape2)) {
    install.packages('reshape2', dep=T)
}
if(!require(entropy)){
	install.packages('entropy', dep=T)
}
if(!require(graph)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("graph")
}
if(!require(BioNet)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("BioNet")
}
if(!require(cluster)){
	install.packages('cluster', dep=T)
}
if(!require(princurve)){
    install.packages('princurve', dep=T)
}
if(!require(lmtest)){
  install.packages('lmtest', dep=T)
}
if(!require(mgcv)){
  install.packages('mgcv', dep=T)
}
if(!require(gridExtra)){
  install.packages('gridExtra', dep=T)
}



require(gridExtra)
require(Biobase)
require(graph)
require(BioNet)
require(entropy)
require(cluster)
require(ggplot2)
require(grid)
require(princurve)
require(splines)
require(mgcv)
require(lmtest)
require(igraph)





###########################################################################################################
####################                                                                   ####################
####################                slice class                                        ####################
####################                                                                   ####################
###########################################################################################################



#' The slice Class
#'
#' The object storing all information associated with the slice analysis, including data, annotations, analyes, etc. 
#'
#' Information is stored in slots. Key slots include:
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{data}:}{\code{"ExpressionSet"}, the scRNA-seq data in ExpressionSet format, encoding the expression matrix, and cell and gene information }
#'    \item{\code{ident}:}{\code{"vector"}, the identity (cluster/state/type) of cells  }
#'    \item{\code{ident.slice}:}{\code{"vector"}, the identity (cluster/state/type) of single cells inferred by SLICE  }
#'    \item{\code{ident.original}:}{\code{"vector"}, the original identity (cluster/state/type) of single cells provided by the user during object construction}
#'    \item{\code{entropy}:}{\code{"vector"}, the entropy of single cells  }
#'    \item{\code{vds}:}{\code{"data.frame"}, a data frame containing the x and y coordinates of cells for visualization  }
#'    \item{\code{rds}:}{\code{"data.frame"}, a data frame containing reduced dimensions of each cell   }    
#'    \item{\code{model}:}{\code{"list"}, a list containing slice inferred lineage model   }    
#'    \item{\code{transitions}:}{\code{"list"}, a list containing slice inferred cell trajectories   }    
#'    \item{\code{profiles}:}{\code{"list"}, a list containing slice inferred cell trajectory dependent gene expression profiles   }    
#'    \item{\code{projname}:}{\code{"character"}, a string describing the analysis  }    
#'      
#'}
#' @name slice
#' @rdname slice
#' @aliases slice-class
#' @exportClass slice

slice <- setClass("slice", slots = 
                     c(cellnames="vector", genenames="vector", exprmatrix="data.frame",
                       data="ANY", ident="factor", 
                       ident.slice="factor", ident.original="factor", 
                       entropy="numeric", entropies="data.frame",
                       vds="data.frame", rds="data.frame",
                       model="ANY", gmodel="ANY", cmodel="ANY", 
                       transitions="ANY", sp.transitions="ANY", pc.transitions="ANY", 
                       profiles="ANY", sp.profiles="ANY", pc.profiles="ANY",
                       projname="character"))

					   
#' construct a slice object from an expression matrix
#'
#' construct and initialize a slice object from an expression matrix
#'
#' @param exprmatrix A numeric matrix of expression values; rows are genes, columns are cells; row names are official NCBI gene symbols, column names are cell ids/names
#' @param cellidentity A vector encoding the original identity of cells (e.g., time information or cell types); if not set, the default is a vector of "cell"
#' @param projname A string describing the analysis; if not set, the default value is "SLICE-" with timestamp information
#' @return the constructed and initialized slice object
#' @export
construct <- function(exprmatrix, cellidentity=NULL, projname=NULL) {
  
    if (is.null(cellidentity)) cellidentity <- rep("cell", dim(exprmatrix)[2])
    if (is.null(projname)) projname <- paste("SLICE-", format(Sys.time(), "%b%d_%H_%M_%S"), "-", sep="")
    
    object <- slice(exprmatrix=exprmatrix, ident.original=cellidentity, projname=projname)
    object@ident <- object@ident.original
    object@cellnames <- colnames(exprmatrix)
    object@genenames <- rownames(exprmatrix)
    
    genes <- data.frame(SYMBOL=object@genenames)
    rownames(genes) <- object@genenames
    cells <- data.frame(Cell=object@cellnames, state=object@ident)
    rownames(cells) <- object@cellnames
    object@data <- new("ExpressionSet", exprs=as.matrix(exprmatrix), phenoData=new("AnnotatedDataFrame", data=cells), featureData=new("AnnotatedDataFrame", data=genes))
    
    object@entropies <- data.frame(FAKE=rep(NA, dim(object@data)[2]))
    rownames(object@entropies) <- object@cellnames
    
    return(object)
}


#' calculate scEntropies of individual cells
#'
#' performing bootstrap or deterministic calculation of scEntropy
#'
#' @param object Slice object
#' @param cluster A data frame containing the functional clustering of genes; contains two columns "GENE" - official NCBI gene symbols and "CLUSTER" - cluster memberships.
#' @param km A symmetric matrix encoding the functional similarity of genes; the row names and column names must be official NCBI gene symbols.
#' @param calculation Choose the "bootstrap" or "deterministic" calculation of scEntropy
#' @param exp.cutoff The expression threshold of a gene to be considered expressed 
#' @param B.size In bootstrap calculation, specifies the size of bootstrap sample; in deterministic calculation, specific the number of top expressed genes from each cell
#' @param B.num The number of bootstrap samples
#' @param clustering.k The number of functional clusters
#' @param random.seed The seed for randomness
#' @return an slice object with calculated entropies in the "entropy" and "entropies" slots
#' @export
setGeneric("getEntropy", function(object, clusters=NULL, km=NULL, calculation="bootstrap", 
                                  exp.cutoff=1, B.size=1000, B.num=1, clustering.k=0, 
                                  random.seed=201602, ...) standardGeneric("getEntropy"))
#' @export
setMethod("getEntropy","slice",
          function(object, clusters=NULL, km=NULL, calculation="bootstrap", 
                       exp.cutoff=1, B.size=1000, B.num=1, clustering.k=0, 
                      random.seed=201602, ...) {
              
              es <- object@data
            
              ret <- NULL
              ret <- scEntropy(exprs(es), clusters=clusters, km=km, calculation=calculation, r=B.num,
                               exp.cutoff=exp.cutoff, n=B.size, p=0, k=clustering.k, cmethod="kmeans", 
                               prior="pr", random.seed=random.seed, context_str=object@projname, verbose=F)
              
              ename <- paste("scEntropy.", calculation, sep="")
              evalues <- ret$entropy
              object <- setEntropy(object, ename=ename, evalues=evalues)
              
              return(object)
          }
)
setGeneric("setEntropy", function(object, ename="scEntropy.bootstrap", evalues=NULL, ...) standardGeneric("setEntropy"))
setMethod("setEntropy","slice",
          function(object, ename="scEntropy.bootstrap", evalues=NULL, ...) {
            enames <- colnames(object@entropies)
            if (is.null(evalues)) {
              eid <- which(enames == ename)
              if (length(eid != 1)) {
                evalues <- object@entropies[, eid]
              } else {
                stop(paste(calculation, " entropy not found", sep))
              }
            } else {
              if (length(evalues) != length(object@ident)) {
                stop("Invalid size of input entropy vector")
              } else {
                if (length(enames)==1 & enames=="FAKE") {
                  object@entropies[, 1] <- evalues
                  colnames(object@entropies)[1] <- ename
                } else {
                  object@entropies[, ename] <- evalues
                }
                
              }
            }
            object@entropy <- evalues
            pData(object@data)$entropy <- evalues
            
            return(object)
          }
)




#' reduce expression space
#'
#' perform dimesion reduction using PCA using the expression profiles of genes satisfied the following criteria:
#'    included in the "gene.use",
#'    had log2 expression variance greater than min.var, and
#'    expressed (>1) in at least "min.cells" number of cells. 
#'
#' @param object slice object
#' @param genes.use The set of genes whose expression will be used for dimension reduction; if NULL, set to all genes
#' @param method The method for dimension reduction; currently, only support "pca"; more options will be integrated soon
#' @param num_dim The number of dimensions to be extracted
#' @param log.base The base of logarithm
#' @param do.center If TRUE, substract row-mean from each value
#' @param do.scale If TRUE, divide by row standard deviation
#' @param use.cor  If TRUE, use Pearson's correlation matrix as input to prcomp function
#' @param min.var  The minimum expression variance for a gene to be included in dimension reduction
#' @param min.cells The minimum number of expressed cells for a gene to be included in dimension reduction 
#' @return updated slice object with reduced dimensions in the "rds" slot, the first two components are also used to set the "vds" slot for visualization.
#' @export
setGeneric("getRDS", function(object, genes.use=NULL, method="pca", num_dim=2, 
                                  log.base=2, do.center=TRUE, do.scale=FALSE, 
                                  use.cor=T, min.var=0.5, min.cells=3, ...) standardGeneric("getRDS"))
#' @export
setMethod("getRDS","slice",
          function(object, genes.use=NULL, method="pca", num_dim=2, 
                   log.base=2, do.center=TRUE, do.scale=FALSE, 
                   use.cor=T, min.var=0.5, min.cells=3, ...) {
            pcs <- NULL
            if (!is.null(genes.use)) {
              genes.use <- genes.use[which(genes.use %in% object@genenames)]
            } else {
              genes.use <- object@genenames
            }
            pcs <- reduceExpressionSpace(object@data[which(object@genenames %in% genes.use), ], method=method, num_dim=num_dim, 
                                         log.base=log.base, do.center=do.center, do.scale=do.scale, 
                                         use.cor=use.cor, min.var=min.var, min.cells=min.cells, 
                                         verbose=FALSE)
            object <- setRDS(object, pcs)
            return(object)
          }
)
setGeneric("setRDS", function(object, rds, vds.x=1, vds.y=2, ...) standardGeneric("setRDS"))
setMethod("setRDS","slice",
          function(object, rds, ...) {
            
            object@vds <- rds[, c(vds.x, vds.y)]
            object@rds <- rds
            
            # for compatibility
            vds.idx <- which(substr(colnames(pData(object@data)), 1, 3) == "Dim")
            if (length(vds.idx)>0) {
              pData(object@data) <- pData(object@data)[, -vds.idx]
            } 
            pData(object@data) <- cbind(pData(object@data), rds)
            
            pData(object@data)$x <- rds[, 1]
            pData(object@data)$y <- rds[, 2]
            
            return(object)
          } 
)
                   




#' Infer entropy-directed cell lineage model 
#' 
#' SLICE provides two independent methods (clustering-based and graph-based) to perform lineage model inference 
#' 
#' The two algorithms both first divide cells into distinct cell clusters (representing cell states or cell types), find the stable state in each cluster, 
#' and then infer the lineage model by constructing a directed graph or a directed minimum spanning tree among stable states.
#' The directions of edges in the model follow entropy reduction (differentiation model) and entropy induction (dedifferentiation model)   
#' 
#' Prerequisite: the "vds", "rds", and "entropy" slots in the input slice object must be set
#'
#' @param object Slice object 
#' @param lm.method Select "clustering" based or "graph" based method to infer lineage model
#' @param model.type The type of models that will be infered: "tree" - directed minimum spanning tree based, "graph" - directed graph based 
#' @param reverse If FALSE, infer differentiation lineage model (edge directions follow entropy reduction), if TRUE, infer dedifferentiation lineage model (edge directions follow entropy induction)
#' @param ss.method The method for defining core cell set for stable state detection: 
#'                   all - all the cells in a cluster constitute the core cell set; 
#'                   top - cells with scEntropy lower than the ss.threshold quantile of all the values in a cluster constitute the core cell set; 
#'                   pcst - cells with scEntropy lower than the ss.threshold quantile of all the values in a cluster constitute the prize nodes, linear prize-collecting steiner tree algorithm is used to approximate an optimal subnetwork, the cells in the subnetwork constitute the core cell set. Stable states are defined as the centroids of the core cell sets.   
#' @param ss.threshold The threshold used when ss.method is "top" or "pcst". Default: 0.25.
#' @param wiring.threshold The threshold for the local wiring procedure in the cell-cell network construction. 
#'                         By default, it is set to the longest edge in the minimum spanning treee.
#'                         Only take effect when lm.method is "grpah"
#' @param community.method The method for network community detection. 
#'                         Most of the community detection methods implemented in the igraph package are supported, 
#'                         including "fast_greedy", "edge_betweenness", "label_prop", "leading_eigen","louvain","spinglass", "walktrap". 
#'                         If this parameter is set to "auto", the algorithm will perform all the community detection methods and select the one that generates the communities with best modularity. 
#'                         Only take effect when lm.method is "graph"
#' @param cluster.method Use "kmeans" or "pam" to divide cells into clusters. Only take effect when lm.method is "clustering" 
#' @param k The number of cell clusters. If NULL, Gap statistic will be used to determine an optimal k.
#' @param k.max The "k.max" parameter of cluster::clusGap(); used when k is NULL.
#' @param B The "B" parameter of cluster::clusGap(); used when k is NULL
#' @param k.opt.method The "method" parameter of cluster::maxSE(); used when k is NULL
#' @return updated slice object with inferred lineage model in the "model" slot
#' @export
setGeneric("getLineageModel", function(object, lm.method="clustering", model.type="tree", reverse=F, ss.method="all", ss.threshold=0.25, # common parameter
                                       wiring.threshold=function(mst) max(mst), community.method="louvain",                # parameters for graph-based method                        
                                       cluster.method="kmeans", k=NULL, k.max=10, B=100, k.opt.method="firstmax",          # parameters for clustering-based method
                                       ...) standardGeneric("getLineageModel"))
#' @export
setMethod("getLineageModel","slice",
          function(object, lm.method="clustering", model.type="tree", reverse=F, ss.method="all", ss.threshold=0.25, # common parameter
                            wiring.threshold=function(mst) max(mst), community.method="louvain",                     # parameters for graph-based method                        
                            cluster.method="kmeans", k=NULL, k.max=10, B=100, k.opt.method="firstmax",               # parameters for clustering-based method
                            ...) {
            
              lmmethods <- c("clustering","graph")
              mid <- pmatch(lm.method, lmmethods)
              if (is.na(mid)) {
                stop("Invalid lm.method. Please select \"cluster\" or \"graph\"")
              }
              lm.method=lmmethods[mid]
              
              lm = NULL
              if (lm.method=="graph") {
                lm <- getLM.graph(object@data, model.type=model.type, wiring.threshold=wiring.threshold, 
                                  community.method=community.method,                                            
                                  ss.method=ss.method, ss.threshold=ss.threshold,                                  
                                  reverse = reverse, do.plot=T, context_str=object@projname)
              } else if (lm.method=="clustering") {
                lm <- getLM.clustering(object@data, model.type=model.type, cluster.method=cluster.method, 
                                       k=k, k.max=k.max, B=B, k.opt.method=k.opt.method, 
                                       ss.method=ss.method, ss.threshold=ss.threshold,                                   
                                       reverse = reverse, do.plot=T, context_str=object@projname) 
              }
              
              
              object <- setLineageModel(object, lm.method=lm.method, lm=lm)
              return(object)
          }
)
setGeneric("setLineageModel", function(object, lm.method="clustering", lm=NULL, ...) standardGeneric("setLineageModel"))
setMethod("setLineageModel","slice", 
          function(object, lm.method="clustering", lm=NULL, ...) {
            if (is.null(lm)) {
              if (lm.method=="clustering") {
                if(is.null(object@cmodel)) {
                  stop("Invalid model")
                } else {
                  lm <- object@cmodel
                }
              } else if (lm.method=="graph") {
                if(is.null(object@gmodel)) {
                  stop("Invalid model")
                } else {
                  lm <- object@gmodel
                }
              } else {
                stop("Invalid lm.method")
              }
            }
            
            if (lm.method=="clustering") {
              object@cmodel <- lm
            } else if (lm.method=="graph") {
              object@gmodel <- lm
            } else {
              stop("Invalid lm.method")
            }
            object@model <- lm
            
            return(object)
          }
)


#' Reconstructing cell trajectories following a specific lineage (start -> end) in the inferred lineage model 
#' 
#' This function implements both shortest-path based and principal-curve based cell transitional path reconstruction.
#' 
#' Prerequisite: the "model" slot of the input slice object must be set
#'
#' @param object Slice object
#' @param method Select either "sp" (shortest-path based) method or "pc" (principal-curve based) method to reconstruct cell trajectories following a lineage in the infered lineage model 
#' @param start The numeric ID of the stable state from which the cell transition starts 
#' @param end The numeric ID of the stable state at which the cell transition ends 
#' @param network Used in the shortest path based method; "mst" means the shortest path will be explored on the minimum spanning tree of cells, while "csn" means on the cell similarity network (i.e., mst + local wiring) 
#' @param NN.threshold Used in the shortest path based method, specifies the distance threshold for determining neighbors of each cell in the path. Default is 0.7 quantile of all edge weights in the network 
#' @param do.trim Used in the principal curve based method. When set to TRUE, only the curve between the start and end stable states will be returned; otherwise, the whole curve will be returned.
#' @return updated slice object with reconstructed cell trajectories in the "transitions" slot.
#' @export
setGeneric("getTrajectories", function(object, method="sp", start, end, 
                                      network="mst", NN.threshold=0.7,  # sp parameters
                                      do.trim=F,                                      # pc parameters
                                      ...) standardGeneric("getTrajectories"))
#' @export
setMethod("getTrajectories","slice", 
          function(object, method="sp", start, end, 
                           network="mst", NN.threshold=0.7, # sp parameters
                           do.trim=F,                                      # pc parameters
                           ...) {
  
              methods <- c("sp","pc")
              mid <- pmatch(method, methods)
              if (is.na(mid)) {
                stop("Invalid method. Please select \"sp\" or \"pc\"")
              }
              method=methods[mid]
              
              model=NULL
              if (is.null(model)) {
                model <- object@model
              }
               
              trajectories <- NULL
              NN.type="all"
              do.plot=T
              if (method=="sp") {
                trajectories <- get.SP.Transitions(object@data, model, start, end, network=network, NN.threshold=NN.threshold, NN.type=NN.type, 
                                   do.plot=T, context_str=object@projname)
              } else if (method=="pc") {
                trajectories <- get.PC.Transitions(object@data, model, start, end, do.trim=do.trim, 
                                   do.plot=T, context_str=object@projname)
              }
              object <- setTrajectories(object, method, trajectories=trajectories)
              return(object)
         }
)
setGeneric("setTrajectories", function(object, method="sp", trajectories=NULL, ...) standardGeneric("setTrajectories"))
setMethod("setTrajectories","slice", 
          function(object, method="sp", trajectories=NULL, ...) {
            if (method=="sp") {
              object@sp.transitions <- trajectories
              object@transitions <- trajectories
            } else if (method=="pc") {
              object@pc.transitions <- trajectories
              object@transitions <- trajectories
            }
            return(object)
          }
)



#' Exract lineage dependent gene expression profiles
#' 
#' Call after getTrajectories function to extract the gene expression profiles in each of the reconstructed cell trajectories
#'
#' @param object Slice object 
#' @param trajectory.type The extract algorithm will be different depending on the type of cell trajectories: sp (shortest-path based) or pc (principal-curve based)
#' @param genes.use If set, only the expression profiles of the specified set of genes will be extracted; otherwise, the expression profiles of all genes will be extracted 
#' @return updated slice object with extracted expression profiles in the "profiles" slot.
#' @export
setGeneric("getProfiles", function(object, trajectory.type="sp", genes.use=NULL, ...) standardGeneric("getProfiles"))
#' @export
setMethod("getProfiles","slice",
          function(object, trajectory.type="sp", genes.use=NULL, ...) {
  
            trajectories=NULL
            
            types <- c("sp","pc")
            mid <- pmatch(trajectory.type, types)
            if (is.na(mid)) {
              stop("Invalid transition.type. Please select \"sp\" or \"pc\"")
            }
            trajectory.type=types[mid]
            
            if (is.null(genes.use)) {
              genes.use <- object@genenames
            } else {
              genes.use <- genes.use[which(genes.use %in% object@genenames)]
            }
            
            if (is.null(trajectories)) {
              if (trajectory.type=="sp") {
                trajectories <- object@sp.transitions
              } else if (trajectory.type=="pc") {
                trajectories <- object@pc.transitions
              }
            }
            
            profiles <- NULL
            if (trajectory.type=="sp") {
              profiles <- get.SPT.Profiles(object@data[genes.use, ], transitions=trajectories$transitions, 
                                           context_str=object@projname, do.plot=T, plot.xlabel=NULL, plot.ylabel=NULL, plot.w=2.8, plot.h=1)
            } else if (trajectory.type=="pc") {
              profiles<- get.PCT.Profiles(object@data[genes.use,], transitions=trajectories, 
                               context_str=object@projname, do.plot=T, plot.xlabel=NULL, plot.ylabel=NULL, plot.w=2.8, plot.h=1)  
            }
            object <- setProfiles(object, trajectory.type=trajectory.type, profiles=profiles)
            return(object)
          }
)
setGeneric("setProfiles", function(object, trajectory.type="sp", profiles, ...) standardGeneric("setProfiles"))
setMethod("setProfiles","slice", 
          function(object, trajectory.type="sp", profiles, ...) {
            if (trajectory.type=="sp") {
              object@sp.profiles <- profiles
              object@profiles <- object@sp.profiles
            } else if (trajectory.type=="pc") {
              object@pc.profiles <- profiles
              object@profiles <- object@pc.profiles
            }
            return(object)
          }
)


#' Detecting lineage dependent differentially expressed genes
#' 
#' given the lineage dependent expression profiles extracted by the getProfiles function, this function first smoothes the expression profile of each gene using cubic regression splines, 
#' fits to two models (M1: changes in the lineage; M0: no change),
#' and then uses likelihood ratio test to compare the goodness of fit of the two models and assigns a p-value to each gene
#' 
#' Genes meet the following criteria will be selected for likelihood ratio test: expression >=thresh.exp in at least thresh.cell cells and variance>=thresh.var
#' Genes with FDR<diff.thresh will be considered significant  
#'
#' @param profiles The lineage dependent gene expression profiles 
#' @param diff.thresh The FDR threshold for significance
#' @param thresh.exp The expression threshold
#' @param thresh.cell The cell number threshold
#' @param thresh.var The variance threshold
#' @return lineage dependent differentially expressed genes and their associated information, e.g., variance, pvalue, fdr, and M1 predicted expression
#' @export  
getDiffGenes <- function(profiles, diff.thresh=0.1, log.base=2, thresh.exp=1, thresh.cell=2, thresh.var=0.5) {
    
	cat("\nDetecting lineage dependent differentially expressed genes\n")
	
    n <- length(profiles)
    
    ret <- list()
     
    for (i in 1:n) {
        
        ret[[i]] <- list()
        cat("\nDetecting differentially expressed genes in transition: ", names(profiles)[i], "\n", sep="")
        
        # obtain the expression data in the transition ss1 -> ss2
        exprs <- profiles[[i]]$exprs
        ptime <- profiles[[i]]$ptime
        cells <- profiles[[i]]$cells
        genes <- profiles[[i]]$genes
        
        if ((log.base %% 1==0) & (log.base>1)) {
            exprs <- log(exprs+1, log.base)
        }
        
        diff <- data.frame(Gene=rownames(exprs), check.names=FALSE)
        rownames(diff) <- diff$Gene
        
        # expressed (>=1) in at least two cells in the lineage
        #num_expressed <- as.numeric(apply(exprs, 1, function(x) length(which(x>=thresh.exp))) )
        
        
        diff$var <- NA
        diff$pvalue <- NA
        diff$fdr <- NA
        diff$significant <- 0
        
        k=NULL
        exprs.fit <- t(apply(exprs, 1, function(x) fit_helper(x, ptime=ptime, k=k)))
        
        fit.pvalues <- exprs.fit[, 1]
        exprs.fit <- exprs.fit[, 2:dim(exprs.fit)[2]]
        if (dim(exprs.fit)[2] == dim(exprs)[2]) {
            colnames(exprs.fit) <- colnames(exprs)
        } else {
            colnames(exprs.fit) <- paste("T",unique(ptime),sep="")
        }
        
        num_expressed <- rowSums(exprs.fit >= thresh.exp)
        diff$num_expressed <- num_expressed
        
        # variance of the fitted expression
        fit.vars <- apply(exprs.fit, 1, var)
        
        diff$var <- fit.vars
        diff$pvalue <- fit.pvalues
        
        exp.idx <- which(num_expressed>=thresh.cell)  
        var.idx <- which(fit.vars>=thresh.var)
        
        diff.idx <- intersect(exp.idx, var.idx)
        
        exprs.fit.diff <- NULL
        
        if (length(diff.idx) >0 ) {
            
            if (length(diff)==1) {
                diff$fdr[diff.idx] <- diff$pvalue[diff.idx]
            }
            
            diffgenes <- diff[diff.idx, ]
            diffgenes$fdr <- p.adjust(diffgenes$pvalue, method="fdr")
            
            diff[as.character(diffgenes$Gene), "fdr"] <- diffgenes$fdr
            
            diff.idx <- which(diffgenes$fdr < diff.thresh)
            
            if (length(diff.idx) > 0) {
                
                cat(paste(length(diff.idx), " diff genes detected (num_expressed>=", thresh.cell, ", var(fit)>=", thresh.var, ", fdr<", diff.thresh, ")\n", sep=""))
                
                diff.genes <- as.character(diffgenes$Gene[diff.idx])
                diff[diff.genes, "significant"] <- 1
                exprs.fit.diff <- exprs.fit[diff.genes, ]
            }
        }
        
        ret[[i]] <- list(exprs.fit=exprs.fit, exprs.fit.diff=exprs.fit.diff, diff=diff)
    }
    
    names(ret) <- names(profiles)
    
    return(ret)
    
}




#' Discover lineage dependent temporal patterns
#' 
#' Given the M1 predicted expression profiles of differentially expressed genes in a transition, 
#' the function uses PAM clustering to divide differentially expressed genes into clusters (patterns)  
#'
#' @param expr The M1 predicted expression profiles of differentially expressed genes in a transition
#' @param do.zscore If TRUE, perform gene based zscore transformation prior to clustering
#' @param c.method The algorithm to find pattern clusters, currently support "pam" and "kmeans" 
#' @param k The number of patterns to be discovered; if NULL, gap statistic (cluster::clusGap) will be used to determine an optimal k; parameters of cluster::clusGap are set as follows: k.max=10, B=50
#' @param d.method The distance method
#' @param do.plot Whether to plot discovered patterns
#' @param plot.filename The file name of the plot; if NULL, plot to the screen 
#' @return An object of class "pam" or "kmeans" representing the clustering results.
#' @export
getPatterns <- function(expr, do.zscore=T, c.method="pam", k=NULL, d.method=function(x) as.dist((1 - cor(t(x)))/2), do.plot=T, plot.filename=NULL) {
    
	cat("\nDiscovering lineage dependent temporal patterns\n")
	
	mat <- expr
    
    if (do.zscore==T) {
        zscore_helper <- function(x) {
            x <- as.numeric(x)
            if (sd(x)==0) {
                x <- rep(0, length(X))
            } else {
                x <- (x-mean(x, na.rm = T))/sd(x, na.rm = T)
            }
            return(x)
        }
        x <- t(apply(mat, 1, function(y) zscore_helper(y)))
        rownames(x) <- rownames(mat)
        colnames(x) <- colnames(mat)
        mat <- x
    }
    
    if (is.null(k)) {
        k.max <- min(dim(expr)[1]-1, 10)
        B=50
        if (c.method=="pam") {
            gskmn <- cluster::clusGap(as.matrix(mat), FUN = pam,  K.max = k.max, B = B)
        } else if (c.method=="kmeans") {
            gskmn <- cluster::clusGap(as.matrix(mat), FUN = kmeans,  K.max = k.max, B = B)
        }
        plot(gskmn, main = paste("Number of temporal patterns determined by\n Gap statistic (FUN=",c.method,", k.max=", k.max, ", B=",B, ")", sep=""))
        
        # default: firstmax rule
        k <- maxSE(f=as.numeric(gskmn$Tab[,"gap"]), SE.f=as.numeric(gskmn$Tab[,"SE.sim"]), method="firstmax")
        abline(v=k, col="blue", lty=2, lwd=3)
        text(k, min(gskmn$Tab[,"gap"]), paste(" k=",k, "\n(firstmax)", sep=""), col="blue", adj = c(-.1, -.1), cex=1.1)
    }
    
    if (c.method == "pam") {
        db <- pam(as.matrix(mat), k)
        cluster <- db$clustering
    } else if (c.method=="kmeans") {
        db <- kmeans(as.matrix(mat), k)
        cluster <- db$cluster
    }
	
    if (do.plot==T) {
      
        plot.xlabel="Transitional Path"
        plot.ylabel="Expression (+1, log2)"
        plot.w=3
        plot.h=4
        context_str=""
      
        dd <- as.data.frame(expr)
        dd$gene <- rownames(expr)
        dd$cluster <- factor(cluster)
        
        times <- seq(0, 1, length.out=dim(expr)[2])
        ptime <- data.frame(cell=colnames(expr), ptime=times)
        
        dd.melt <- melt(dd, id.vars = c("gene", "cluster"))
        dd.melt <- merge(dd.melt, ptime, by.x="variable", by.y="cell")
        dd.melt$cluster <- paste("Pattern ", dd.melt$cluster, sep="")
        
        g <- ggplot(dd.melt) + facet_wrap("cluster") 
        
        g <- g + geom_line(aes(x = ptime, y = value, group=gene), color="grey80", alpha=0.5, size=I(0.1)) 
        
        g <- g + stat_summary(aes(x = ptime, y = value, group = 1), 
                              fun.data = mean_sdl, color = "black", fill = "black", 
                              alpha = 0.1, size = I(1), geom = "smooth")
        
        g <- g + theme(axis.text.x = element_text(angle = 0, hjust = 0)) + xlab(plot.xlabel) + ylab(plot.ylabel)
        
        g <- g + SLICE_theme_opts()
        
        g <- g + theme(strip.background = element_rect(colour = 'white', fill = 'white')) + 
            theme(panel.border = element_blank(), axis.line.x = element_line(size=0.2), axis.line.y = element_line(size=0.2)) +
            theme(axis.ticks = element_line(size = 0.2)) +
            theme(legend.position="none") +
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
            theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
        
        
        if (is.null(plot.filename)) {
            print(g)
        } else {
            ggsave(file=paste(plot.filename, ".pdf", sep=""), width=plot.w, height=plot.h)
        }
    }
    
    return(db)
}



#' Visualize and compare entropies
#' 
#' Visualize entropies using boxplot and compare using regression  
#'
#' @param object Slice object containing the entropies in object@entropies
#' @export
setGeneric("plotEntropies", function(object, ...) standardGeneric("plotEntropies"))
#' @export
setMethod("plotEntropies","slice",
          function(object, ...) {
              en <- dim(object@entropies)[2]
              ee <- data.frame(object@entropies, ident=object@ident)
              
              enames <- colnames(object@entropies)
              
              theme_border <- function() {
                  theme(panel.border=element_rect(color="grey", fill=NA, size=0.5))
              }
              
              pdf(file=paste(object@projname, "entropies.pdf", sep=""), width=4*en, height=3*en)
            
              gs <- list()
              for (i in 1:en) {
                  for (j in 1:en) {
                      if (i == j) {
                          gs[[(i-1)*en + j]] <- ggplot(data=ee) + geom_boxplot(aes_string(x="ident", y=enames[i], col="ident")) + SLICE_theme_opts() + theme_border()  
                      } else {
                          gs[[(i-1)*en + j]] <- ggplot(data=ee, aes_string(x=enames[i], y=enames[j], col="ident")) + geom_point() + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, size=0.5) + SLICE_theme_opts() + theme_border()
                      }
                  }
              }
             
              grid.arrange(grobs=gs, ncol=en)
              
              dev.off()
          }
)





#######################################################
#              Detailed Implementations               #
#######################################################


scEntropy <- function(exp.m, clusters=NULL, km=NULL, calculation="bootstrap", r=1,
                      exp.cutoff=1, n=1000, p=0, k=0, cmethod="kmeans", prior="pr",
                      random.seed=201602, context_str="", verbose=F) {
  
  if (is.null(clusters) & is.null(km)) {
    stop("Please provide either a clustering of genes or a similarity matrix of genes")
  }
  
  calculation.methods <- c("bootstrap","deterministic")
  sid <- pmatch(calculation, calculation.methods)
  if (is.na(sid)) {
    stop("Please choose bootstrap or deterministic calculation of scEntropy.")
  }
  calculation <- calculation.methods[sid]
  
  if (!is.null(random.seed)) {
    set.seed(random.seed)
  }
  
  cmethods <- c("kmeans","pam")
  sid <- pmatch(cmethod, cmethods)
  if (is.na(sid)) {
    stop("Please choose kmeans or pam for gene clustering")
  }
  cmethod <- cmethods[sid]
  
  
  cat("\nPerforming ", calculation, " calculation of scEntropy\n", sep="")
  
  ncells <- dim(exp.m)[2]
  
  if ((r > 0) & (r %% 1==0)) { 
    
    ee <- data.frame(r1=rep(NA, ncells)) # storing entropy in each iteration
    rownames(ee) <- colnames(exp.m)
    
    for (rr in 1:r) { # for each iteration
      
      cat("\rIteration:", rr, sep="")
      
      genelists <- NULL
      
      # select genes per cell
      if (calculation=="deterministic") {
        genelists <- expressed.top(exp.m, exp.cutoff=exp.cutoff, k=n, p=p)
      } else if (calculation=="bootstrap") {
        genelists <- expressed.random(exp.m, exp.cutoff=exp.cutoff, k=n, p=p)
      } else {
        stop("Invalid scEntropy calculation method.")
      }
      
      file.prefix <- context_str
      
      genes <- unique(as.character(unlist(genelists))) # all genes selected
      
      gclus <- NULL
      if (is.null(clusters)) { # if gene clustering is not provided
        
        genes <- genes[which(genes %in% rownames(km))] # keep genes with annotation
        sim <- km[as.character(genes), as.character(genes)]
        sim.dist <- 1-sim
        sim.dist <- as.dist(sim.dist)
        
        if (k<2) {
          kk = max(2, floor(sqrt(dim(sim)[1]/2)))
        } else {
          kk <- k
        }
        
        if (cmethod=="pam") {
          gclus <- pam(sim.dist, kk, cluster.only=TRUE)
        } else if (cmethod=="kmeans") {
          gclus <- kmeans(sim.dist, centers=kk, iter.max=10, nstart=1)$cluster
        } else {
          stop("Please use \"pam\" or \"kmeans\" in cmethod.\n")
        }
        gclus <- data.frame(GENE=names(gclus), CLUSTER=as.numeric(gclus))
        
      } else {
        gclus <- clusters
        
        # only keep genes in the clustering scheme
        genes <- genes[which(genes %in% as.character(gclus$GENE))]
        
        gclus <- gclus[which(as.character(gclus$GENE) %in% genes), ]
        gclus$CLUSTER <- factor(gclus$CLUSTER)
      }
      
      if (verbose) {
        write.table(gclus, file=paste(file.prefix, rr, "-1-gclus.0.txt", sep=""), sep="\t", col.names=T, row.names=F)
      }
      
      gclus <- gclus[order(gclus$CLUSTER), ]
      gclus.tbl <- table(gclus$CLUSTER)
      gclus.summary <- data.frame(CLUSTER=as.numeric(names(gclus.tbl)), COUNT=as.numeric(gclus.tbl)) # make sure CLUSTER is numeric vector, 12/08/2016
      gclus.summary$PERCENTAGE <- gclus.summary$COUNT/sum(gclus.summary$COUNT)
      
      gclus.summary <- gclus.summary[order(gclus.summary$CLUSTER), ]
      
      if (verbose) {
        write.table(gclus, file=paste(file.prefix, rr, "-1-gclus.1.txt", sep=""), sep="\t", col.names=T, row.names=F)
        write.table(gclus.summary, file=paste(file.prefix, rr, "-1-gclus.2.txt", sep=""), sep="\t", col.names=T, row.names=F)
      }
      
      if (verbose) {
        print(table(gclus$CLUSTER))
      }
      
      g.exp.m <- exp.m[as.character(gclus$GENE), ] # keep expression profiles with cluster information
      
      g.exp.m.b <- g.exp.m                         # binarize
      g.exp.m.b[g.exp.m.b<=exp.cutoff] <- 0
      g.exp.m.b[g.exp.m.b>0] <- 1
      
      
      if (verbose) {
        write.table(data.frame(SYMBOL=rownames(g.exp.m), CLUSTER=gclus$CLUSTER, g.exp.m), file=paste(file.prefix, rr, "-2-exp.m.txt", sep=""), sep="\t", row.names=F, col.names=T)
      }
      
      y <- aggregate(g.exp.m.b, by=list(gclus$CLUSTER), sum) 
      y <- y[order(y[,1]), ]                                 
      
      if (verbose) {
        write.table(y, file=paste(file.prefix, rr, "-3-exp.cnts.txt", sep=""), sep="\t", row.names=F, col.names=T)
      }
      
      y.clusters <- y[, 1]
      y <- y[, 2:dim(y)[2]] # each value represents the number of expressed genes in each cluster in each cell
      
      
      entropy <- rep(NA, ncells)


      if (prior=="pr") { 
        a <- as.numeric(gclus.summary$PERCENTAGE)
      } else if (prior=="no") { 
		a <- rep(0, dim(gclus.summary)[1])
	  }
      
      for (i in 1:ncells) {
        entropy[i] <- entropy::entropy.Dirichlet(y[,i], a=a)
      }
      
      ee[, rr] <- as.numeric(entropy)
      rr <- rr + 1
    }
    
  } else {
    stop("invalid value for r")
  }
  cat("\n")
  return(list(entropy=apply(ee, 1, mean), entropies=ee))
}


expressed.all <- function(exp.m, exp.cutoff=1) {
  genelists <- list()
  n <- dim(exp.m)[2]
  cells <- colnames(exp.m)
  genes <- rownames(exp.m)
  
  for (i in 1:n) {
    i.idx <- which(exp.m[,i] > exp.cutoff)
    if (length(i.idx) > 0) {
      genelists[[i]] <- genes[i.idx]
    } else {
      genelists[[i]] <- NA
    }
    names(genelists)[i] <- cells[i]
  }
  
  return(genelists)
}


expressed.top <- function(exp.m, exp.cutoff=1, k=0, p=0.2) {
  genelists <- list()
  n <- dim(exp.m)[2]
  cells <- colnames(exp.m)
  genes <- rownames(exp.m)
  
  for (i in 1:n) {
    i.idx <- which(exp.m[,i] > exp.cutoff)
    if (length(i.idx) > 0) {
      if (!is.null(k) & k>1) {
        if (length(i.idx) <= k) {
          genelists[[i]] <- genes[i.idx]
        } else {
          i.exp <- exp.m[i.idx, i]
          genelists[[i]] <- names(i.exp[order(i.exp, decreasing=T)[1:k]])
        }
      } else if (p>0 & p<=1) {
        kk <- max(floor(length(i.idx)*p), 2)
        i.exp <- exp.m[i.idx, i]
        genelists[[i]] <- names(i.exp[order(i.exp, decreasing=T)[1:kk]])
      }
      # genelists[[i]] <- genes[i.idx]
    } else {
      genelists[[i]] <- NA
    }
    names(genelists)[i] <- cells[i]
  }
  
  return(genelists)
}


expressed.random <- function(exp.m, exp.cutoff=1, k=0, p=0.5) {
  genelists <- list()
  n <- dim(exp.m)[2]
  cells <- colnames(exp.m)
  genes <- rownames(exp.m)
  
  genelists <- expressed.all(exp.m, exp.cutoff=exp.cutoff)
  
  allgenes <- sort(unique(as.character(unlist(genelists))))
  
  rndgenes <- NULL
  if ((k>2) & (k %% 1==0)) {
    if (length(allgenes) < k) {
      rndgenes <- allgenes
    } else {
      rndgenes <- sample(allgenes, k, replace=FALSE)
    }
  } else {
    kk <- max(10, floor(length(allgenes) * p))
    rndgenes <- sample(allgenes, kk, replace=FALSE)
  }
  
  genelists <- list()
  for (i in 1:n) {
    genelists[[i]] <- rndgenes
    names(genelists)[i] <- cells[i]
  }
  
  return(genelists)
}


reduceExpressionSpace <- function(es, method="pca", num_dim=2, 
                                  log.base=2, do.center=TRUE, do.scale=FALSE, 
                                  use.cor=T, min.var=0.5, min.cells=3, 
                                  verbose=FALSE) {
  
  cat("\nPerforming dimension reduction\n")
  
  if (!(method == "pca")) {
    stop("Please set method to pca. Other methods will be supported soon!")
  }
  
  if (!(log.base %% 1 == 0) | log.base<1) {
    stop("Invalid log.base. please select an integer greater than zero")
  }
  
  N <- dim(fData(es))[1] # total number of features
  
  exp.m <- exprs(es)
  
  # abundancy critertia: 
  temp <- exp.m
  temp <- log(temp+1, base=2)
  
  if (is.null(min.cells)) {
    min.cells <- dim(temp)[2]*0.3
  }
  idx.expressed <- which(rowSums(temp>1)>=min.cells)
  idx.variable <- which(apply(temp, 1, var) > min.var)
  
  idx <- intersect(idx.expressed, idx.variable)
  
  rm(temp)
  if (length(idx) >= 3) {
    if (verbose) {
      cat(length(idx), " profiles used for dimension reduction\n", sep="")
    }
  } else {
    stop("Insufficent genes for dimension reduction.")
  }
  
  exp.m <- exp.m[idx, ]
  
  
  if (log.base > 1) {
    exp.m <- log(exp.m+1, log.base) # to avoid zero, increase by 1, and then log
  }
  
  if (do.center==TRUE | do.scale==TRUE) {
    exp.m <- t(scale(t(exp.m), center=do.center, scale=do.scale))
  }
  
  pcs <- NULL
  
  if (use.cor==TRUE) {
    mat <- cor(exp.m)
  } else {
    mat <- t(exp.m)
  }
  rm(exp.m)
  
  pca <- prcomp(mat, retx = TRUE, center = F, scale. = F)
  pcs <- pca$x[, 1:num_dim]
  
  
  pcs <- as.data.frame(pcs, check.names=FALSE)
  colnames(pcs) <- paste("Dim", 1:num_dim, sep="")
  
  return(pcs)
}


getLM.graph <- function(es, model.type="tree", 
                        wiring.threshold=function(mst) max(mst), 
                        community.method="louvain",                                            
                        ss.method="all", ss.threshold=0.25,                                  
                        reverse = FALSE, 
                        do.plot=T, context_str="") {
  
  ssmethods <- c("all", "top", "pcst")
  ssid <- pmatch(ss.method, ssmethods)
  if (is.na(ssid)) {
    stop("Invalid stable state definition method.")
  }
  ss.method <- ssmethods[ssid]
  
  if (is.null(pData(es)$state)) {
    stop("Missing column pData(es)$state, which encodes the type/state/cluster information of cells.")
  }
  
  if (is.null(pData(es)$entropy)) {
    stop("Missing column pData(es)entropy, which encodes the entropy of cells.")
  }
  
  if (do.plot==TRUE) {
    if (is.null(pData(es)$x) | is.null(pData(es)$y)) {
      stop("do.plot is set to TRUE, columns pData(es)$x and pData(es)$y are required for visualization.")
    }
  }
  
  dim.cols.idx <- which(substr(colnames(pData(es)), 1, 3)=="Dim")
  if (length(dim.cols.idx)<2) {
    stop("At least two columns with prefix \"Dim\" in pData(es) are required.")
  }
  dim.cols <- colnames(pData(es))[dim.cols.idx]
  
  cat("\nInferring lineage model using graph-based method\n")
  
  
  if (do.plot) {
    pdf.w <- 8
    pdf.h <- 6
    
    dd <- pData(es)[, c("state", "entropy","x", "y")]
    dd$Cluster <- dd$state
    dd$Entropy <- dd$entropy/max(dd$entropy)
    dd$z <- dd$entropy
    
    g <- ggplot(data=dd, aes(x=x, y=y)) + geom_point(aes(color=Cluster, size=Entropy)) + SLICE_theme_opts() + ggtitle("Cells") + labs(x="PC1", y="PC2") #+ geom_text(label=rownames(dd), size=2, hjust=-0.2, vjust=-0.2) ;
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    ggsave(filename=paste(context_str, "lineage-g.step1-cells.pdf", sep=""), width=pdf.w, height=pdf.h) 
    print(g)
  }
  
  
  ################################################## 
  # 1. generating the cell similarity network
  ##################################################
  
  reducedDims <- t(pData(es)[, dim.cols.idx])
  
  nets <- getCellSimilarityNetwork(reducedDims, dist.method="euclidean", wiring.threshold=wiring.threshold)
  
  mst <- nets$mst
  gp <- nets$csn
  
  if (do.plot) { # minimum spanning tree
    
    # visualize the similarity network
    edge_list <- as.data.frame(get.edgelist(mst))  
    colnames(edge_list) <- c("source", "target")
    
    dd <- pData(es)[, c("state", "entropy","x", "y")]
    dd$color <- dd$state
    dd$z <- dd$entropy
    if(is.null(dd$CELL)) {
      dd.df <- cbind(CELL=rownames(dd), dd)
    } else {
      dd.df <- dd
    }
    
    edge.df <- merge(dd.df, edge_list, by.x="CELL", by.y="source", all=T)
    edge.df <- edge.df[,c("CELL","state", "x","y","z","target")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "target")
    edge.df <- merge(edge.df, dd.df, by.x="target", by.y="CELL", all=T)
    edge.df <- edge.df[,c("source", "source.state", "source.x", "source.y", "source.z", "target", "state", "x", "y", "z")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "target","target.state","target.x", "target.y", "target.z")
    
    edge.df$Cluster=edge.df$source.state
    edge.df$Entropy=edge.df$source.z
    g <- ggplot(data=edge.df, aes(x=source.x, y=source.y)) + ggtitle("Minimum Spanning Tree") + labs(x="PC1", y="PC2") + SLICE_theme_opts();
    g <- g + geom_segment(aes_string(xend="target.x", yend="target.y"), size=I(0.5), linetype="solid", col=I("grey60"), na.rm=TRUE)
     g <- g + geom_point(aes(colour=Cluster, size=Entropy),  na.rm=TRUE) 
    #g <- g + geom_point(aes(colour=Entropy), size=5,  na.rm=TRUE) 
    # g <- g + geom_text(aes_string(label="source"), size=2, hjust=-0.2, vjust=-0.2, na.rm=TRUE)
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    ggsave(filename=paste(context_str, "lineage-g.step2-mst_of_cells.pdf", sep=""), width=pdf.w, height=pdf.h) #, dpi=dpi, compression="lzw")
    
    # g
    print(g)
    
  }
  
  if (do.plot) { # cell similarity network
    
    # visualize the similarity network
    edge_list <- as.data.frame(get.edgelist(gp))  
    colnames(edge_list) <- c("source", "target")
    
    dd <- pData(es)[, c("state", "entropy","x", "y")]
    dd$color <- dd$state
    dd$z <- dd$entropy
    if(is.null(dd$CELL)) {
      dd.df <- cbind(CELL=rownames(dd), dd)
    } else {
      dd.df <- dd
    }
    
    edge.df <- merge(dd.df, edge_list, by.x="CELL", by.y="source", all=T)
    edge.df <- edge.df[,c("CELL","state", "x","y","z","target")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "target")
    edge.df <- merge(edge.df, dd.df, by.x="target", by.y="CELL", all=T)
    edge.df <- edge.df[,c("source", "source.state", "source.x", "source.y", "source.z", "target", "state", "x", "y", "z")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "target","target.state","target.x", "target.y", "target.z")
    
    edge.df$Cluster=edge.df$source.state
    edge.df$Entropy=edge.df$source.z
    g <- ggplot(data=edge.df, aes(x=source.x, y=source.y)) + ggtitle("Cell-cell network") + labs(x="PC1", y="PC2") + SLICE_theme_opts();
    g <- g + geom_segment(aes_string(xend="target.x", yend="target.y"), size=I(0.5), linetype="solid", col=I("grey60"), na.rm=TRUE)
    g <- g + geom_point(aes(colour=Cluster, size=Entropy),  na.rm=TRUE) 
    #g <- g + geom_point(aes(colour=Entropy), size=5, na.rm=TRUE) 
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    ggsave(filename=paste(context_str, "lineage-g.step3-cell_network.pdf", sep=""), width=pdf.w, height=pdf.h) #, dpi=dpi, compression="lzw")
    
    # g
    print(g) 
  }
  
  ################################################## 
  # 2. identify cell clusters
  ##################################################
  
  cmm <- getCellCommunities(gp, method=community.method)
  wc.df <- data.frame(CELL=names(cmm), WC=paste("C", as.numeric(cmm), sep=""))  
  
  pData(es)$slice.state <- paste("C", as.numeric(cmm), sep="")
  
  if (do.plot) {
    
    # visualize the similarity network
    edge_list <- as.data.frame(get.edgelist(gp))  
    colnames(edge_list) <- c("source", "target")
    
    dd <- pData(es)[, c("state", "entropy","x", "y")]
    dd$color <- dd$state
    dd$z <- dd$entropy
    if(is.null(dd$CELL)) {
      dd.df <- cbind(CELL=rownames(dd), dd)
    } else {
      dd.df <- dd
    }
    
    edge.df <- merge(dd.df, edge_list, by.x="CELL", by.y="source", all=T)
    edge.df <- edge.df[,c("CELL","state", "x","y","z","target")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "target")
    edge.df <- merge(edge.df, dd.df, by.x="target", by.y="CELL", all=T)
    edge.df <- edge.df[,c("source", "source.state", "source.x", "source.y", "source.z", "target", "state", "x", "y", "z")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "target","target.state","target.x", "target.y", "target.z")
    
    edge.df$STATE=edge.df$source.state
    edge.df$ENTROPY=edge.df$source.z
    
    edge.df.wc <- merge(edge.df, wc.df, by.x="source", by.y="CELL")
    
    edge.df.wc$Cluster=edge.df.wc$WC
    edge.df.wc$Entropy=edge.df.wc$source.z
    g <- ggplot(data=edge.df.wc, aes(x=source.x, y=source.y)) + ggtitle("Cell clusters") + labs(x="PC1", y="PC2") + SLICE_theme_opts();
    g <- g + geom_segment(aes_string(xend="target.x", yend="target.y"), size=I(0.5), linetype="solid", col=I("grey60"), na.rm=TRUE, alpha=.5)
    g <- g + geom_point(aes(colour=Cluster, size=Entropy), na.rm=TRUE)
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    ggsave(filename=paste(context_str, "lineage-g.step4-cell_clusters.pdf", sep=""), width=pdf.w, height=pdf.h) #, dpi=dpi, compression="lzw")
    
    print(g)
  }
  
  
  ################################################## 
  # 2. detect core cells
  ##################################################
  
  # detecting core cells
  
  cells.df <- pData(es)[, c("x","y", dim.cols, "entropy", "slice.state")]
  cells.df$slice.realcell <- 1
  cells.df$slice.stablestate <- "NA"
  cells.df$entropy <- cells.df$entropy/max(cells.df$entropy)
  
  states <- sort(unique(as.character(cells.df$slice.state)))
  
  for (state in states) {
    s.cells <- rownames(cells.df)[which(cells.df$slice.state==state)]
    if (length(s.cells) > 1) {
      s.vids <- which(V(gp)$name %in% s.cells)
      s.gp <- induced.subgraph(gp, s.vids)
      
      s.vertices <- V(s.gp)$name
      
      s.dd <- cells.df[s.vertices, ]
      s.scores <- 1-s.dd$entropy
      names(s.scores) <- rownames(s.dd)
      s.score.threshold <- quantile(as.numeric(s.scores), probs=(1-ss.threshold))
      
      s.scores[which(s.scores<s.score.threshold)] <- 0
      
      s.ss.idx <- which(s.scores>0)
      
      # TODO: minimum number of cells
      
      if (length(s.ss.idx)>0) {
        
        if (ss.method=="pcst") {
          # to use BioNet package, convert igraph object to graphNEL object
          s.gp.2 <- s.gp
          E(s.gp.2)$weight <- E(s.gp)$weight/max(E(s.gp)$weight) # scale the distance to [0,1]
          
          s.gpNEL <- NULL
          s.gpNEL <- igraph.to.graphNEL(s.gp.2)
          
          # library(BioNet)
          s.gpNEL <- rmSelfLoops(s.gpNEL)
          s.module <- runFastHeinz(s.gpNEL, s.scores)
          s.ss <- nodes(s.module)
        } else if (ss.method=="top") {
          s.ss <- names(s.scores)[which(as.numeric(s.scores)>0)]
        } else if (ss.method=="all") {
          s.ss <- names(s.scores)
        }
      } else {
        i.ss <- names(i.scores)
      }
    } else {
      s.ss <- s.cells
    }
    cells.df[s.ss, "slice.stablestate"] <- state
  }
  
  
  states <- sort(unique(as.character(cells.df$slice.state)))
  
  for (state in states) {
    s.ss.cells <- rownames(subset(cells.df, slice.realcell==1 & slice.stablestate==state))
    s.ss.df <- data.frame()
    if (length(s.ss.cells)==1) {
      s.ss.df <- data.frame(cells.df[s.ss.cells, ], check.names=FALSE)
      s.ss.df$slice.realcell <- 0
    } else if (length(s.ss.cells)>1) {
      s.ss.cells.df <- cells.df[which(rownames(cells.df) %in% s.ss.cells), ]
      s.ss.df <- t(data.frame(apply(s.ss.cells.df[, c("x","y",dim.cols,"entropy")], 2, mean), check.names=FALSE))
      s.ss.df <- data.frame(s.ss.df, slice.state=state, slice.realcell=0, slice.stablestate=state, check.names=FALSE)
    }
    rownames(s.ss.df) <- paste("slice.ss.", state, sep="")
    
    cells.df <- rbind(cells.df, s.ss.df)
  }
  
  write.table(cbind(CELL=rownames(cells.df), cells.df), file=paste(context_str, "lineage-g.step4-cell_clusters.txt", sep=""), sep="\t", col.names=T, row.names=F)
  
  
  ss.cells.df <- cells.df[which(cells.df$slice.realcell==0), ]
  if (model.type=="graph") {
    stateGraph <- getStateGraph(as.matrix(t(ss.cells.df[, dim.cols])), type="graph")
  } else { # by default, get mst 
    stateGraph <- getStateGraph(as.matrix(t(ss.cells.df[, dim.cols])), type="mst")
  }
  
  # add pseudo_cells to the network and connect them with stable state cells
  gp.wSSC <- gp    
  min.edge.weight <- min(E(gp)$weight)
  gp.wSSC <- add_vertices(gp.wSSC, dim(ss.cells.df)[1], attr=list(name=rownames(ss.cells.df)))
  for (i in 1:length(states)) {
    i.ss <- rownames(subset(cells.df, slice.stablestate==states[i] & slice.realcell==1))
    i.ssc <- paste("slice.ss.", states[i], sep="")
    if (length(i.ss)>0) {
      for (j in i.ss) {
        i.ssc.vid <- which(V(gp.wSSC)$name %in% i.ssc)
        j.vid <- which(V(gp.wSSC)$name %in% j)
        gp.wSSC <- igraph::add_edges(gp.wSSC, c(i.ssc.vid, j.vid), attr=list(weight=ed(cells.df[which(rownames(cells.df)==i.ssc), c("x","y")], cells.df[which(rownames(cells.df)==j), c("x","y")])))
      }
    }
  }
  
  
  
  # lineage model
  lineageModel <- igraph::as.directed(stateGraph)
  
  edge.df <- as.data.frame(igraph::get.edgelist(lineageModel))  
  edge.df$include <- 1
  for (ei in 1:dim(edge.df)[1]) {
    v1.id <- which(rownames(ss.cells.df) == as.character(edge.df$V1[ei]))
    v2.id <- which(rownames(ss.cells.df) == as.character(edge.df$V2[ei]))
    
    if (reverse==TRUE) {
      if (ss.cells.df$entropy[v1.id] > ss.cells.df$entropy[v2.id]) {
        edge.df$include[ei] <- 0
      }
    } else { # by default, keep entropy reduction path
      if (ss.cells.df$entropy[v1.id] < ss.cells.df$entropy[v2.id]) {
        edge.df$include[ei] <- 0
      }
    }
  }
  todelete <- which(edge.df$include == 0) 
  if (length(todelete)>0) {
    lineageModel <- igraph::delete.edges(lineageModel, todelete)
  }
  
  
  if (do.plot) {
    
    # visualize the similarity network
    edge_list <- as.data.frame(get.edgelist(gp.wSSC))  
    colnames(edge_list) <- c("source", "target")
    
    dd <- cells.df[, c("entropy","x", "y")]
    dd$state <- cells.df$slice.state 
    dd$z <- dd$entropy
    dd$ss <- cells.df$slice.stablestate
    dd$cmm <- cells.df$slice.state
    
    if(is.null(dd$CELL)) {
      dd.df <- cbind(CELL=rownames(dd), dd)
    } else {
      dd.df <- dd
    }
    
    edge.df <- merge(dd.df, edge_list, by.x="CELL", by.y="source", all=T)
    edge.df <- edge.df[,c("CELL","state", "x","y","z","cmm", "ss","target")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "source.cmm", "source.ss", "target")
    edge.df <- merge(edge.df, dd.df, by.x="target", by.y="CELL", all=T)
    edge.df <- edge.df[,c("source", "source.state", "source.x", "source.y", "source.z", "source.cmm", "source.ss", "target", "state", "x", "y", "z", "cmm", "ss")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "source.cmm", "source.ss", "target","target.state","target.x", "target.y", "target.z", "target.cmm", "target.ss")
    
    edge.df$STATE=edge.df$source.state
    edge.df$ENTROPY=edge.df$source.z
    edge.df$STABLE_STATE <- edge.df$source.ss
    edge.df$COMMUNITY <- edge.df$source.cmm
    
    g <- ggplot(data=edge.df, aes(x=source.x, y=source.y)) + ggtitle("Cell clusters and stable states") + labs(x="PC1", y="PC2") + SLICE_theme_opts();
    g <- g + geom_segment(aes_string(xend="target.x", yend="target.y"), size=I(0.5), linetype="solid", col=I("grey60"), na.rm=TRUE, alpha=0.5)
    #g <- g + geom_point(aes(colour=COMMUNITY, size=ENTROPY), na.rm=TRUE, alpha=0.1) 
    
    g <- g + geom_point(data=subset(cells.df, slice.realcell==1), aes(x=x, y=y, col=slice.state, size=entropy), alpha=0.9) 
    g <- g + geom_point(data=cells.df[which(cells.df$slice.realcell==0), ], aes(x=x, y=y, size=entropy), col="grey40")
    #g <- g + geom_point(data=subset(cells.df, slice.realcell==1 & slice.stablestate == "NA" ), aes(x=x, y=y, col=slice.state, size=entropy), alpha=0.2)
    
    
    #g <- g + geom_point(aes(colour=STATE, size=ENTROPY), na.rm=TRUE) # + scale_color_manual(values = cluster_colors)
    
    # g <- g + geom_text(aes_string(label="source"), size=2, hjust=-0.2, vjust=-0.2, na.rm=TRUE)
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    ggsave(filename=paste(context_str, "lineage-g.step5-stable_states.pdf", sep=""), width=pdf.w, height=pdf.h )#, dpi=dpi, compression="lzw")
    
    print(g)
  }
  
  
  if (do.plot) {
    
    # visualize the similarity network
    edge_list <- as.data.frame(get.edgelist(gp.wSSC))  
    colnames(edge_list) <- c("source", "target")
    
    dd <- cells.df[, c("entropy","x", "y")]
    dd$state <- cells.df$slice.state 
    dd$z <- dd$entropy
    dd$ss <- cells.df$slice.stablestate
    dd$cmm <- cells.df$slice.state
    
    if(is.null(dd$CELL)) {
      dd.df <- cbind(CELL=rownames(dd), dd)
    } else {
      dd.df <- dd
    }
    
    edge.df <- merge(dd.df, edge_list, by.x="CELL", by.y="source", all=T)
    edge.df <- edge.df[,c("CELL","state", "x","y","z","cmm", "ss","target")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "source.cmm", "source.ss", "target")
    edge.df <- merge(edge.df, dd.df, by.x="target", by.y="CELL", all=T)
    edge.df <- edge.df[,c("source", "source.state", "source.x", "source.y", "source.z", "source.cmm", "source.ss", "target", "state", "x", "y", "z", "cmm", "ss")]
    colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "source.cmm", "source.ss", "target","target.state","target.x", "target.y", "target.z", "target.cmm", "target.ss")
    
    edge.df$STATE=edge.df$source.state
    edge.df$ENTROPY=edge.df$source.z
    edge.df$STABLE_STATE <- edge.df$source.ss
    edge.df$COMMUNITY <- edge.df$source.cmm
    
    g <- ggplot(data=edge.df, aes(x=source.x, y=source.y)) + ggtitle("Inferred Lineage Model") + labs(x="PC1", y="PC2") + SLICE_theme_opts();
    g <- g + geom_segment(aes_string(xend="target.x", yend="target.y"), size=I(0.5), linetype="solid", col=I("grey60"), na.rm=TRUE, alpha=0.5)
    #g <- g + geom_point(aes(colour=COMMUNITY, size=ENTROPY), na.rm=TRUE, alpha=0.1) 
    
    g <- g + geom_point(data=subset(cells.df, slice.realcell==1), aes(x=x, y=y, col=slice.state, size=entropy)) 
    g <- g + geom_point(data=cells.df[which(cells.df$slice.realcell==0), ], aes(x=x, y=y, size=entropy), col="grey40")
    #g <- g + geom_point(data=subset(cells.df, slice.realcell==1 & slice.stablestate == "NA" ), aes(x=x, y=y, col=slice.state, size=entropy), alpha=0.2)
    
    edge.df2 <- as.data.frame(igraph::get.edgelist(lineageModel))  
    edge.df2$src.x <- edge.df2$src.y <- edge.df2$dst.x <- edge.df2$dst.y <- 0
    for (ei in 1:dim(edge.df2)[1]) {
      src.id <- which(rownames(ss.cells.df) == as.character(edge.df2$V1[ei]))
      dst.id <- which(rownames(ss.cells.df) == as.character(edge.df2$V2[ei]))
      edge.df2$src.x[ei] <- ss.cells.df$x[src.id]
      edge.df2$src.y[ei] <- ss.cells.df$y[src.id]
      edge.df2$dst.x[ei] <- ss.cells.df$x[dst.id]
      edge.df2$dst.y[ei] <- ss.cells.df$y[dst.id]
    }
    g <- g + geom_segment(data=edge.df2, aes(x=src.x, y=src.y, xend=dst.x, yend=dst.y), size=I(2), linetype="solid", col=I("black"), alpha=0.6, arrow=arrow(), na.rm=TRUE)
    
    #g <- g + geom_point(aes(colour=STATE, size=ENTROPY), na.rm=TRUE) # + scale_color_manual(values = cluster_colors)
    
    # g <- g + geom_text(aes_string(label="source"), size=2, hjust=-0.2, vjust=-0.2, na.rm=TRUE)
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    ggsave(filename=paste(context_str, "lineage-g.step6-lineage_model.pdf", sep=""), width=pdf.w, height=pdf.h )#, dpi=dpi, compression="lzw")
    
    print(g)
    
  }
  
  # lineage model
  ret <- list(mst=mst, csn=gp, csn.wSS= gp.wSSC, cells.df=cells.df, stateGraph=stateGraph, lineageModel=lineageModel)
  # class(ret) <- "slice.gmodel"
  return(ret)
}


get.SP.Transitions <- function(es, model, start, end,  network="mst", NN.threshold=0.8, NN.type="all", do.plot=T, context_str="") {
  
  if (is.null(model$lineageModel) | length(V(model$lineageModel))<2 ) {
    stop("Not enough states")
  }
  
  validStates <- 1:length(V(model$lineageModel))
  if (!(start %in% validStates)) {
    stop("Invalid start")
  }
  if (!(end %in% validStates)) {
    stop("Invalid end state")
  }
  
  g <- model$lineageModel
  paths <- igraph::all_simple_paths(g, from=V(g)[start], to=V(g)[end], mode="out")
  
  transition.paths <- list()
  cells.df <- model$cells.df
  dim.cols.idx <- which(substr(colnames(cells.df), 1, 3)=="Dim")
  if (length(dim.cols.idx)<2) {
    stop("At least two columns with prefix \"Dim\" in pData(es) are required.")
  }
  dim.cols <- colnames(cells.df)[dim.cols.idx]
  
  cat("\nReconstructing shortest-path based cell trajectories following the inferred lineage C", start, "->C", end, "\n", sep="")
  
  if (network=="mst") {
    net <- getCellSimilarityNetwork(as.matrix(t(cells.df[, dim.cols])), wiring.threshold = NULL)
    net <- net$mst
  } else if (network == "csn") {
    net <- model$csn.wSS
  } else {
    stop("Invalid network")
  }
  
  transitions <- list()
  transitions.wNN <- list()
  NN.dist <- as.numeric(quantile(E(net)$weight, probs=NN.threshold))
  # NN.dist <- as.numeric(max(E(net)$weight)*NN.threshold)
  
  if (length(paths)<=0) {
    warning(paste("no transition path from ", V(g)$name[start], " to ", V(g)$name[end], sep=""))
  } else {
    for (i in 1:length(paths)) {
      i.start <- as.character(V(g)$name[start])
      i.end <- as.character(V(g)$name[end])
      i.path <- paths[[i]]$name
      i.path.states <- cells.df$slice.state[which(rownames(cells.df) %in% i.path)]
      
      
      j <- 1
      k <- j+1
      
      jk.start <- i.path[j]
      jk.end <- i.path[k]
      
      jk.transition <- igraph::get.shortest.paths(net, from=jk.start, to=jk.end)
      
      i.transition.df <- data.frame(CELL=as.character(jk.transition$vpath[[1]]$name), check.names=FALSE)
      i.transition.df$PTIME <- 1:dim(i.transition.df)[1]
      i.transition.df$CellType <- "real"
      i.transition.df$CellType[which(substr(i.transition.df$CELL, 1,9)=="slice.ss.")] <- "pseudo"
      
      j <- j+1
      k <- k+1
      
      while(k<=length(i.path)) {
        
        jk.start <- i.path[j]
        jk.end <- i.path[k]
        
        temp <- igraph::get.shortest.paths(net, from=jk.start, to=jk.end)
        
        temp.df <- data.frame(CELL=as.character(temp$vpath[[1]]$name[-1]), check.names=FALSE)
        temp.df$PTIME <- (dim(i.transition.df)[1]+1):(dim(i.transition.df)[1]+dim(temp.df)[1])
        temp.df$CellType <- "real"
        temp.df$CellType[which(substr(temp.df$CELL, 1,9)=="slice.ss.")] <- "pseudo"
        
        i.transition.df <- rbind(i.transition.df, temp.df)
        
        j <- j+1
        k <- k+1
      }
      
      
      i.transition.wNN.df <- i.transition.df
      i.cells.n <- dim(i.transition.df)[1]
      for (k in 1:i.cells.n) {
        k.cell <- i.transition.df$CELL[k]
        k.cell.vid <- which(V(net)$name %in% k.cell)
        if (NN.type == "NN") {
          k.neighbors <- igraph::neighbors(net, k.cell.vid)$name
        } else {
          k.neighbors <- igraph::V(net)$name
        }
        
        k.nn <- c()
        for (m in k.neighbors) {
          if (!(m==k.cell)) {
            m.k.ed <- Inf
            m.k.ed <- ed(cells.df[which(rownames(cells.df)==m), c("x","y")], cells.df[which(rownames(cells.df)==k.cell), c("x","y")])
            if (m.k.ed <= NN.dist) {
              k.nn <- c(k.nn, m)
            }
          }
          
        }
        if (length(k.nn) > 0) {
          k.df <- NULL
          k.df <- data.frame(CELL=k.nn, PTIME=i.transition.df$PTIME[k], CellType="real")
          i.transition.wNN.df <- rbind(i.transition.wNN.df, k.df)
        }
      }
      
      i.transition.wNN.df$CellType[which(substr(i.transition.wNN.df$CELL, 1, 9)=="slice.ss.")] <- "pseudo"
      
      #i.transition.wNN.df <- i.transition.wNN.df[which(i.transition.wNN.df$CellType == "real"), ]
      
      i.transition.wNN.df <- i.transition.wNN.df[order(i.transition.wNN.df$PTIME, decreasing=FALSE), ]
      # i.transition.wNN.df$CellType <- "real"
      
      i.pseudotime <- i.transition.wNN.df[which(i.transition.wNN.df$CellType == "real"), c("CELL","PTIME")]
      
      
      
      transitions[[i]] <- list(start=i.start, end=i.end, path=i.path,
                               i.path.states=i.path.states,
                               i.transition= i.transition.df,
                               i.transition.wNN = i.transition.wNN.df,
                               i.pseudotime = i.pseudotime
      )
      names(transitions)[i] <- paste(i.start, "_to_", i.end, "_path", i, sep="")
      
      if (do.plot) {
        
        pdf.w <- 8
        pdf.h <- 6
        
        # visualize the transitions
        edge_list <- as.data.frame(get.edgelist(net))  
        colnames(edge_list) <- c("source", "target")
        
        dd <- cells.df[, c("entropy","x", "y")]
        dd$state <- cells.df$slice.state
        dd$ss <- cells.df$slice.stablestate
        dd$cmm <- cells.df$slice.state
        dd$z <- dd$entropy
        if(is.null(dd$CELL)) {
          dd.df <- cbind(CELL=rownames(dd), dd)
        } else {
          dd.df <- dd
        }
        
        edge.df <- merge(dd.df, edge_list, by.x="CELL", by.y="source", all=T)
        edge.df <- edge.df[,c("CELL","state", "x","y","z", "cmm", "ss", "target")]
        colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "source.cmm", "source.ss", "target")
        edge.df <- merge(edge.df, dd.df, by.x="target", by.y="CELL", all=T)
        edge.df <- edge.df[,c("source", "source.state", "source.x", "source.y", "source.z", "source.cmm", "source.ss", "target", "state", "x", "y", "z", "cmm", "ss")]
        colnames(edge.df) <- c("source", "source.state", "source.x", "source.y", "source.z", "source.cmm", "source.ss", "target","target.state","target.x", "target.y", "target.z", "target.cmm", "target.ss")
        
        edge.df$STATE=edge.df$source.state
        edge.df$Entropy=edge.df$source.z
        edge.df$Cluster <- factor(edge.df$source.cmm)
        edge.df$STABLE_STATE <- edge.df$source.ss
        
        g <- ggplot(data=edge.df, aes(x=source.x, y=source.y)) + ggtitle(paste("Transition: ", i.start, " -- ", i.end, sep="")) + labs(x="PC1", y="PC2") + SLICE_theme_opts();
        g <- g + geom_segment(aes_string(xend="target.x", yend="target.y"), size=I(0.5), linetype="solid", col=I("grey60"), na.rm=TRUE, alpha=0.5)
        g <- g + geom_point(aes(colour=Cluster, size=Entropy), na.rm=TRUE) + scale_shape_manual(values = c(16, 17, 15, 18, 21, 22, 23, 24, 25,  1,2,0,5,6,7,11,14))
        # g <- g + geom_point(aes(colour=STATE, size=ENTROPY), na.rm=TRUE) # + scale_color_manual(values = cluster_colors)
        # g <- g + geom_point(data=edge.df[which(substr(edge.df$STATE, 1, 4)=="ssc_"),], aes(colour=STABLE_STATE)) + scale_color_manual(values="black")
        g <- g + geom_point(data=cells.df[which(cells.df$slice.realcell==0), ], aes(x=x, y=y, size=entropy), col="grey30")
        
        i.transition.df <- merge(i.transition.df, dd.df, by.x="CELL", by.y="CELL")
        i.transition.df <- i.transition.df[order(i.transition.df$PTIME),]
        
        g <- g + geom_path(aes(x=x, y=y), alpha=I(0.5), color=I("orange3"), size=I(5), data=i.transition.df, na.rm=TRUE) #+ geom_point(aes(x=x, y=y, color=CLUSTER), size=I(2), data=mfb.df, na.rm=TRUE)
        
        for (j in 1:dim(i.transition.df)[1]) {
          dat <- circleFun(as.numeric(i.transition.df[j,c("x", "y")]), NN.dist)
          # g <- g + geom_path(data=dat, aes(x=x, y=y), col=I('limegreen'), size=I(0.8), linetype=2, alpha=0.5)
        }
        g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
        
        ggsave(filename=paste(context_str, "path-sp-", names(transitions)[i], ".pdf", sep=""), width=pdf.w, height=pdf.h) #, dpi=dpi, compression="lzw")
        
      }
    }
  }
  
  ret <- list(net=net, cells.df=cells.df, transitions=transitions)
  #class(ret) <- "slice.sp.transitions"
  return(ret)
}


get.SPT.Profiles <- function(es, transitions, context_str="", 
                             do.plot=T, plot.xlabel=NULL, plot.ylabel=NULL, plot.w=2.8, plot.h=1) {    
    
  cat("\nExtracting lineage dependent gene expression profiles\n")
  
  transition.profiles <- list()
  
  for (i in 1:length(transitions)) {
    
    # generate expression data along the transition
    i.transition <- transitions[[i]]
    i.pseudotime <- i.transition$i.pseudotime
    i.es <- es[, as.character(i.pseudotime$CELL)]
    i.exp <- as.data.frame(exprs(i.es), check.names=FALSE)
    colnames(i.exp) <- paste("T", i.pseudotime$PTIME, "_", i.pseudotime$CELL, sep="")
    
    i.gexp <- data.frame(CELL=i.pseudotime$CELL, PTIME=i.pseudotime$PTIME, t(i.exp), check.names=FALSE)
    
    transition.profiles[[i]] <- list(exprs=i.exp, cells=colnames(i.exp), genes=rownames(i.exp), ptime=i.pseudotime$PTIME)
    names(transition.profiles)[i] <- names(transitions)[i]
    
    # write.table(data.frame(RID=rownames(i.gexp), i.gexp, check.names=FALSE), file=paste(context_str, "lineage-sp-profiles-", names(transitions)[i], ".txt", sep=""), sep="\t", col.names=T, row.names=F)
    
    if(is.null(plot.xlabel)) {
      plot.xlabel = names(transitions)[i]
    }
    if (is.null(plot.ylabel)) {
      plot.ylabel = "Expression (+1, log2)"
    }
    
    if (do.plot) {
      
      if (dim(fData(es))[1] > 20) {
        
        warning("Too many number of profiles (>20) for visualization. Plots will not be generated!")
        
      } else {
        i.exp <- log(i.exp+1,2)
        
        dd <- data.frame(gene=rownames(i.exp), i.exp, check.names=FALSE)
        cell.ptime <- data.frame(cell=colnames(i.exp), ptime=i.pseudotime$PTIME, check.names = FALSE)
        
        # cell.ptime$ptime <- (cell.ptime$ptime-min(cell.ptime$ptime))/(max(cell.ptime$ptime)-min(cell.ptime$ptime))
        
        dd.melt <- melt(dd, id.vars = c("gene"))
        dd.melt <- merge(dd.melt, cell.ptime, by.x="variable", by.y="cell")
        
        dd.melt$gene <- factor(dd.melt$gene, levels = as.character(rownames(fData(es))))
        
        i.exp.avg <- aggregate(x=t(i.exp), by=list(PTIME=i.pseudotime$PTIME), FUN="mean")
        i.exp.avg <- i.exp.avg[, 2:dim(i.exp.avg)[2]]
        i.exp.avg <- t(i.exp.avg)
        colnames(i.exp.avg) <- paste("v", 1:dim(i.exp.avg)[2], sep="")
        cell.ptime.2 <- data.frame(cell=colnames(i.exp.avg), ptime=1:dim(i.exp.avg)[2])
        dd.avg <- data.frame(gene=rownames(i.exp.avg), i.exp.avg, check.names=FALSE)
        dd.avg.melt <- melt(dd.avg, id.vars=c("gene"))
        dd.avg.melt <- merge(dd.avg.melt, cell.ptime.2, by.x="variable", by.y="cell")
        
        fmt<-function(x){format(x,nsmall = 2,scientific = FALSE)}
        
        g <- ggplot(dd.melt) + facet_wrap("gene", ncol=1, scales="free_y") 
        
        if (FALSE) {
          g <- g + stat_summary(aes(x = ptime, y = value, group = 1), 
                                fun.data = mean_cl_normal, fun.args = list(mult = 1), color = "grey80", fill = "grey80", 
                                size = I(1), geom = c("ribbon"))
        }
        g <- g + stat_summary(aes(x = ptime, y = value, group = 1), 
                              fun.data = mean_cl_normal, fun.args = list(mult = 1), color = "grey40", fill = "grey40", 
                              size = I(0.8), geom = c("point"))
        
        g <- g + stat_smooth(data=dd.avg.melt, aes(x=ptime, y=value, group=1), method="loess", colour="royalblue4", se=FALSE)
        
        g <- g + scale_y_continuous(labels=fmt, limits=c(0, NA))
        
        g <- g + theme(axis.text.x = element_text(angle = 0, hjust = 0)) + xlab(plot.xlabel) + ylab(plot.ylabel)
        
        g <- g + SLICE_theme_opts()
        
        g <- g + theme(strip.background = element_rect(colour = 'white', fill = 'white')) + 
          theme(panel.border = element_blank(), axis.line.x = element_line(size=0.2), axis.line.y = element_line(size=0.2)) +
          theme(axis.ticks = element_line(size = 0.2)) +
          theme(legend.position="none") +
          theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
          theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
        
        n <- dim(i.exp.avg)[1]
        
        fig.filename <- paste(context_str, "path-sp-", names(transitions)[i], "-profiles.pdf", sep="")
        
        ggsave(fig.filename, width=plot.w, height=plot.h*n)
        
        print(g)
        
      }
    }
    
  }
  
  #class(transition.profiles) <- "slice.sp.profiles"
  return(transition.profiles)
}


getLM.clustering <- function(es, model.type="tree", cluster.method="kmeans", 
                             k=NULL, k.max=10, B=100, k.opt.method="firstmax", 
                             ss.method="top", ss.threshold=0.25,                                   
                             reverse = FALSE, 
                             do.plot=T, context_str="") {
  
  cat("\nInferring lineage model using clustering-based method\n")
  
  # if (is.null(rds.cols)) {
  dim.cols.idx <- which(substr(colnames(pData(es)), 1, 3)=="Dim")
  if (length(dim.cols.idx)<2) {
    stop("At least two columns with prefix \"Dim\" in pData(es) are required.")
  }
  dim.cols <- colnames(pData(es))[dim.cols.idx]
  # }
  
  if (do.plot) { 
    pdf.w <- 8
    pdf.h <- 6
    
    dd <- pData(es)[, c("state", "entropy","x", "y")]
    dd$Cluster <- dd$state
    dd$Entropy <- dd$entropy/max(dd$entropy)
    dd$z <- dd$entropy
    
    g <- ggplot(data=dd, aes(x=x, y=y)) + geom_point(aes(color=Cluster, size=Entropy)) + SLICE_theme_opts() + ggtitle("Cells") + labs(x="PC1", y="PC2") #+ geom_text(label=rownames(dd), size=2, hjust=-0.2, vjust=-0.2) ;
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    ggsave(filename=paste(context_str, "lineage-c.step1-cells.pdf", sep=""), width=pdf.w, height=pdf.h) 
    print(g)
  }
  
  
  ################################################## 
  # 1. detecting cell states
  ##################################################
  
  # if (do.clustering==TRUE) {
  
  reducedDims <- (pData(es)[, dim.cols.idx])
  
  if (is.null(k)) {
    if (cluster.method=="pam") {
      gskmn <- cluster::clusGap(as.matrix(reducedDims), FUN = pam,  K.max = k.max, B = B)
    } else if (cluster.method=="kmeans") {
      gskmn <- cluster::clusGap(as.matrix(reducedDims), FUN = kmeans,  K.max = k.max, B = B)
    }
    plot(gskmn, main = paste("Number of Cellular States determined by\n Gap statistic (FUN=",cluster.method,", k.max=", k.max, ", B=",B, ")", sep=""))
    
    # default: firstmax rule
    k <- maxSE(f=as.numeric(gskmn$Tab[,"gap"]), SE.f=as.numeric(gskmn$Tab[,"SE.sim"]), method=k.opt.method)
    abline(v=k, col="blue", lty=2, lwd=3)
    text(k, min(gskmn$Tab[,"gap"]), paste(" k=",k, "\n(", k.opt.method, ")", sep=""), col="blue", adj = c(-.1, -.1), cex=1.1)
  }
  
  if (cluster.method == "pam") {
    db <- pam(as.matrix(reducedDims), k)
    pData(es)$slice.state <- paste("C",as.numeric(db$clustering),sep="")
  } else if (cluster.method=="kmeans") {
    db <- kmeans(as.matrix(reducedDims), k)
    pData(es)$slice.state <- paste("C",as.numeric(db$cluster),sep="")
  }
  
  # }
  
  
  if (do.plot) {
    pdf.w <- 8
    pdf.h <- 6
    
    dd <- pData(es)[, c("slice.state", "entropy","x", "y")]
    dd$color <- dd$slice.state
    dd$z <- dd$entropy
    
    g <- ggplot(data=dd, aes(x=x, y=y)) + geom_point(aes(color=color, size=entropy)) + SLICE_theme_opts() + ggtitle("Cell clusters") + labs(x="PC1", y="PC2") #+ geom_text(label=rownames(dd), size=2, hjust=-0.2, vjust=-0.2) ;
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    ggsave(filename=paste(context_str, "lineage-c.step2-cell_clusters.pdf", sep=""), width=pdf.w, height=pdf.h) 
    print(g)
  }
  
  
  # detect core cells and stable states
  
  cells.df <- pData(es)[, c("x","y", dim.cols, "entropy","slice.state")]
  cells.df$slice.realcell <- 1
  cells.df$slice.stablestate <- "NA"
  cells.df$entropy <- cells.df$entropy/max(cells.df$entropy)
  
  
  states <- sort(unique(as.character(cells.df$slice.state)))
  for (state in states) {
    s.cells <- rownames(cells.df)[which(cells.df$slice.state == state)]
    s.entropies <- cells.df$entropy[which(cells.df$slice.state == state)]
    
    s.entropies.thresh <- quantile(s.entropies, probs=(ss.threshold))
    
    s.ss.idx <- which(s.entropies <= s.entropies.thresh)
    
    if (ss.method=="pcst") {
      # TODO
    }
    
    s.ss.cells <- s.cells[s.ss.idx]
    
    cells.df$slice.stablestate[which(rownames(cells.df) %in% s.ss.cells)] <- state
    
    s.ss.cells.df <- cells.df[which(cells.df$slice.stablestate == state), ]
    
    s.ss.df <- t(data.frame(apply(s.ss.cells.df[, c("x","y",dim.cols,"entropy")], 2, mean)))
    rownames(s.ss.df) <- paste("slice.ss.", state, sep="")
    s.ss.df <- data.frame(s.ss.df, slice.state=state, slice.realcell=0, slice.stablestate=state)
    
    cells.df <- rbind(cells.df, s.ss.df)
    
  }
  
  if (do.plot) {
    g <- ggplot() + ggtitle("Cell clusters and stable states") + labs(x="PC1", y="PC2")
    g <- g + geom_point(data=subset(cells.df, slice.realcell==1), aes(x=x, y=y, col=slice.state, size=entropy)) 
    #g <- g + geom_point(data=subset(cells.df, slice.realcell==1 & slice.stablestate == "NA" ), aes(x=x, y=y, col=slice.state, size=entropy))#, alpha=0.2)
	g <- g + geom_point(data=cells.df[which(cells.df$slice.realcell==0), ], aes(x=x, y=y, size=entropy), col="black")
    g <- g + SLICE_theme_opts()
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    ggsave(filename=paste(context_str, "lineage-c.step3-stable_states.pdf", sep=""), width=pdf.w, height=pdf.h) 
    
    print(g)
  }
  
  ss.cells.df <- cells.df[which(cells.df$slice.realcell==0), ]
  if (model.type=="graph") {
    stateGraph <- getStateGraph(as.matrix(t(ss.cells.df[, dim.cols])), type="graph")
  } else { # by default, get mst 
    stateGraph <- getStateGraph(as.matrix(t(ss.cells.df[, dim.cols])), type="mst")
  }
  
  
  write.table(cbind(CELL=rownames(cells.df), cells.df), file=paste(context_str, "lineage-c.step3-cell_clusters.txt", sep=""), sep="\t", col.names=T, row.names=F)
  
  
  # lineage model
  lineageModel <- igraph::as.directed(stateGraph)
  
  edge.df <- as.data.frame(igraph::get.edgelist(lineageModel))  
  edge.df$include <- 1
  for (ei in 1:dim(edge.df)[1]) {
    v1.id <- which(rownames(ss.cells.df) == as.character(edge.df$V1[ei]))
    v2.id <- which(rownames(ss.cells.df) == as.character(edge.df$V2[ei]))
    
    if (reverse==TRUE) {
      if (ss.cells.df$entropy[v1.id] > ss.cells.df$entropy[v2.id]) {
        edge.df$include[ei] <- 0
      }
    } else { # by default, keep entropy reduction path
      if (ss.cells.df$entropy[v1.id] < ss.cells.df$entropy[v2.id]) {
        edge.df$include[ei] <- 0
      }
    }
  }
  todelete <- which(edge.df$include == 0) 
  if (length(todelete)>0) {
    lineageModel <- igraph::delete.edges(lineageModel, todelete)
  }
  
  
  if (do.plot) {
    g <- ggplot() + ggtitle("Inferred Lineage Model") + labs(x="PC1", y="PC2")
    g <- g + geom_point(data=subset(cells.df, slice.realcell==1 & slice.stablestate != "NA" ), aes(x=x, y=y, col=slice.state, size=entropy)) 
    g <- g + geom_point(data=cells.df[which(cells.df$slice.realcell==0), ], aes(x=x, y=y, size=entropy), col="black")
    g <- g + geom_point(data=subset(cells.df, slice.realcell==1 & slice.stablestate == "NA" ), aes(x=x, y=y, col=slice.state, size=entropy))
    
    edge.df <- as.data.frame(igraph::get.edgelist(lineageModel))  
    edge.df$src.x <- edge.df$src.y <- edge.df$dst.x <- edge.df$dst.y <- 0
    for (ei in 1:dim(edge.df)[1]) {
      src.id <- which(rownames(ss.cells.df) == as.character(edge.df$V1[ei]))
      dst.id <- which(rownames(ss.cells.df) == as.character(edge.df$V2[ei]))
      edge.df$src.x[ei] <- ss.cells.df$x[src.id]
      edge.df$src.y[ei] <- ss.cells.df$y[src.id]
      edge.df$dst.x[ei] <- ss.cells.df$x[dst.id]
      edge.df$dst.y[ei] <- ss.cells.df$y[dst.id]
    }
    g <- g + geom_segment(data=edge.df, aes(x=src.x, y=src.y, xend=dst.x, yend=dst.y), size=I(2), linetype="solid", col=I("black"), alpha=0.6, arrow=arrow(), na.rm=TRUE)
    
    
    g <- g + SLICE_theme_opts()
    g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
    
    
    ggsave(filename=paste(context_str, "lineage-c.step4-lineage_model.pdf", sep=""), width=pdf.w, height=pdf.h) 
    
    print(g)
    
  }
  
  ret <- list(es=es, cells.df=cells.df, stateGraph=stateGraph, lineageModel=lineageModel)
  class(ret) <- "slice.cmodel"
  return(ret)
}


get.PC.Transitions <- function(es, model, start, end, do.trim=F, do.plot=T, context_str="") {
  
  if (is.null(model$lineageModel) | length(V(model$lineageModel))<2 ) {
    stop("Not enough states")
  }
  
  validStates <- 1:length(V(model$lineageModel))
  if (!(start %in% validStates)) {
    stop("Invalid start")
  }
  if (!(end %in% validStates)) {
    stop("Invalid end state")
  }
  
  g <- model$lineageModel
  paths <- igraph::all_simple_paths(g, from=V(g)[start], to=V(g)[end], mode="out")
  
  transition.paths <- list()
  cells.df <- model$cells.df
  dim.cols.idx <- which(substr(colnames(cells.df), 1, 3)=="Dim")
  if (length(dim.cols.idx)<2) {
    stop("At least two columns with prefix \"Dim\" in pData(es) are required.")
  }
  dim.cols <- colnames(cells.df)[dim.cols.idx]
  
  cat("\nReconstructing principal-curve based cell trajectories following the inferred lineage C", start, "->C", end, "\n", sep="")
    
  if (length(paths)<=0) {
    warning(paste("no transition path from ", V(g)$name[start], " to ", V(g)$name[end], sep=""))
  } else {
    for (i in 1:length(paths)) {
      i.start <- as.character(V(g)$name[start])
      i.end <- as.character(V(g)$name[end])
      i.path <- paths[[i]]$name
      i.path.states <- cells.df$slice.state[which(rownames(cells.df) %in% i.path)]
      i.cells.df <- cells.df[which(cells.df$slice.state %in% i.path.states), ]
      i.cells.w <- (1-as.numeric(i.cells.df$entropy))
      i.pcurve.fit <- principal.curve(as.matrix(i.cells.df[, dim.cols]), smoother="smooth.spline", w=i.cells.w, df=length(i.path)+1)
      
      i.pseudotime <- as.data.frame(i.pcurve.fit$s[i.pcurve.fit$tag, ]) # correct direction
      s.id <- which(rownames(i.pseudotime)==i.start)
      e.id <- which(rownames(i.pseudotime)==i.end)
      if (s.id > e.id) { 
        if (do.trim==T) {
          i.pseudotime <- i.pseudotime[e.id:s.id, ]
        }
        i.pseudotime <- i.pseudotime[dim(i.pseudotime)[1]:1, ]
      } else {
        if (do.trim==T) {
          i.pseudotime <- i.pseudotime[s.id:e.id, ]
        }
      }
      i.pseudotime <- i.pseudotime[which(!(rownames(i.pseudotime)) %in% i.path),] # remove pseudotime cells
      i.pseudotime$ptime <- seq(0, 1, length.out = dim(i.pseudotime)[1]) # add uniform pseudotime
      
      transition.paths[[i]] <- list(start=i.start, end=i.end, path=i.path,
                                    i.path.states=i.path.states,
                                    i.cells.df=i.cells.df, i.cells.w=i.cells.w,
                                    i.pcurve.fit=i.pcurve.fit, 
                                    i.pseudotime=i.pseudotime)
      names(transition.paths)[i] <- paste(i.start, "_to_", i.end, "_path", i, sep="")
    }
    
    if (do.plot) {
      pdf.w <- 8
      pdf.h <- 6
      
      for (i in 1:length(transition.paths)) {
        g <- ggplot() + ggtitle("Cell transitional path") + labs(x="PC1", y="PC2")
        g <- g + ggtitle(names(transition.paths)[i])
        g <- g + geom_point(data=subset(cells.df, slice.realcell==1 & slice.stablestate != "NA" ), aes(x=x, y=y, col=slice.state, size=entropy)) 
        g <- g + geom_point(data=cells.df[which(cells.df$slice.realcell==0), ], aes(x=x, y=y, size=entropy), col="black")
        g <- g + geom_point(data=subset(cells.df, slice.realcell==1 & slice.stablestate == "NA" ), aes(x=x, y=y, col=slice.state, size=entropy))
        
        g <- g + geom_path(data=data.frame(transition.paths[[i]]$i.pseudotime), aes(x=Dim1, y=Dim2), alpha=0.6, col="orange3", size=4, arrow=arrow(ends="last", type = "closed", length = unit(0.25, "inches")))
        
        g <- g + SLICE_theme_opts()
        g <- g + theme(axis.line.x = element_line(size=0.5), axis.line.y=element_line(size=0.5)) 
        
        
        ggsave(filename=paste(context_str, "path-pc-", names(transition.paths)[i], ".pdf", sep=""), width=pdf.w, height=pdf.h) #, dpi=dpi, compression="lzw")
        
        print(g)
       
      }
    }
  }
  
  #class(transition.paths) <- "slice.pc.transitions"
  return(transition.paths)
}


get.PCT.Profiles <- function(es, transitions, context_str="", 
                             do.plot=T, plot.xlabel=NULL, plot.ylabel=NULL, plot.w=2.8, plot.h=1) {
  
  cat("\nExtracting lineage dependent gene expression profiles\n")
  
  if(length(transitions)<=0) {
    stop("Empty transitions")
  }
  
  transition.profiles <- list()
  
  for (i in 1:length(transitions)) {
    
    # generate expression data along the transition
    i.transition <- transitions[[i]]
    
    i.pseudotime <- i.transition$i.pseudotime
    
    i.es <- es[, as.character(rownames(i.pseudotime))]
    i.exp <- as.data.frame(exprs(i.es), check.names=FALSE)
    #colnames(i.exp) <- paste("T", i.transition$PTIME, "_", i.transition$CELL, sep="")
    
    i.gexp <- data.frame(CELL=as.character(rownames(i.pseudotime)), PTIME=i.pseudotime$ptime, t(i.exp), check.names=FALSE)
    
    transition.profiles[[i]] <- list(exprs=i.exp, cells=colnames(i.exp), genes=rownames(i.exp), ptime=i.pseudotime$ptime)
    names(transition.profiles)[i] <- names(transitions)[i]
    
    #write.table(data.frame(RID=rownames(i.gexp), i.gexp), file=paste(context_str, "lineage-trajectory.", names(transitions)[i], ".profiles.txt", sep=""), sep="\t", col.names=T, row.names=F)
    
    if(is.null(plot.xlabel)) {
      plot.xlabel = names(transitions)[i]
    }
    if (is.null(plot.ylabel)) {
      plot.ylabel = "Expression (+1, log2)"
    }
    
    if (do.plot) {
      
      if (dim(fData(es))[1] > 20) {
        
        warning("Too many number of profiles (>20) for visualization. Plots will not be generated!")
        
      } else {
        
        
        dd <- data.frame(gene=rownames(i.exp), log(i.exp+1,2), check.names=FALSE)
        cell.ptime <- data.frame(cell=colnames(i.exp), ptime=i.pseudotime$ptime, check.names=FALSE)
        
        dd.melt <- melt(dd, id.vars = c("gene"))
        dd.melt <- merge(dd.melt, cell.ptime, by.x="variable", by.y="cell")
        
        dd.melt$gene <- factor(dd.melt$gene, levels = as.character(rownames(i.es)))
        
        fmt<-function(x){format(x,nsmall = 2,scientific = FALSE)}
        
        g <- ggplot(dd.melt) + facet_wrap("gene", ncol=1, scales="free_y") 
        
        
        g <- g + stat_summary(aes(x = ptime, y = value, group = 1), 
                              fun.data = mean_cl_normal, fun.args = list(mult = 1), color = "grey40", fill = "grey40", 
                              size = I(1), geom = c("point"))
        
        g <- g + stat_smooth(data=dd.melt, aes(x=ptime, y=value, group=1), method="lm", formula=y~ns(x, 3), colour="royalblue4", se=FALSE)
        
        
        g <- g + scale_y_continuous(labels=fmt, limits=c(0, NA))
        
        g <- g + theme(axis.text.x = element_text(angle = 0, hjust = 0)) + xlab(plot.xlabel) + ylab(plot.ylabel)
        
        g <- g + SLICE_theme_opts()
        
        g <- g + theme(strip.background = element_rect(colour = 'white', fill = 'white')) + 
          theme(panel.border = element_blank(), axis.line.x = element_line(size=0.2), axis.line.y = element_line(size=0.2)) +
          theme(axis.ticks = element_line(size = 0.2)) +
          theme(legend.position="none") +
          theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
          theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
        
        print(g)
        
        n <- dim(i.exp)[1]
        
        fig.filename <- paste(context_str, "path-pc-", names(transitions)[i], "-profiles.pdf", sep="")
        
        ggsave(fig.filename, width=plot.w, height=plot.h*n)
        
      }
    }
    
  }
  
  #class(transition.profiles) <- "slice.pc.profiles"
  return(transition.profiles)
}


fitSmoothModel <- function(x, y=NULL, k=NULL, model=1) {
  
  if (is.null(y)) {
    y <- seq(0,1,length.out=length(x))
  } 
  if (is.null(k)) {
    k <- length(unique(y))
  }
  
  dd <- data.frame(V1=as.numeric(x), V2=y)
  if (model==1) {
    m <- mgcv::gam(V1 ~ s(V2, k=4, bs = "cr"), data = dd, method = "REML")
  } else {
    m <- mgcv::gam(V1 ~ 1, data=dd, method="REML")
  }
  return(m)
}


fit_helper <- function(x, ptime, k) {
  
  #dd <- data.frame(V=as.numeric(x), ptime=ptime)
  if (var(as.numeric(x))==0) {
    pvalue <- 1
    dd <- data.frame(V=as.numeric(x), ptime=ptime)
    dd <- unique(dd)
    m1.expr <- as.numeric(dd$V)
  } else {
    m1 <- fitSmoothModel(x=as.numeric(x), y=ptime, k=k, model=1)
    m0 <- fitSmoothModel(x=as.numeric(x), y=ptime, k=k, model=0)
    lr.ret <- lmtest::lrtest(m1, m0)
    pvalue <- lr.ret[2, "Pr(>Chisq)"]
    dd <- data.frame(V=as.numeric(predict(m1)), ptime=ptime)
    dd <- unique(dd)
    m1.expr <- as.numeric(dd$V)
  }
  
  return (c(pvalue, m1.expr))
}


getCellSimilarityNetwork <- function(x, dist.method="euclidean", wiring.threshold=function(x) max(x)*1) {
  
  
  #######################################################################
  ##############           input validation            ##################
  #######################################################################
  
  dms <- c("euclidean", "pearson", "spearman")
  d.id <- pmatch(dist.method, dms)
  if (is.na(d.id)) {
    stop("invalid distance measure.")
  }
  dist.method <- dms[d.id]
  
  #######################################################################
  ##############         cell pairwise distance        ##################
  #######################################################################
  
  # cell pairwise distance
  A <- NULL
  
  if (class(x)=="dist") {
    A <- as.matrix(x)
  } else if (class(x)=="matrix") {
    if (dist.method=="pearson" | dist.method=="spearman") {
      A <- (1-cor(x, method=dist.method))/2
    } else {
      cdist <- dist(t(x), method=dist.method)
      A <- as.matrix(cdist)
    }
  } else {
    stop("invalid format. x shall be of class \"dist\" or \"matrix\".")
  }
  
  #######################################################################
  ##############  constructing cell similarity graph   ##################
  #######################################################################
  
  # minimum spanning tree
  mst <- NULL 
  mst <- graph.adjacency(A, mode="undirected", weighted=TRUE)
  mst <- minimum.spanning.tree(mst)
  
  csn <- mst
  
  if (!is.null(wiring.threshold)) {
    threshold <- wiring.threshold(E(csn)$weight)
    
    if (threshold>0) {
      # local wiring
      AE <- melt(A)
      colnames(AE) <- c("source", "target", "weight")
      AE <- AE[which(AE$weight<=threshold & AE$weight>0),]
      for (i in 1:dim(AE)[1]) {
        v1 <- as.character(AE$source[i])
        v2 <- as.character(AE$target[i])
        w <- as.numeric(AE$weight[i])
        if(are.connected(csn, v1, v2) == FALSE) {
          csn <- add.edges(csn, c(v1, v2), attr=list(weight=w))
        }
      }
    }
  }
  
  return(list(mst=mst, csn=csn))
}


getCellCommunities <- function(csn, method="louvain") {
  
  ms <- c("auto", "fast_greedy", "edge_betweenness", "label_prop", "leading_eigen","louvain","spinglass","walktrap")
  m.id <- pmatch(method, ms)
  if (is.na(m.id)) {
    stop("Invalid Community Detection Method.")
  }
  method <- ms[m.id]
  
  func <- "cluster_walktrap"
  if (method=="auto") {
    
    functions <- paste("cluster_", ms[-1], sep="")
    cmmms <- rep(-1, length(functions))
    for (i in 1:length(functions)) {
      cm <- NULL
      cm <- do.call(functions[i], args=list(graph=csn, weights=E(csn)$weight))
      cmm <- membership(cm)
      cmmms[i] <- modularity(csn, cmm)
    }
    func <- functions[which.max(cmmms)]
    
  } else {
    func <- paste("cluster_", method, sep="")
  }
  
  cm <- NULL
  cmm <- NULL
  cm <- do.call(func, args=list(graph=csn, weights=E(csn)$weight))
  cmm <- igraph::membership(cm)
  
  return(cmm)
}



getStateGraph <- function(x, dist.method="euclidean", type="mst") {
  
  
  #######################################################################
  ##############           input validation            ##################
  #######################################################################
  
  dms <- c("euclidean", "pearson", "spearman")
  d.id <- pmatch(dist.method, dms)
  if (is.na(d.id)) {
    stop("invalid distance measure.")
  }
  dist.method <- dms[d.id]
  
  #######################################################################
  ##############         cell pairwise distance        ##################
  #######################################################################
  
  # cell pairwise distance
  A <- NULL
  
  if (class(x)=="dist") {
    A <- as.matrix(x)
  } else if (class(x)=="matrix") {
    if (dist.method=="pearson" | dist.method=="spearman") {
      A <- (1-cor(x, method=dist.method))/2
    } else {
      cdist <- dist(t(x), method=dist.method)
      A <- as.matrix(cdist)
    }
  } else {
    stop("invalid format. x shall be of class \"dist\" or \"matrix\".")
  }
  
  #######################################################################
  ##############  constructing cell similarity graph   ##################
  #######################################################################
  
  csn <- NULL 
  
  if (type=="mst") {
    csn <- graph.adjacency(A, mode="undirected", weighted=TRUE)
    csn <- minimum.spanning.tree(csn)
  } else if (type=="graph") {
    csn <- graph.adjacency(A, mode="undirected", weighted=TRUE)
  }
  
  return(csn)
}


plotViolin <- function(x, y, ncol=NULL, nrow=2, scale="free",
                       do.violin=T, 
                       cols=NULL, point.size=1, 
                       x.label="Group", y.label="Expression (+1, log2)", 
                       x.label.size=10, y.label.size=10, 
                       axis.size=8, fig.label.size=12) {
  x.melt <- melt(x, id.vars=c("V"))
  x.melt <- merge(x.melt, y, by.x="variable", by.y="V1")
  
  x.melt$V <- factor(x.melt$V, levels=as.character(x$V))
  
  fmt<-function(x){format(x,nsmall = 2,scientific = FALSE)}
  
  g <- ggplot(x.melt)
  
  if (!(is.null(nrow)) & !(is.null(ncol))) {
    g <- g + facet_wrap("V", nrow=nrow, ncol=ncol, scale="free")
  }
  if (!(is.null(ncol))) {
    g <- g + facet_wrap("V", ncol=ncol, scale="free")
  }
  if (!(is.null(nrow))) {
    g <- g + facet_wrap("V", nrow=nrow, scale="free")
  }
  
  g <- g + geom_point(aes(x=V2, y=value, col=V2), position="jitter", size=point.size)
  g <- g + scale_y_continuous(labels=fmt, limits=c(0, NA))
  
  if (do.violin==T) {
    g <- g + geom_violin(aes(x=V2, y=value, fill=V2), alpha=0.3, scale="width",adjust=0.75)
  }
  
  if (!is.null(cols)) {
    g <- g + scale_fill_manual(values=cols) 
  }
  if(!(is.null(x.label))) {
    g <- g + xlab(x.label)
  }
  if(!(is.null(y.label))) {
    g <- g + ylab(y.label)
  } 
  
  g <- g + theme(strip.background = element_rect(colour = 'white', fill = 'white'), strip.text=element_text(size=fig.label.size, face="italic")) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) + 
    theme(legend.key = element_blank()) +
    theme(axis.text = element_text(size=axis.size)) +
    theme(axis.line.y = element_line(size = 0.2)) +
    theme(axis.line.x= element_line(size = 0.2)) 
  
  print(g)
  
}


ed <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

SLICE_theme_opts <- function() {
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) + 
    theme(legend.key = element_blank())
}











