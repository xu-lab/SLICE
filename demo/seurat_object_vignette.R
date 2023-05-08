library(Seurat)
library(SLICE)

# In this vignette we use pbmc data to demonstrate computing cell entropies 
# obtained from a Seurat object. The Seurat object is created from cellranger
# output

###########################################################################
#####                                                               #######
#####               Loading data and preprocessing                  #######
#####                                                               #######
###########################################################################

#change the following path string to the folder address where pbmc cellranger
# output count matrix is stored
counts_folder = "."

# Load the PBMC dataset
# and follow the pre-process steps from Seurat tutorial
pbmc <- CreateSeuratObject(counts = Read10X(data.dir = counts_folder, gene.column = 2),
                           project = "pbmc3k", min.cells = 3, 
                           min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)

# find clusters following the code in Seurat tutorial
# to use as cellidentity to create slice object
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc$celltypes = factor(unlist(lapply(pbmc$seurat_clusters, function(x){new.cluster.ids[[x]]})))

#run UMAP and DimPlot to check everything worked as expected
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "celltypes") + NoLegend()

# if scale=TRUE the scaled highly variable genes are used, if not the raw counts are used
sc <- construct_from_seurat(obj=pbmc,                # the Seurat object to get counts from
                            cellidentity="celltypes",       # cellidentities which can be a meta.data column or character vector
                            scale=TRUE,              # if the scaled counts or raw counts should be used
                            projname=NULL)           # project name to be passed to slice object

###########################################################################
#####                                                               #######
#####     Measuring cell differentiation state using scEntropy      #######
#####                                                               #######
###########################################################################


# loading the pre-computed human gene-gene Kappa similarity matrix
# The matrix will be loaded into a variable named "km"
data(hs_kappasim)

# bootstrap calculation of scEntropy
sc <- getEntropy(sc, km=kh,                             # use the pre-computed kappa similarity matrix of human genes
                 calculation="bootstrap",               # choose the bootstrap calculation
                 B.num=100,                             # 100 iterations
                 exp.cutoff=1,                          # the threshold for expressed genes
                 B.size=1000,                           # the size of bootstrap sample
                 clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters  
                 random.seed=201602)                    # set the random seed to reproduce the results in the paper

plotEntropies(sc)
