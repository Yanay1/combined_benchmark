#! /usr/bin/env Rscript
# SEURAT VERSION 3.0 PRE RELEASE
library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
verbose = TRUE
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  print("Running on 25K subsample")
  filename <- "home/1M_neurons_matrix_subsampled_25K.h5"
} else if (length(args) >= 1) {
  # default output file
  filename <- args[1]
  if (length(args) == 2) {
    verbose <- as.logical(args[2])
  }
}

start <- Sys.time()

# Load the Neuron dataset
print("Reading, making names unique and filtering:")
now <- Sys.time()
neuron.data <- Read10X_h5(filename)
rownames(neuron.data) <- make.names(neuron.data[,1], unique = TRUE)

neuron <- CreateSeuratObject(counts = neuron.data, min.cells = 3, min.features = 200, project = "10X_neuron")
print("Reading, making names unique and filtering:")
print(Sys.time() - now)

print("Log Normalizing Data")
now <- Sys.time()
neuron <- NormalizeData(object = neuron, normalization.method = "LogNormalize", scale.factor = 1e3, verbose = verbose)
print("Log Normalizing Data Time:")
print(Sys.time() - now)

print("Finding Variable Features")
now <- Sys.time()
neuron <- FindVariableFeatures(object = neuron, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), verbose = verbose)
print("Finding Variable Features Time:")
print(Sys.time() - now)

print("Scaling Data")
now <- Sys.time()
neuron <- ScaleData(object = neuron, features = rownames(x = neuron), verbose = verbose)
print("Scaling Data Time:")
print(Sys.time() - now)

print("Running PCA")
now <- Sys.time()
neuron <- RunPCA(object = neuron, features = VariableFeatures(object = neuron), verbose = verbose)
print("Runing PCA Time:")
print(Sys.time() - now)

print("Finding Neighbors")
now <- Sys.time()
neuron <- FindNeighbors(object = neuron, dims = 1:10, verbose = verbose)
print("Finding Neighbors Time:")
print(Sys.time() - now)

print("Finding Clusters (Louvain)")
now <- Sys.time()
neuron <- FindClusters(object = neuron, resolution = 0.6, verbose = verbose)
print("Finding Clusters (Louvain) Time:")
print(Sys.time() - now)

print("Running UMAP")
now <- Sys.time()
neuron <- RunUMAP(neuron, reduction.use = "pca", dims = 1:10, verbose = verbose)
print("Running UMAP Time:")
print(Sys.time() - now)

print("Final Object")
neuron

print("Total Benchmark Time:")
print(Sys.time() - start)
