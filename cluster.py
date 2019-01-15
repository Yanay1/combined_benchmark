#! /usr/bin/env python
import sys
import numpy as np
import pandas as pd
import scanpy.api as sc
import datetime



verbose = 3

if (len(sys.argv) == 1):
  print("Running on 25K subsample")
  filename = "home/1M_neurons_matrix_subsampled_25K.h5"
elif len(sys.argv) > 1:
  # default output file
  filename = sys.argv[1]
  if len(sys.argv) == 3:
    verbose = int(sys.argv[2])


start = datetime.datetime.now()

sc.settings.verbosity = verbose             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

print("Reading, making names unique and filtering:")
now = datetime.datetime.now()
adata = sc.read_10x_h5(filename)

adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print("Reading, making names unique and filtering Time:", datetime.datetime.now() - now)

print("Log Normalizing Data")
now = datetime.datetime.now()
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e3)
sc.pp.log1p(adata)
print("Log Normalizing Data Time:", datetime.datetime.now() - now)

adata.raw = adata
print("Finding Variable Features")
now = datetime.datetime.now()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Finding Variable Features Time:", datetime.datetime.now() - now)
adata = adata[:, adata.var['highly_variable']]

print("Scaling Data")
now = datetime.datetime.now()
sc.pp.scale(adata, max_value=10)
print("Scaling Data Time:", datetime.datetime.now() - now)

print("Running PCA")
now = datetime.datetime.now()
sc.tl.pca(adata, svd_solver='arpack')
print("Running PCA Time:", datetime.datetime.now() - now)

print("Finding Neighbors")
now = datetime.datetime.now()
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
print("Finding Neighbors Time:", datetime.datetime.now() - now)

print("Finding Clusters (Louvain)")
now = datetime.datetime.now()
sc.tl.louvain(adata)
print("Finding Clusters (Louvain) Time:", datetime.datetime.now() - now)

print("Running UMAP")
now = datetime.datetime.now()
sc.tl.umap(adata)
print("Running UMAP Time:", datetime.datetime.now() - now)

print("Final Object:")
print(adata)

print("Total Benchmark Time:", datetime.datetime.now() - start)
