# combined_benchmark
**Combined Benchmarking Docker and Scripts**

Build directory for yanay/combined_benchmark docker image, hosted publicly on docker hub

This docker container has Seurat 3.0 prerelease and Scanpy 1.3.7 (Latest version) as well as a number of other libraries such as igraph, louvain and umap-learn

**Build Instructions**
*Note: this repo does not contain the h5 files added in this docker image
  # To download the 1.3 million neurons
  cd .../combined_benchmark
  curl -O http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5
  

  bash build_docker
