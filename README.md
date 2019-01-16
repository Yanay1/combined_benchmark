# Combined Benchmarking Docker and Scripts

*Build directory for yanay/combined_benchmark docker image, hosted publicly on docker hub*

This docker container has Seurat 3.0 prerelease and Scanpy 1.3.7 (Latest version) as well as a number of other libraries such as igraph, louvain and umap-learn

## Build Instructions
*Note: this repo does not contain the h5 files added in this docker image*
*To build this image you need to generate them, or just remove the lines where they are added to the dockerfile*

  cd .../combined_benchmark
  chmod a+x subsample.py # make this file executable locally
  # To download the 1.3 million neurons and subsample (this can be done locally)
  curl -O http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5
  subsample.py input_file.h5 output_file.h5 cells_to_subsample_to (optional genome arg)
  # example to generate the 1M sample
  subsample.py 1M_neurons_filtered_gene_bc_matrices_h5.h5 1M_neurons_matrix_subsampled_1M.h5 1000000 mm10

**If you don't want to add these quite large files to your version**
  # remove these lines from the docker file
  ADD 1M_neurons_matrix_subsampled_1M.h5 /home/1M_neurons_matrix_subsampled_1M.h5
  ADD 1M_neurons_matrix_subsampled_500K.h5 /home/1M_neurons_matrix_subsampled_500K.h5
  ADD 1M_neurons_matrix_subsampled_250K.h5 /home/1M_neurons_matrix_subsampled_250K.h5
  ADD 1M_neurons_matrix_subsampled_100K.h5 /home/1M_neurons_matrix_subsampled_100K.h5
  ADD 1M_neurons_matrix_subsampled_25K.h5 /home/1M_neurons_matrix_subsampled_25K.h5
  ADD 1M_neurons_filtered_gene_bc_matrices_h5.h5 /home/1M_neurons_filtered_gene_bc_matrices_h5.h5

### Build Steps
  cd .../combined_benchmark
  # first set version in VERSION.txt and change tag name in build_docker and push_docker
  bash build_docker

### Push Steps
  cd .../combined_benchmark
  # first set version in VERSION.txt and change tag name in build_docker and push_docker
  # after building
  bash push_docker

## Usage
To use the image in google cloud VM make sure to add yanay/combined_benchmark as the container when creating the VM.

First, start the docker container
  # In the ssh window of the VM or locally
  docker run -it yanay/combined_benchmark bash

To benchmark Seurat
  seurat_benchmark.R input_h5 optional_verbose_setting
  # verbose setting true or false (display loading bar during steps)

To Benchmark Scanpy
  scanpy_benchmark.py input_h5 optional_verbose_setting
  # verbose options 0, 1, 2, 3: errors (0), warnings (1), info (2), hints (3)

To create more subsamples
  subsample.py input_file.h5 output_file.h5 cells_to_subsample_to (optional genome arg)

## The Benchmark
| Step        | Seurat Command           | Scanpy Command  |
| ------------- |-------------| -----|
| Read the h5 file, make names unique and filter matrix to min 3 cells per gene and min 200 genes per cell      | Read10X_h5, make.names, CreateSeuratObject, | sc.read_10x_h5, var_names_make_unique, filter_cells, filter_genes |
| Log Normalize the Data to 1e3 counts| NormalizeData |   sc.pp.normalize_per_cell, sc.pp.log1p |
| Find variable features | FindVariableFeatures      |    highly_variable_genes |
| Scale data | ScaleData |    sc.pp.scale |
| Run PCA | RunPCA      |    sc.tl.pca |
| Find neighbors | FindNeighbors      |    sc.pp.neighbors |
| Find Clusters (Louvain) | FindClusters  |    sc.tl.louvain |
| Run UMAP | RunUMAP  | sc.tl.umap |
| View the Final object | neuron  | adata |

At the start of each step, the step name will be printed.

At the end of each step, the time it took will be printed

At the end of the benchmark a summary of the final object and the total benchmark time will be printed.
