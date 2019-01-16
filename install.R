#! /usr/bin/env R
install.packages('devtools',repos = "http://cran.us.r-project.org")
devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
install.packages("hdf5r")
