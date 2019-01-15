#! /usr/bin/env R
Sys.setenv(GITHUB_PAT="0a99d5f8d571514bab5c0ef7caaa6231f3e16ab2")
install.packages('devtools',repos = "http://cran.us.r-project.org")
devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
install.packages("hdf5r")