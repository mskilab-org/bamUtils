#!/bin/bash

Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); biocLite("Rsamtools"); biocLite("openssl");
biocLite("rtracklayer");biocLite("GenomicAlignments");library(openssl);library(rtracklayer);library(GenomicAlignments);
biocLite("BSgenome.Hsapiens.UCSC.hg19");biocLite("GenomicRanges"); install.packages("devtools"); devtools::install_github("mskilab/gUtils"); 
devtools::install_github("jimhester/covr");devtools::install_github("mskilab/skitools");library(skitools);'
