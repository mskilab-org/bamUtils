#!/bin/bash

Rscript -e 'install.packages("devtools"); library(devtools); install_github("mskilab/gUtils"); install_github("jimhester/covr"); source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); biocLite("GenomicRanges")'
