#!/bin/bash

Rscript -e 'source("https://bioconductor.org/biocLite.R"); 
            biocLite("BiocInstaller"); 
            biocLite("Rsamtools");
            library(Rsamtools);
  
            biocLite("rtracklayer");
            library(rtracklayer);
            biocLite("GenomicAlignments");
            library(GenomicAlignments);
            biocLite("BSgenome.Hsapiens.UCSC.hg19");
            biocLite("GenomicRanges"); 
            install.packages("devtools"); 
            devtools::install_github("mskilab/gUtils"); 
            devtools::install_github("jimhester/covr");'
