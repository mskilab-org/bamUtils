#!/usr/bin/env bash

## create example BAMs to be used in tests using HCC1143BL

## /gpfs/commons/groups/imielinski_lab/data/CellLines/HCC1143/WGS/HCC1143BL.final.bam

source ~/.bash_profile


module unload samtools/1.1
module load samtools/1.6


samtools view -H HCC1143BL.final.bam > header.sam
samtools view HCC1143BL.final.bam | head -n 10000 | cat header.sam - | samtools view -Sb - > smallHCC1143BL.bam

samtools index smallHCC1143BL.bam


## samtools view -f 0x02 -F 0x10 ./tests/testthat/small.bam -q 30 

samtools view -f 0x02 -F 0x10 HCC1143BL.final.bam -q 30 | cat header.sam - | samtools view -Sb - > HCC1143BL.filtered.bam


samtools view HCC1143BL.filtered.bam | head -n 20000 | cat header.sam - | samtools view -Sb - > smallHCC1143BL.filtered.bam


samtools index smallHCC1143BL.filtered.bam
### add MD tags to BAM and retrieve reads via 

## use GATK FASTA for hg19
## /gpfs/commons/home/biederstedte-934/DB/GATK/human_g1k_v37_decoy.fasta



## samtools calmd 
samtools calmd smallHCC1143BL.filtered.bam /gpfs/commons/home/biederstedte-934/DB/GATK/human_g1k_v37_decoy.fasta > smallHCC1143BL.filtered.MD.sam

samtools view -Sb  smallHCC1143BL.filtered.MD.sam >  smallHCC1143BL.filtered.MD.bam

samtools index smallHCC1143BL.filtered.MD.bam





