#!/usr/bin/env bash

## create example BAMs to be used in tests using HCC1143BL

## /gpfs/commons/groups/imielinski_lab/data/CellLines/HCC1143/WGS/HCC1143BL.final.bam

source ~/.bash_profile


module unload samtools/1.1
module load samtools/1.6


samtools view -H HCC1143BL.final.bam > header.sam
samtools view HCC1143BL.final.bam | head -n 10000 | cat header.sam - | samtools view -Sb - > smallHCC1143BL.bam

samtools index smallHCC1143BL.bam

### add MD tags to BAM and retrieve reads via 

## use GATK FASTA for hg19
## /gpfs/commons/home/biederstedte-934/DB/GATK/human_g1k_v37_decoy.fasta

## samtools calmd 
samtools calmd