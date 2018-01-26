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


## create tumor/normal pairs with FASTA and VCF

### FASTA

samtools faidx human_g1k_v37_decoy.fasta

samtools faidx human_g1k_v37_decoy.fasta chr1 > your_subset_file.fa

cat chr1_human_g1k_v37_decoy.fasta | head -n 200000 > chr1_human_g1k_v37_decoy.subset.fasta

samtools faidx  chr1_human_g1k_v37_decoy.subset.fasta


## TUMOR

### we have a small version of HCC1143BL
## now get HCC1143

samtools view -H HCC1143.final.bam  > HCC1143.header.sam
samtools view HCC1143.final.bam  | head -n 10000 | cat  HCC1143.header.sam - | samtools view -Sb - > HCC1143.final.subset.bam 

samtools index HCC1143.final.subset.bam 

### VCF

## module load vcftools/0.1.14
## vcf-validator HCC1143_nygc.snowman.somatic.sv.vcf   ## buggy

## use this
grep -E '^(#|1[[:space:]])' HCC1143_nygc.snowman.somatic.sv.vcf  > chrom1.vcf 

## vcf-validator chrom1.vcf    ## obviously still buggy

##module load  htslib/1.5
##module load  bcftools/1.4.1 

## htslib commands
##bgzip -c HCC1143_nygc.snowman.somatic.sv.vcf  > HCC1143_nygc.snowman.somatic.sv.vcf.gz
##tabix -p vcf HCC1143_nygc.snowman.somatic.sv.vcf.gz

##bcftools filter  HCC1143_nygc.snowman.somatic.sv.vcf.gz -r 1

## [E::vcf_parse_format] Invalid character '.' in 'PL' FORMAT field at 1:6525379

## /gpfs/commons/groups/imielinski_lab/projects/CellLines/Flow/Snowman/HCC1143_nygc/HCC1143_nygc.snowman.somatic.sv.vcf

## R
## > library(skidb)
## > library(data.table)
## > totalvcf = read_vcf('')
## > chr1 = totalvcf[seqnames(totalvcf) == 1]
## > length(chr1) ## 60 
## 








