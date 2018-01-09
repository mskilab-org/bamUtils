## tests on mskilab/bamUtils/tests/testthat/toy.bam and index

library(bamUtils)
library(testthat)  ## remove

context("test bamUtils on fake BAM, 'small.bam' and index 'small.bam.bai' ")

## 'small.bam' created via chromosome 22 in HCC1143 BL
## $ samtools view -H HCC1143_BL_phased_possorted.bam > header_HCC1143_BL_phased_possorted.sam
## $ samtools view HCC1143_BL_phased_possorted.bam  | grep "chr22" | cat  header_HCC1143_BL_phased_possorted.sam  - | samtools view -Sb - > chr22.unique.bam  ## keeps header
## $ samtools view -h chr22.unique.bam | head -n 100000 > small.sam 
## $ samtools view -S -b small.sam > small.bam
## $ samtools index small.bam

example_bam = 'small.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
example_bai = 'small.bam.bai' 

## example_bam = './tests/testthat/small.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
## example_bai = './tests/testthat/small.bam.bai'

### MUST CHECK 'pairs.grl' issue with 'get.mate.gr()'
### MUST CHECK stripstrand
### if (stripstrand){
###    strand(intervals) = '*'
### }



test_that('read.bam', {
	## test 'all' FLAG
	expect_error(read.bam(example_bam))  ## custom error message, all=FALSE by default
    expect_equal(length(read.bam(example_bam, all=TRUE)), 97383)   ## read in entire BAM
    expect_true(is(read.bam(example_bam, all=TRUE), 'GenomicRanges'))   ## check that output is "GRanges" if as.data.table == FALSE (default)
    expect_equal(ncol(read.bam(example_bam, all=TRUE, as.data.table=TRUE)), 16) ## data.table == NULL, should have ncol == 16
    expect_equal(nrow(read.bam(example_bam, all=TRUE, as.data.table=TRUE)), 97383) ## data.table == NULL, should have ncol == 16
    expect_true(is(read.bam(example_bam, all=TRUE, as.data.table=TRUE), 'data.table'))  
    ## check GRanges for 'intervals' whereby intervals to retrieve overlap BAM reads
    expect_equal(length(read.bam(example_bam, all=TRUE, intervals = GRanges('chr1:10075-10100'))), 2) 
    expect_equal(read.bam(example_bam, all=TRUE, intervals = GRanges('chr1:10075-10100'))$flag[2], 113)  ## same results with all=TRUE
    expect_equal(read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'))$flag[2], 113)  ## same results with all=FALSE
    expect_match(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'))$cigar[2], '55M4D69M3S')
    ## check GRanges for 'intervals' whereby intervals to retrieve DO NOT overlap BAM reads
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:500-750'))), 0) 
    ## check 'stripstrand' 10K to 27K
    ## issues here: 
    ##expect_match(as.character(strand(read.bam(example_bam, all=TRUE, intervals = GRanges('chr1:10000-27000', strand = "+"), stripstrand = FALSE)[1])), '-')
    ##expect_match(as.character(strand(read.bam(example_bam, all=TRUE, stripstrand = FALSE)[6])), '+')
    ## issue with 'what'
    ## read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'), what='qwidth')
    ## read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'), what='MD')  ## error, https://www.rdocumentation.org/packages/Rsamtools/versions/1.24.0/topics/BamInput
    ## verbose
    ## read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'), verbose=TRUE))
    ## check 'tag' works correctly
    expect_true(('R1' %in% colnames(read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'), tag = 'R1', as.data.table=TRUE))))
    expect_error(('nonsense_tag' %in% colnames(read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'), tag = 'nonsense_tag', as.data.table=TRUE))))
    ## checking Rsamtools::scanBamFlag() flags
    ## 'isPaired'
    expect_error(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isPaired = 'foo'))
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isPaired = TRUE)), 2)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isPaired = FALSE)), 0)
    expect_equal(length(read.bam(example_bam, all = TRUE, isPaired = TRUE)), 97383)
    expect_equal(length(read.bam(example_bam, all = TRUE, isPaired = FALSE)), 0)
    ## 'isProperPair'
    expect_error(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isProperPair = 'foo'))
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isProperPair = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isProperPair = FALSE)), 2)
    expect_equal(length(read.bam(example_bam, all = TRUE, isProperPair = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, all = TRUE, isProperPair = FALSE)), 97383)
    ## 'isUnmappedQuery'
    expect_error(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isUnmappedQuery = 'foo'))
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isUnmappedQuery = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isUnmappedQuery = FALSE)), 2)
    expect_equal(length(read.bam(example_bam, all = TRUE, isUnmappedQuery = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, all = TRUE, isUnmappedQuery = FALSE)), 97383)
    ## 'hasUnmappedMate'
    expect_error(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), hasUnmappedMate = 'foo'))
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), hasUnmappedMate = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), hasUnmappedMate = FALSE)), 2)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = TRUE)), 18439)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = FALSE)), 78944)
    ## 'isNotPassingQualityControls' 
    expect_error(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isNotPassingQualityControls = 'foo'))
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isNotPassingQualityControls = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isNotPassingQualityControls = FALSE)), 2)
    expect_equal(length(read.bam(example_bam, all = TRUE, isNotPassingQualityControls = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, all = TRUE, isNotPassingQualityControls = FALSE)), 97383)
    ## 'isDuplicate' 
    expect_error(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isDuplicate = 'foo'))
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isDuplicate = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('chr1:10075-10100'), isDuplicate = FALSE)), 2)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = TRUE)), 18439)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = FALSE)), 78944)
    ## issue with 'ignore.indels'
    ## length(read.bam(example_bam, all = TRUE,  ignore.indels = TRUE))
    ## Error in duplicated((y[[1]][y[[1]] == M])) : object 'M' not found

})






##test_that('bamflag', 
##
##
##	)

