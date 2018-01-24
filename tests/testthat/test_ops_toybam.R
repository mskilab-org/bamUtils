## tests on mskilab/bamUtils/tests/testthat/toy.bam and index

library(bamUtils)
library(gUtils)  ## remove
##library(testthat)  ## remove

context("test bamUtils on fake BAM, 'small.bam' and index 'small.bam.bai' ")

## 'small.bam' created via chromosome 22 in HCC1143 BL
## $ samtools view -H HCC1143_BL_phased_possorted.bam > header_HCC1143_BL_phased_possorted.sam
## $ samtools view HCC1143_BL_phased_possorted.bam  | grep "chr22" | cat  header_HCC1143_BL_phased_possorted.sam  - | samtools view -Sb - > chr22.unique.bam  ## keeps header
## $ samtools view -h chr22.unique.bam | head -n 100000 > small.sam 
## $ samtools view -S -b small.sam > small.bam
## $ samtools index small.bam

##example_bam = 'small.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
##example_bai = 'small.bam.bai' 

example_bam = 'smallHCC1143BL.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
example_bai = 'smallHCC1143BL.bam.bai' 

small_MD_bam = 'smallHCC1143BL.filtered.MD.bam'
small_MD_bai = 'smallHCC1143BL.filtered.MD.bam.bai'

## example_bam = './tests/testthat/smallHCC1143BL.bam'   ### all tests below are specific to this BAM, and will fail otherwise 
## example_bai = './tests/testthat/smallHCC1143BL.bam.bai'

### MUST CHECK 'pairs.grl' issue with 'get.mate.gr()'
### MUST CHECK stripstrand
### if (stripstrand){
###    strand(intervals) = '*'
### }



test_that('read.bam', {
    
    ## default
    expect_error(read.bam(example_bam))  ## Error in read.bam(example_bam) : Must provide non empty interval list
    ## test 'all' FLAG
    expect_equal(length(read.bam(example_bam, all=TRUE)), 4999)   ## read in entire BAM
    expect_equal(length(read.bam(small_MD_bam, all=TRUE)), 19089)     
    expect_false(is(read.bam(example_bam, all=TRUE), 'GenomicRanges'))   ## check that output is "GRanges" if as.data.table == FALSE (default)
    expect_true(is(read.bam(example_bam, all=TRUE), 'GRangesList')) 
    expect_equal(ncol(read.bam(example_bam, all=TRUE, as.data.table=TRUE)), 16) ## data.table == NULL, should have ncol == 16
    expect_equal(nrow(read.bam(example_bam, all=TRUE, as.data.table=TRUE)), 9998) 
    expect_true(is(read.bam(example_bam, all=TRUE, as.data.table=TRUE), 'data.table'))  
    ## check GRanges for 'intervals' whereby intervals to retrieve overlap BAM reads
    ## check 'chr##' vs. '##' error
    expect_error(read.bam(example_bam, all=TRUE, intervals = GRanges('chr1:10075-10100'))) ## Error in GenomeInfoDb:::getDanglingSeqlevels()
    expect_equal(length(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))), 1027) 
    expect_equal(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]$flag[1], 163)  ## same results with all=TRUE
    ## expect_equal(read.bam(example_bam, all=FALSE, intervals = GRanges('1:10075-10100'))$flag[2], 113)  ## same results with all=FALSE
    expect_match(read.bam(example_bam, intervals = GRanges('1:10075-10100'))[[1]]$cigar[1], '2S6M1I63M1I8M1I8M11S')
    ## check GRanges for 'intervals' whereby intervals to retrieve DO NOT overlap BAM reads
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:500-750'))), 0) 
    ## check 'stripstrand' 10K to 27K
    expect_match(as.character( strand(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10000-27000', strand = "+"), stripstrand = FALSE)[[1]][1])), '+')
    expect_match(as.character( strand(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10000-27000', strand = "+"), stripstrand = FALSE)[[1]][2])), '-')    
    ## issue with 'what'
    ## read.bam(example_bam, all=FALSE, intervals = GRanges('1:10075-10100'), what='qwidth')
    ## read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'), what='MD')  ## error, https://www.rdocumentation.org/packages/Rsamtools/versions/1.24.0/topics/BamInput
    ## verbose
    ## read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'), verbose=TRUE))
    ## check 'tag' works correctly
    expect_true(('R1' %in% colnames(read.bam(example_bam, all=FALSE, intervals = GRanges('1:10075-10100'), tag = 'R1', as.data.table=TRUE))))
    expect_error(('nonsense_tag' %in% colnames(read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'), tag = 'nonsense_tag', as.data.table=TRUE))))
    #### checking Rsamtools::scanBamFlag() flags
    ## 'isPaired'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isPaired = 'foo'))
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isPaired = TRUE)), 1027)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isPaired = FALSE)), 0)
    expect_equal(length(read.bam(example_bam, all = TRUE, isPaired = TRUE)), 4999)
    expect_equal(length(read.bam(example_bam, all = TRUE, isPaired = FALSE)), 0)
    ## 'isProperPair'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isProperPair = 'foo'))   ### Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isProperPair = TRUE)), 896)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isProperPair = FALSE)), 131)
    expect_equal(length(read.bam(example_bam, all = TRUE, isProperPair = TRUE)), 4343)
    expect_equal(length(read.bam(example_bam, all = TRUE, isProperPair = FALSE)), 656)
    ## 'isUnmappedQuery'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isUnmappedQuery = 'foo')) ## Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isUnmappedQuery = TRUE)), 55)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isUnmappedQuery = FALSE)), 975)
    expect_equal(length(read.bam(example_bam, all = TRUE, isUnmappedQuery = TRUE)), 365)
    expect_equal(length(read.bam(example_bam, all = TRUE, isUnmappedQuery = FALSE)), 4724)
    ## 'hasUnmappedMate'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), hasUnmappedMate = 'foo')) ## Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), hasUnmappedMate = TRUE)), 3)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), hasUnmappedMate = FALSE)), 1027)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = TRUE)), 90)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = FALSE)), 4999)
    ## 'isNotPassingQualityControls' 
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isNotPassingQualityControls = 'foo'))  ## Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isNotPassingQualityControls = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isNotPassingQualityControls = FALSE)), 1027)
    expect_equal(length(read.bam(example_bam, all = TRUE, isNotPassingQualityControls = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, all = TRUE, isNotPassingQualityControls = FALSE)), 4999)
    ## 'isDuplicate' 
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isDuplicate = 'foo'))  ## Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isDuplicate = TRUE)), 147)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isDuplicate = FALSE)), 1027)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = TRUE)), 90)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = FALSE)), 4999)
    ## issue with 'ignore.indels'
    ## length(read.bam(example_bam, all = TRUE,  ignore.indels = TRUE))  ### Error in .Call2("explode_cigar_ops", cigar, ops, PACKAGE = "GenomicAlignments") :  'cigar[1]' is NA
    ## Error in duplicated((y[[1]][y[[1]] == M])) : object 'M' not found
    ## > length(read.bam(small_MD_bam , all = TRUE,  ignore.indels = TRUE))
     ## Error in duplicated((y[[1]][y[[1]] == M])) : object 'M' not found

})


### bam.cov.gr


### bam.cov.tile


### get.pairs.grl --- currently non-exported
##  missing code, get.mate.gr() in get.pairs.grl
## Needed for read.bam() as well




test_that('count.clips', {
    ## check errors
    expect_error(count.clips('foo'))
    expect_error(count.clips(example_bam)) ## Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.
    ## default
    expect_equal(count.clips(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]])$right.clips[1], 11)
    expect_equal(count.clips(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]])$right.clips[2], 0)
    expect_equal(count.clips(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]])$left.clips[1], 2)
    expect_equal(count.clips(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]])$left.clips[2], 0)
})


## varbase




test_that('splice.cigar', {
    ## check input
    expect_error(splice.cigar('foo'))
    expect_error(splice.cigar(example_bam)) ##  Error: Reads must be either GRanges, GRangesList, GappedAlignments, or data.table object. Please see documentation for details.
    ## check default
    expect_equal(length(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))), 2054) 
    expect_equal(length(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]), 9)  
    expect_equal(length(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[2]]), 0)  
    expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$type[1], 'S') 
    expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$type[2], 'M') 
    expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$type[3], 'I')
    expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$type[4], 'M')
    expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]])[1], 0)
    expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]])[2], 6)
    expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]])[3], 0)
    expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]])[4], 63)
    expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$fid[1], 1)
    expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$fid[2], 1)
    expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$qname[1], 'C1Y1JACXX130321:7:1101:15151:82244')
    ## fast = FALSE
    expect_equal(length(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast = FALSE)[[2]]), 0)
    expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]])[1], 6)
    expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]])[2], 0)
    expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$type[1], 'M')
    expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$type[2], 'I')
    expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$fid[1], 1)
    expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$fid[2], 1)
    expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$rid[1], 1)
    expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$rid[2], 1)
    ## use.D = FALSE
    expect_true('riid' %in% colnames(as.data.frame(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), use.D = FALSE)[[1]])))  ## should be a 'riid' column in GRanges
    expect_true('riid' %in% colnames(as.data.frame(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), use.D = FALSE)[[2]])))
    ## tests for 'rem.soft' and 'get.seq'
})




test_that('bamflag', {
    ## isPaired  1 0
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isPaired[1], 1)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isPaired[2], 0)
    ## isProperPair  1 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isProperPair[1], 1)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isProperPair[2], 0)
    ## isUnmappedQuery  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isUnmappedQuery[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isUnmappedQuery[2], 0)
    ## hasUnmappedMate  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$hasUnmappedMate[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$hasUnmappedMate[2], 0)
    ## isMinusStrand  1 1 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMinusStrand[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMinusStrand[2], 0)
    ## isMateMinusStrand  1 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMateMinusStrand[1], 1)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMateMinusStrand[2], 0)
    ## isFirstMateRead  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isFirstMateRead[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isFirstMateRead[2], 0)
    ## isSecondMateRead  1 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSecondMateRead[1], 1)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSecondMateRead[2], 0)
    ## isNotPrimaryRead  0 0
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPrimaryRead[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPrimaryRead[2], 0)
    ## isNotPassingQualityControls  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPassingQualityControls[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPassingQualityControls[2], 0)
    ## isDuplicate  0 0 
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isDuplicate[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isDuplicate[2], 0)
    ## isSupplementary  0 0
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSupplementary[1], 0)
    expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSupplementary[2], 0)

})




test_that('bamtag', {
    ## check errors
    expect_error(bamtag('foo'))   ### Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.
    expect_error(bamtag(example_bam))    ### Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.
    ## default
    ### FAILS AT GRANGELIST: WHY? cf read.bam()
    ## expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100')))[1], 'ST-K00126:3:H5TL3BBXX:1:1127:17310:39893_1__') 
    expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]])[1], 'C1Y1JACXX130321:7:1101:15151:82244_2__')   
    ## 'secondary' == TRU#
    expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], secondary = TRUE)[1], 'C1Y1JACXX130321:7:1101:15151:82244_2__0')
    ## 'gr.string' == TRUE
    expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], gr.string = TRUE)[1], 'C1Y1JACXX130321:7:1101:15151:82244_2_1:10034-10117+_')
    expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], gr.string = TRUE)[2], 'C1Y1JACXX130321:7:1101:15151:82244_2_1:10160-10260-_')
    ## 'secondary' == TRUE & 'gr.string' == TRUE
    expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], secondary = TRUE, gr.string = TRUE)[1], 'C1Y1JACXX130321:7:1101:15151:82244_2_1:10034-10117+_0')
    expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], secondary = TRUE, gr.string = TRUE)[2], 'C1Y1JACXX130321:7:1101:15151:82244_2_1:10160-10260-_0')

})




test_that('countCigar', {

    expect_warning(countCigar(example_bam))   ### should be warning message: 'Warning message: In countCigar(example_bam) : NAs introduced by coercion'
    expect_equal(dim(countCigar(example_bam))[1], 1)
    expect_equal(dim(countCigar(example_bam))[2], 4)
    expect_match(colnames(countCigar(example_bam))[1], "D")
    expect_match(colnames(countCigar(example_bam))[2], "I")
    expect_match(colnames(countCigar(example_bam))[3], "M")
    expect_match(colnames(countCigar(example_bam))[4], "S")
    expect_equal(countCigar(example_bam)[1], 0)
    expect_equal(countCigar(example_bam)[2], 0)
    expect_equal(countCigar(example_bam)[3], 0)
    expect_equal(countCigar(example_bam)[4], 0)

})




##test_that('is.paired.end', {
##    expect_equal(as.logical(is.paired.end(example_bam)), TRUE)
##    expect_equal(as.logical(is.paired.end('foo')), NA)   ### error checking, should return NA
##})


