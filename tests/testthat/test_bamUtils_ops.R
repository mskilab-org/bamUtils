
library(bamUtils)
library(testthat)
library(gUtils)
context("test bamUtils on fake BAM, 'small.bam' and index 'small.bam.bai' ")
example_bam = system.file("extdata", 'smallHCC1143BL.bam', package = "bamUtils")   ### all tests below are specific to this BAM, and will fail otherwise
example_bai = system.file("extdata", 'smallHCC1143BL.bam.bai', package = "bamUtils")
small_MD_bam = system.file("extdata", 'smallHCC1143BL.filtered.MD.bam', package = "bamUtils")
small_MD_bai = system.file("extdata", 'smallHCC1143BL.filtered.MD.bam.bai', package = "bamUtils")
tumor_bam = system.file("extdata", 'HCC1143.final.subset.bam', package = "bamUtils")
tumor_bai = system.file("extdata", 'HCC1143.final.subset.bam.bai', package = "bamUtils")
small_reference = system.file("extdata", 'chr1_human_g1k_v37_decoy.subset.fasta',package = "bamUtils")
somatic_vcf = system.file("extdata", 'chrom1.vcf', package = "bamUtils")
normalbam = system.file("extdata", 'smallA4AD_BL.bam', package = "bamUtils")
tumorbam = system.file("extdata", 'smallA4AD_tum.bam', package = "bamUtils")
mafpath = system.file("extdata", 'snv.annotated.A4AD.maf', package = "bamUtils")
noindexbam = system.file("extdata", 'bam_noindex.bam', package = "bamUtils")

message("Tests starting...")

test_that('read.bam', {
    ## default
    ##   Cannot open BAM. A valid BAM for 'bam' must be provided.    
    expect_error(read.bam('fake_bam.txt'))
    expect_error(read.bam(example_bam))  ## Error in read.bam(example_bam) : Must provide non empty interval list
    expect_error(read.bam(noindexbam, all=TRUE)) ## Error in value[[3L]](cond) : valid 'index' file required

    ## test 'all' FLAG
    tmp = read.bam(example_bam, all=TRUE)
    expect_equal(length(tmp), 4999)   ## read in entire BAM
    expect_true(inherits(read.bam(example_bam, all=TRUE), 'GRangesList'))

    tmp = read.bam(example_bam, all=TRUE, as.data.table = TRUE)
    expect_equal(ncol(tmp), 16) ## data.table == NULL, should have ncol == 16
    expect_equal(nrow(tmp), 9998)
    expect_true(inherits(tmp, 'data.table'))
    
    tmp = read.bam(small_MD_bam, all=TRUE)
    expect_equal(length(tmp), 19089)

    ## ignroe.indels
    tmp = read.bam(example_bam, all=TRUE, ignore.indels=TRUE)
    expect_equal(length(tmp), 4999)
    
    ## check GRanges for 'intervals' whereby intervals to retrieve overlap BAM reads
    ## check 'chr##' vs. '##' error
    expect_error(read.bam(example_bam, all=TRUE, intervals = GRanges('chr1:10075-10100'))) ## ERROR in GenomeInfoDb:::getDanglingSeqlevels()
    
    tmp = read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))
    expect_equal(length(tmp), 987)
    expect_equal(tmp[[1]]$flag[1], 163)  ## same results with all=TRUE
    expect_match(tmp[[1]]$cigar[1], '2S6M1I63M1I8M1I8M11S')
    
    ## check GRanges for 'intervals' whereby intervals to retrieve DO NOT overlap BAM reads
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:500-750'))), 0)
    
    ## check 'stripstrand' 10K to 27K
    tmp = read.bam(example_bam, all=TRUE, intervals = GRanges('1:10000-27000', strand = "+"), stripstrand = FALSE)
    expect_match(as.character( strand(tmp[[1]][1])), '+')
    expect_match(as.character( strand(tmp[[1]][2])), '-')
    
    ## verbose
    expect_equal(length(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'), verbose=TRUE)), 987)
    

    expect_equal(length(read.bam(example_bam, all=TRUE, intervals = grl2, verbose=TRUE)), 0)
    expect_equal(length(read.bam(example_bam, all=TRUE, intervals = as.data.frame(gr2dt(GRanges('1:10075-10100'))), verbose=TRUE)), 987)
    
    ## check 'tag' works correctly
    expect_true(('R1' %in% colnames(read.bam(example_bam, all=FALSE, intervals = GRanges('1:10075-10100'), tag = 'R1', as.data.table=TRUE, verbose=TRUE))))
    expect_error(('nonsense_tag' %in% colnames(read.bam(example_bam, all=FALSE, intervals = GRanges('chr1:10075-10100'), tag = 'nonsense_tag', as.data.table=TRUE))))
    #### checking Rsamtools::scanBamFlag() flags
    ## 'isPaired'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isPaired = 'foo'))
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isPaired = TRUE)), 987)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isPaired = FALSE)), 0)
    expect_equal(length(read.bam(example_bam, all = TRUE, isPaired = TRUE)), 4999)
    expect_equal(length(read.bam(example_bam, all = TRUE, isPaired = FALSE)), 0)

    
    ## 'isProperPair'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isProperPair = 'foo'))   ### Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isProperPair = TRUE)), 896)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isProperPair = FALSE)), 91)
    expect_equal(length(read.bam(example_bam, all = TRUE, isProperPair = TRUE)), 4343)
    expect_equal(length(read.bam(example_bam, all = TRUE, isProperPair = FALSE)), 656)
    ## 'isUnmappedQuery'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isUnmappedQuery = 'foo')) ## Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isUnmappedQuery = TRUE)), 14)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isUnmappedQuery = FALSE)), 975)
    expect_equal(length(read.bam(example_bam, all = TRUE, isUnmappedQuery = TRUE)), 365)
    expect_equal(length(read.bam(example_bam, all = TRUE, isUnmappedQuery = FALSE)), 4724)
    ## 'hasUnmappedMate'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), hasUnmappedMate = 'foo')) ## Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), hasUnmappedMate = TRUE)), 3)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), hasUnmappedMate = FALSE)), 986)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = TRUE)), 90)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = FALSE)), 4999)
    ## 'isNotPassingQualityControls'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isNotPassingQualityControls = 'foo'))  ## Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isNotPassingQualityControls = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isNotPassingQualityControls = FALSE)), 987)
    expect_equal(length(read.bam(example_bam, all = TRUE, isNotPassingQualityControls = TRUE)), 0)
    expect_equal(length(read.bam(example_bam, all = TRUE, isNotPassingQualityControls = FALSE)), 4999)
    ## 'isDuplicate'
    expect_error(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isDuplicate = 'foo'))  ## Error in !is.na(x) && x : invalid 'y' type in 'x && y'
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isDuplicate = TRUE)), 147)
    expect_equal(length(read.bam(example_bam, intervals = GRanges('1:10075-10100'), isDuplicate = FALSE)), 987)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = TRUE)), 90)
    expect_equal(length(read.bam(example_bam, all = TRUE, hasUnmappedMate = FALSE)), 4999)
    ## issue with 'ignore.indels'
    ## length(read.bam(example_bam, all = TRUE,  ignore.indels = TRUE))  ### Error in .Call2("explode_cigar_ops", cigar, ops, PACKAGE = "GenomicAlignments") :  'cigar[1]' is NA
    ## Error in duplicated((y[[1]][y[[1]] == M])) : object 'M' not found
    ## > length(read.bam(small_MD_bam , all = TRUE,  ignore.indels = TRUE))
    ## Error in duplicated((y[[1]][y[[1]] == M])) : object 'M' not found

})


### bam.cov.gr
test_that('bam.cov.gr', {

    ## Cannot open BAM. A valid BAM for 'bam' must be provided.
    expect_error(bam.cov.gr('fake_bam.txt'))
    ## if (missing(bam) | missing(intervals)){
    expect_error(bam.cov.gr())
    expect_error(bam.cov.gr(example_bam, intervals=NULL))  ##  Error: Granges of intervals to retrieve 'intervals' must be in the format 'GRanges'. Please see documentation for details.
    ## if (!is.null(bai)){
    expect_error(bam.cov.gr(noindexbam, intervals = GRanges('1:10075-10100')))
    ## intervals
    expect_equal(width(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'))), 26)
    expect_match(as.character(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'))$file), 'smallHCC1143BL.bam')
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'))$records), 1065)
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'))$nucleotides), 107565)
    ## all
    ## verbose
    expect_equal(width(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), verbose=TRUE)), 26)
    ## count.all
    expect_equal(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), count.all=TRUE)$records, 1331)
    expect_equal(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), count.all=TRUE)$nucleotides, 134431)
    ## isPaired
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isPaired = FALSE)$records), 0)
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isPaired = FALSE)$nucleotides), 0)
    ## isProperPair
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isProperPair = FALSE)$records), 77)
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isProperPair = FALSE)$nucleotides), 7777)
    ## isUnmappedQuery
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isUnmappedQuery = TRUE)$records), 0)
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isUnmappedQuery = TRUE)$nucleotides), 0)
    ## hasUnmappedMate
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), hasUnmappedMate = TRUE)$records), 0)
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), hasUnmappedMate = TRUE)$nucleotides), 0)
    ## isNotPassingQualityControls turned. See documentation for Rsamtools::scanBamFlag(). (default == NA)
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isNotPassingQualityControls = TRUE)$records), 0)
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isNotPassingQualityControls= TRUE)$nucleotides), 0)
    ## isDuplicate
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isDuplicate = TRUE)$records), 98)
    expect_equal(as.integer(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), isDuplicate = TRUE)$nucleotides), 9898)
    ## mc.cores
    expect_equal(width(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), mc.cores = 2)), 26)
    ## chunksize
    expect_equal(width(bam.cov.gr(example_bam, intervals = GRanges('1:10075-10100'), chunksize = 1)), 26)
    ##

})



test_that('bam.cov.tile', {

bam.cov.tile(example_bam, window=1e7)
##     ## Cannot open BAM. A valid BAM for 'bam.file' must be provided.
##     expect_error(bam.cov.tile('fake_bam.txt'))
##     ## expect_equal(length(bam.cov.tile(example_bam)), 31018086)  ## Travis complains: cannot allocate vector of size 3.4 Gb
##     ## window
##     expect_equal(length(bam.cov.tile(example_bam, window=1e7)), 382)
##     ## chunksize
##     expect_equal(length(bam.cov.tile(example_bam, chunksize=5e4, window=1e6)), 3173)
##     ## min.map
##     expect_equal(bam.cov.tile(example_bam, window=1e7, verbose=FALSE, min.map=60)[1]$count, 14)
##     expect_equal(bam.cov.tile(example_bam, window=1e7, verbose=FALSE, min.map=15)[1]$count, 932)
##     ## verbose
##     expect_equal(length(bam.cov.tile(example_bam, window=1e7, verbose=FALSE)), 382)
##     ## max.tlen
##     expect_equal( max(bam.cov.tile(example_bam, window=1e7, verbose=FALSE, max.tlen = 1e7)$count), 129)
##     expect_equal(max(bam.cov.tile(example_bam, window=1e7, verbose=FALSE, max.tlen = 100)$count), 0)
##     ## st.flag
##     expect_equal(max(bam.cov.tile(example_bam, window=1e7, verbose=FALSE, max.tlen = 1e6)$count), 129)
##     expect_equal(max(bam.cov.tile(example_bam, window=1e7, verbose=FALSE, max.tlen = 1e6, st.flag = '')$count), 359)
##     ## fragments
##     expect_equal(max(bam.cov.tile(example_bam, window=1e7, verbose=FALSE, max.tlen = 1e6, st.flag = '', fragments=FALSE)$count), 360)
##     ## regions
##     ## do.gc
## #    expect_equal(max(bam.cov.tile(example_bam, window=1e7, verbose=FALSE, max.tlen = 1e6, st.flag = '', do.gc=TRUE)$count), 359)
##     ## midpoint
##     expect_equal(length(bam.cov.tile(example_bam, window=1e7, verbose=FALSE, min.map=60, midpoint = FALSE)), 382)
##     ## bam.cov.tile
## #    expect_equal(bam.cov.tile(example_bam, window=1e7, verbose=TRUE, min.map=60)[1]$count, 14)


})







##
## get.mate.gr()
test_that('get.mate.gr', {

    ## expect_equal(width(get.mate.gr(read.bam(example_bam, all=TRUE)[[1]])[1]), 101)
    ## expect_equal(width(get.mate.gr(read.bam(example_bam, all=TRUE)[[1]])[2]), 101)
    ## expect_match(get.mate.gr(read.bam(example_bam, all=TRUE)[[1]])[1]$qname, 'C1Y1JACXX130321:7:1101:13020:82300')
    ## expect_equal(get.mate.gr(read.bam(example_bam, all=TRUE)[[1]])[1]$mapq, 6)
    ## expect_match(get.mate.gr(read.bam(example_bam, all=TRUE)[[1]])[2]$qname, 'C1Y1JACXX130321:7:1101:13020:82300')
    ## expect_equal(get.mate.gr(read.bam(example_bam, all=TRUE)[[1]])[2]$mapq, 29)

})



## get.pairs.grl
test_that('get.pairs.grl', {

    ## expect_error(get.pairs.grl(read.bam(example_bam, all=TRUE))) ## Error in (function (classes, fdef, mtable)  :  unable to find an inherited method for function ‘granges’ for signature ‘"GRangesList"’
    ## expect_true(inherits(get.pairs.grl(read.bam(example_bam, all=TRUE)[[1]], verbose=TRUE), 'GRangesList'), TRUE)
    ## ## pairs.grl.split
    ## expect_true(inherits(get.pairs.grl(read.bam(example_bam, all=TRUE)[[1]], pairs.grl.split=FALSE, verbose=TRUE), 'GRanges'), TRUE)
    ## ## verbose
    ## expect_true(inherits(get.pairs.grl(read.bam(example_bam, all=TRUE)[[1]], verbose=TRUE), 'GRangesList'))
    ## expect_true(inherits(get.pairs.grl(read.bam(example_bam, all=TRUE)[[1]], verbose=TRUE, pairs.grl.split=FALSE), 'GRanges'), TRUE)
    ## ## if (inherits(reads, 'GappedAlignmentPairs')){
    ## gapp = readGAlignmentPairs(example_bam)


})




test_that('count.clips', {
    ## ## check errors
    ## expect_error(count.clips('foo'))
    ## expect_error(count.clips(example_bam)) ## Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.
    ## ## default
    ## expect_equal(count.clips(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]])$right.clips[1], 11)
    ## expect_equal(count.clips(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]])$right.clips[2], 0)
    ## expect_equal(count.clips(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]])$left.clips[1], 2)
    ## expect_equal(count.clips(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]])$left.clips[2], 0)
    ## ## check 'if (length(reads) == 0)'
    ## expect_equal(length(count.clips(GRanges())), 0)
})



## varbase
test_that('varbase', {
    ## ## Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.
    ## expect_error(varbase('foobar'))
    ## ## default
    ## expect_true(inherits(varbase(read.bam(example_bam, all=TRUE)[[1]]), 'GRangesList'))
    ## expect_equal(length(varbase(read.bam(example_bam, all=TRUE)[[1]])), 2)
    ## expect_match(as.character(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$varbase)[1], 'TAACCCCAACCAAAACCGCCCAACCCTAA')
    ## expect_match(as.character(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$varbase)[2], 'C')
    ## expect_equal(as.integer(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$varlen[1]), 29)
    ## expect_equal(as.integer(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$varlen[2]), 1)
    ## expect_equal(as.character(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$type)[1], 'S')
    ## expect_equal(as.character(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$type)[2], 'X')
    ## expect_equal(length(varbase(read.bam(example_bam, all=TRUE)[[1]])), 2)
    ## expect_match(as.character(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$col)[1], '#FFC0CBE6')
    ## expect_match(as.character(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$col)[2], 'blue')
    ## expect_equal(length(varbase(read.bam(example_bam, all=TRUE)[[1]])), 2)
    ## expect_match(as.character(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$border)[1], '#FFC0CBE6')
    ## expect_match(as.character(varbase(read.bam(example_bam, all=TRUE)[[1]])[[1]]$border)[2], 'blue')
    ## ## soft
    ## expect_match(varbase(read.bam(example_bam, all=TRUE)[[1]], soft=FALSE)[[1]]$varbase, 'C')
    ## expect_equal(as.integer(varbase(read.bam(example_bam, all=TRUE)[[1]], soft=FALSE)[[1]]$varlen), 1)
    ## expect_match(varbase(read.bam(example_bam, all=TRUE)[[1]], soft=FALSE)[[1]]$type, 'X')
    ## expect_match(varbase(read.bam(example_bam, all=TRUE)[[1]], soft=FALSE)[[1]]$col, 'blue')
    ## expect_match(varbase(read.bam(example_bam, all=TRUE)[[1]], soft=FALSE)[[1]]$border, 'blue')
    ## ## verbose
    ## expect_equal(length(varbase(read.bam(example_bam, all=TRUE)[[1]], verbose=FALSE)), 2)
    ## ## if (inherits(reads, 'GRangesList')){
    ## expect_equal(length(varbase(GRangesList(read.bam(example_bam, all=TRUE, intervals = GRanges('1:5075-18800'))))[[2]]), 10)
    ## ## else if (inherits(reads, 'data.frame')){
    ## expect_equal(length(varbase(gr2dt(read.bam(example_bam, all=TRUE, intervals = GRanges('1:5075-18800'))[[2]]))[[1]]), 5)
    ## expect_equal(varbase(gr2dt(read.bam(example_bam, all=TRUE, intervals = GRanges('1:5075-18800'))[[2]]))[[1]][2]$varbase, 'C')

})



## gr2dt(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))


## as.data.frame(gr2dt(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[2]]))


test_that('splice.cigar', {
    ## check input
    ## expect_error(splice.cigar('foo'))
    ## expect_error(splice.cigar(example_bam)) ##  Error: Reads must be either GRanges, GRangesList, GappedAlignments, or data.table object. Please see documentation for details.
    ## ## check default
    ## expect_equal(length(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))), 2054)
    ## expect_equal(length(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]), 9)
    ## expect_equal(length(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[2]]), 0)
    ## expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$type[1], 'S')
    ## expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$type[2], 'M')
    ## expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$type[3], 'I')
    ## expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$type[4], 'M')
    ## expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]])[1], 0)
    ## expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]])[2], 6)
    ## expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]])[3], 0)
    ## expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]])[4], 63)
    ## expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$fid[1], 1)
    ## expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$fid[2], 1)
    ## expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))[[1]]$qname[1], 'C1Y1JACXX130321:7:1101:15151:82244')
    ## ## fast = FALSE
    ## expect_equal(length(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast = FALSE)[[2]]), 0)
    ## expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]])[1], 6)
    ## expect_equal(width(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]])[2], 0)
    ## expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$type[1], 'M')
    ## expect_match(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$type[2], 'I')
    ## expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$fid[1], 1)
    ## expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$fid[2], 1)
    ## expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$rid[1], 1)
    ## expect_equal(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), fast=FALSE)[[1]]$rid[2], 1)
    ## ## use.D = FALSE
    ## expect_true('riid' %in% colnames(as.data.frame(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), use.D = FALSE)[[1]])))  ## should be a 'riid' column in GRanges
    ## expect_true('riid' %in% colnames(as.data.frame(splice.cigar(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')), use.D = FALSE)[[2]])))
    ## ## tests for 'rem.soft' and 'get.seq'
    ## ## check 'if (nreads==0){'
    ## expect_true(inherits(splice.cigar(GRanges()), 'GRangesList'))
    ## expect_equal(length(splice.cigar(GRanges())), 0)
    ## ##  if (inherits(reads, 'GRangesList')){
    ## expect_equal(splice.cigar(GRangesList(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))))[[3]]$type, 'M')
    ## ## if (is.data.frame(reads)){
    ## ##   unable to find an inherited method for function ‘values’ for signature ‘"data.table"’
    ## ## if (inherits(reads, 'GRangesList')){
    ## expect_equal(length(splice.cigar(reads= read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100')))), 2054)
    ## ## GRangesList
    ## ### Error: Reads must have cigar and seq fields specified. Please see documentation for details.
    ## expect_error(splice.cigar(grl2))
    ## data.frame
    ## expect_equal(length(splice.cigar(as.data.frame(gr2dt(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))


})




test_that('bamflag', {
    ## ## isPaired  1 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isPaired[1], 1)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isPaired[2], 0)
    ## ## isProperPair  1 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isProperPair[1], 1)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isProperPair[2], 0)
    ## ## isUnmappedQuery  0 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isUnmappedQuery[1], 0)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isUnmappedQuery[2], 0)
    ## ## hasUnmappedMate  0 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$hasUnmappedMate[1], 0)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$hasUnmappedMate[2], 0)
    ## ## isMinusStrand  1 1
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMinusStrand[1], 0)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMinusStrand[2], 0)
    ## ## isMateMinusStrand  1 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMateMinusStrand[1], 1)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isMateMinusStrand[2], 0)
    ## ## isFirstMateRead  0 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isFirstMateRead[1], 0)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isFirstMateRead[2], 0)
    ## ## isSecondMateRead  1 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSecondMateRead[1], 1)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSecondMateRead[2], 0)
    ## ## isNotPrimaryRead  0 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPrimaryRead[1], 0)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPrimaryRead[2], 0)
    ## ## isNotPassingQualityControls  0 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPassingQualityControls[1], 0)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isNotPassingQualityControls[2], 0)
    ## ## isDuplicate  0 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isDuplicate[1], 0)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isDuplicate[2], 0)
    ## ## isSupplementary  0 0
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSupplementary[1], 0)
    ## expect_equal(as.data.frame(bamflag(read.bam(example_bam, all=TRUE, intervals = GRanges('1:10075-10100'))[[1]]))$isSupplementary[2], 0)

})




test_that('bamtag', {
    ## ## check errors
    ## expect_error(bamtag('foo'))   ### Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.
    ## expect_error(bamtag(example_bam))    ### Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.
    ## ## default
    ## ### FAILS AT GRANGELIST: WHY? cf read.bam()
    ## ## expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100')))[1], 'ST-K00126:3:H5TL3BBXX:1:1127:17310:39893_1__')
    ## expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]])[1], 'C1Y1JACXX130321:7:1101:15151:82244_2__')
    ## ## 'secondary' == TRU#
    ## expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], secondary = TRUE)[1], 'C1Y1JACXX130321:7:1101:15151:82244_2__0')
    ## ## 'gr.string' == TRUE
    ## expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], gr.string = TRUE)[1], 'C1Y1JACXX130321:7:1101:15151:82244_2_1:10034-10117+_')
    ## expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], gr.string = TRUE)[2], 'C1Y1JACXX130321:7:1101:15151:82244_2_1:10160-10260-_')
    ## ## 'secondary' == TRUE & 'gr.string' == TRUE
    ## expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], secondary = TRUE, gr.string = TRUE)[1], 'C1Y1JACXX130321:7:1101:15151:82244_2_1:10034-10117+_0')
    ## expect_equal(bamtag(read.bam(example_bam, all = TRUE, intervals = GRanges('1:10075-10100'))[[1]], secondary = TRUE, gr.string = TRUE)[2], 'C1Y1JACXX130321:7:1101:15151:82244_2_1:10160-10260-_0')

})

test_that('countCigar', {
  ## cigar = c('43S58M', '66S35M', '31S70M', NA, '58S43M', '24M1D30M47S', '46M55S', '101M', '51S50M', '28M1I28M44S', '27S69M1D5M')
  ## res = countCigar(cigar)
  ## expect_equal(as.vector(res), c(0, 0, 0, NA, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, NA, 0, 0, 0, 0, 0, 1, 0, 58, 35, 70, NA, 43, 54, 46, 101, 50, 56, 74, 43, 66, 31, NA, 58, 47, 55, 0, 51, 44, 27, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0))
  ## expect_equal(dim(res), c(11, 5))
  ## expect_equal(colnames(res), c('D', 'I', 'M', 'S', 'N')) 
})


test_that('is.paired.end', {
    ## expect_true(as.logical(is.paired.end(example_bam)))
    ## expect_true(as.logical(is.paired.end(small_MD_bam)))
    ## expect_equal(as.logical(is.paired.end('foo')), NA)   ### error checking, should return NA
    ## expect_equal(as.logical(is.paired.end(NA)), NA)
})

test_that('chunk', {
    ## expect_equal(dim(chunk(2, 10, 1, length.out=4))[1], 4)
    ## expect_equal(dim(chunk(2, 10, 1, length.out=4))[2], 2)
    ## ## to = NULL
    ## expect_equal(chunk(2, NULL, 1, length.out=4)[, 1], c(1, 1, 1, 2))
    ## expect_equal(chunk(2, NULL, 1, length.out=4)[, 2], c(0, 0, 1, 2))
})

### varcount
## min.mapq integer Minimal mapping quality at which to compute bases
## max.baseq integer Minimal base qualitya t which to compute bases
## max.depth integer Maximum read depth to consider
## indel boolean Flag whether to consider indels (default FALSE)
test_that('varcount', {

    ## expect_equal(varcount(small_MD_bam, gr= GRanges('1:10075-10100'))$counts[1], 17)
    ## expect_equal(varcount(small_MD_bam, gr= GRanges('1:10075-10100'))$counts[2], 0)
    ## expect_equal(varcount(small_MD_bam, gr= GRanges('1:10075-10100'))$counts[3], 0)
    ## expect_equal(varcount(small_MD_bam, gr= GRanges('1:10075-10100'))$counts[4], 0)
    ## expect_equal(varcount(small_MD_bam, gr= GRanges('1:10075-10100'))$counts[5], 0)
    ## expect_equal(width(varcount(small_MD_bam, gr= GRanges('1:10075-10100'))$gr), 1)
    ## expect_equal(varcount(example_bam, gr= GRanges('1:10075-10100'))$counts[1], 511)
    ## expect_equal(varcount(example_bam, gr= GRanges('1:10075-10100'))$counts[2], 6)
    ## expect_equal(varcount(example_bam, gr= GRanges('1:10075-10100'))$counts[3], 4)
    ## expect_equal(varcount(example_bam, gr= GRanges('1:10075-10100'))$counts[4], 0)
    ## expect_equal(varcount(example_bam, gr= GRanges('1:10075-10100'))$counts[5], 0)
    ## ## min.mapq
    ## expect_equal(varcount(small_MD_bam, gr= GRanges('1:10075-10100'), min.mapq = 45)$counts[1], 8)
    ## ## min.baseq
    ## expect_equal(varcount(example_bam, gr= GRanges('1:10075-10100'), min.baseq = 1)$counts[1], 529)
    ## expect_equal(varcount(example_bam, gr= GRanges('1:10075-10100'), min.baseq = 1)$counts[2], 12)
    ## ## max.depth
    ## expect_equal(varcount(example_bam, gr= GRanges('1:10075-10100'), max.depth = 1)$counts[1], 75)
    ## expect_equal(varcount(example_bam, gr= GRanges('1:10075-10100'), max.depth = 100)$counts[1], 160)
    ## ## indel
    ## expect_equal(varcount(small_MD_bam, gr= GRanges('1:10075-10100'), indel=TRUE)$gr, NULL)
    ## ## if (any(!file.exists(bami)))
    ## ##  Error: one or more BAM file indices missing
    ## expect_error(varcount(noindexbam, gr= GRanges('1:10075-10100')))
    ## ## try two BAMs
    ## expect_equal(length(varcount(c(normalbam, tumorbam), gr= GRanges('1:10075-10100'))$counts), 10)

})



## read_vcf()
## read_vcf = function(fn, gr = NULL, hg = 'hg19', geno = NULL, swap.header = NULL, verbose = FALSE, add.path = FALSE, tmp.dir = '~/temp/.tmpvcf', ...)
##test_that('read_vcf', {
#    ## error
#    expect_error(read_vcf('foobar'))
#    ## default
#    expect_equal(length(read_vcf(somatic_vcf)), 60)
#    expect_equal(length(seqnames(seqinfo(read_vcf(somatic_vcf)))), 84)
#    ## gr  gr= GRanges('1:10075-10100')
#    ## hg
##    expect_match(unique(as.data.frame(seqinfo(read_vcf(somatic_vcf, hg='hg12345')))$genome), 'hg12345')
#    ## geno
#    ## swap.header
#    expect_equal(length(seqnames(seqinfo(read_vcf(somatic_vcf, swap.header='/Users/ebiederstedt/bamUtils/tests/testthat/new_header.vcf')))), 2)
#    ## verbose
#    expect_equal(length(read_vcf(somatic_vcf, verbose=TRUE)), 60)
#    ## check 'if (!file.exists(swap.header))'
#    expect_error(read_vcf(somatic_vcf, swap.header='foobar'))
#
#})



##test_that('write_vcf', {
    ##expect_error(write_vcf(read_vcf(somatic_vcf), filename = './foo.vcf'), NA)  ### just check it runs
    ##
##})



test_that('mafcount', {

    ## expect_equal(length(mafcount(tumorbam, chunk.size = 1e5, maf = dt2gr(fread(mafpath)))), 54103)
    ## ## include normal BAM
    ## ## if (!is.null(norm.bam)){
    ## expect_equal(length(mafcount(tumorbam, normalbam, chunk.size = 1e5, maf = dt2gr(fread(mafpath)))), 54103)
    ## ##  if (is.data.frame(maf)){
    ## ##  strange error here...
    ## ## expect_equal(length(mafcount(tumorbam, normalbam, chunk.size = 1e5, maf = as.data.frame(dt2gr(fread(mafpath))))), 54103)
    ## ## if (is.null(maf$Tumor_Seq_Allele1)){
    ## ## if (is.null(maf$Tumor_Seq_Allele1)){
    ## ## if (is.null(maf$Reference_Allele)){
    ## ## if (is.null(maf$Reference_Allele)){
    ## ## if (!all(is.character(maf$Tumor_Seq_Allele1))){
    ## ## if (!all(is.character(maf$Reference_Allele))){
    ## ## if (is.null(maf$Reference_Allele) | is.null(maf$Tumor_Seq_Allele1)){
    ## maffoo = dt2gr(fread(mafpath))
    ## maffoo$Tumor_Seq_Allele1 = NULL
    ## maffoo$Reference_Allele = NULL
    ## expect_error(mafcount(tumorbam, chunk.size = 1e5, maf = maffoo))
    ## ## Error: Cannot locate variant columns in input GRanges, please check input to make sure it either has standard VCF ALT / REF columns or MAF file columns specifying alt and ref allele


})





