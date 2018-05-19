#' @import GenomicRanges
#' @import GenomicAlignments
#' @import data.table
#' @import Rsamtools
#' @import gUtils

 
#' @name read.bam
#' @title Read BAM file into GRanges or data.table
#' @description 
#'
#' Wrapper around Rsamtools BAM scanning functions
#' By default, returns GRangesList of read pairs for which <at least one> read lies in the supplied interval
#' 
#' @param bam string Input BAM file. Advisable to make BAM a BamFile instance instead of a plain string, so that the index does not have to be reloaded.
#' @param bai string Input BAM index file.
#' @param intervals GRanges of intervals to retrieve. If left unspecified with 'all = TRUE', will try to pull down entire BAM file  
#' @param gr Granges (default = intervals)
#' @param all boolean Flag to read in all of BAM as a GRanges via `si2gr(seqinfo())` (default = FALSE)
#' @param pairs.grl boolean Flag if TRUE will return GRangesList of read pairs in which at least one read falls in the supplied interval (default = FALSE)
#' @param stripstrand boolean Flag to ignore strand information on the query intervals (default = TRUE)
#' @param what vector What fields to pull down from BAM. (default = \code{scanBamWhat()})
#' @param verbose boolean verbose flag (default = FALSE)
#' @param tag vector Additional tags to pull down from the BAM (e.g. 'R2')
#' @param isPaired boolean Flag indicates whether unpaired (FALSE), paired (TRUE), or any (NA) read should be returned. See documentation for Rsamtools::scanBamFlag(). (default = NA)
#' @param isProperPair boolean Flag indicates whether improperly paired (FALSE), properly paired (TRUE), or any (NA) read should be returned. A properly paired read is defined by the alignment algorithm and might, e.g., represent reads aligning to identical reference sequences and with a specified distance. See documentation for Rsamtools::scanBamFlag(). (default = NA)
#' @param isUnmappedQuery boolean Flag indicates whether unmapped (TRUE), mapped (FALSE), or any (NA) read should be returned. See documentation for Rsamtools::scanBamFlag(). (default = NA)
#' @param hasUnmappedMate boolean Flag indicates whether reads with mapped (FALSE), unmapped (TRUE), or any (NA) mate should be returned. See documentation for Rsamtools::scanBamFlag(). (default = NA)
#' @param isNotPassingQualityControls boolean Flag indicates whether reads passing quality controls (FALSE), reads not passing quality controls (TRUE), or any (NA) read should be returned. See documentation for Rsamtools::scanBamFlag(). (default = NA)
#' @param isDuplicate boolean Flag indicates that un-duplicated (FALSE), duplicated (TRUE), or any (NA) reads should be returned. 'Duplicated' reads may represent PCR or optical duplicates. See documentation for Rsamtools::scanBamFlag(). (default = FALSE)
#' @param pairs.grl.split boolean Return reads as GRangesList. Controls whether get.pairs.grl() does split (default = TRUE)
#' @param as.data.table boolean Return reads in the form of a data.table rather than GRanges/GRangesList (default = FALSE)
#' @param ignore.indels boolean Flag messes with cigar to read BAM with indels removed. Useful for breakpoint mapping on contigs (default = FALSE)
#' @param ... futher arguments passed to Rsamtools::scanBamFlag()
#' @return Reads in one of GRanges, GRangesList or data.table
#' @export
read.bam = function(bam, intervals = NULL, gr = intervals, all = FALSE,
                    bai = NULL,
                    pairs.grl = TRUE,   ## if TRUE will return GRangesList of read pairs for whom at least one read falls in the supplied interval
                                        ##  paired = F, # if TRUE, will used read bam gapped alignment pairs warning: will throw out pairs outside of supplied window
                                        ##  gappedAlignment = T, # if false just read alignments using scanbam
                    stripstrand = TRUE,
                    what = scanBamWhat(),
                    unpack.flag = FALSE,  ## will add features corresponding to read flags
                    verbose = FALSE,
                    tag = NULL,
                    isPaired = NA, ## if these features are NA, then reads satisfying both T and F will be returned
                    isProperPair = NA,
                    isUnmappedQuery = NA,
                    hasUnmappedMate = NA,
                    isNotPassingQualityControls = NA,
                    isDuplicate = FALSE,
                    isValidVendorRead = TRUE,
                    pairs.grl.split = TRUE,  ## return pairs as grl, rather than GRanges .. controls whether get.pairs.grl does split (t/c rename to pairs.grl.split)
                    as.data.table = FALSE, ## returns reads in the form of a data table rather than GRanges/GRangesList
                    ignore.indels = FALSE, ## messes with cigar to read BAM with indels removed. Useful for breakpoint mapping on contigs
                    size.limit = 1e6,
                    ... ## passed to scanBamFlag (
                    )
{

    ## check that the BAM is valid
    check_valid_bam = readChar(gzfile(bam, 'r'), 4)
    if (!identical(check_valid_bam, 'BAM\1')){
        stop("Cannot open BAM. A valid BAM for 'bam_file' must be provided.")
    }

    if (!inherits(bam, 'BamFile'))
    {
        if (is.null(bai))
        {
            if (file.exists(bai <- gsub('.bam$', '.bai', bam))){
                bam = BamFile(bam, bai)
            } else if (file.exists(bai <- paste(bam, '.bai', sep = ''))){
                bam = BamFile(bam, bai)
            } else{
                bam = BamFile(bam)
            }
        } else{
            bam = BamFile(bam, index = bai)
        }
    }
    ## if intervals unspecified will try to pull down entire bam file (CAREFUL)
    if (length(intervals)==0){
        intervals = NULL
    }

    if (is.null(intervals)){
        intervals = gr
    }

    if (is.null(intervals))
    {
        if (all){
            intervals = si2gr(seqinfo(bam))
        } else{
            stop('Error: Must provide non empty interval list')
        }
    }

    if (class(intervals) == 'data.frame'){
        intervals = seg2gr(intervals);
    }

    if (inherits(intervals, 'GRangesList')){
        intervals = unlist(intervals);
    }

    if (stripstrand){
        strand(intervals) = '*'
    }

    intervals = reduce(intervals);

    now = Sys.time();

    if (pairs.grl){
        paired = FALSE
    }

    flag = scanBamFlag(isPaired = isPaired, isProperPair = isProperPair, isUnmappedQuery = isUnmappedQuery,
                       hasUnmappedMate = hasUnmappedMate, isNotPassingQualityControls = isNotPassingQualityControls,
                       isDuplicate = isDuplicate, ...)

    tag = unique(c('MD', 'MQ', tag))

    param = ScanBamParam(which = gr.fix(intervals, bam, drop = TRUE), what = what, flag = flag, tag = tag)

    if (verbose){
        cat('Reading bam file\n')
    }

  if (class(bam) == 'BamFile'){
        out <- suppressWarnings(scanBam(bam, param=param))
    } else{
        out <- suppressWarnings(scanBam(bam, index=bai, param=param))
    }

    if (verbose) {
        print(Sys.time() - now)
        print('BAM read. Making into data.frame')
    }

    out <- out[sapply(out, function(x) length(x$qname)>0)]
    ## names(out[[1]])
    ## [1] "qname"  "flag"   "rname"  "strand" "pos"    "qwidth" "mapq"   "cigar" 
    ## [9] "mrnm"   "mpos"   "isize"  "seq"    "qual"   "tag"  

    if (length(out)>0)
    {
        if (verbose) {
            print(Sys.time() - now)
            print('combining lists')
        }

        out <- as.data.frame(rbindlist(lapply(out, function(x){
            
            x <- c(x[-match('tag', names(x))], x$tag)

            x <- x[sapply(x, length)>0]
            conv <- which(!(sapply(x, class) %in% c('integer', 'numeric', 'character')))
            x[conv] <- lapply(x[conv], as.character)

            for (t in tag){
                if (!(t %in% names(x))){
                    x[[t]] = rep(NA, length(x$qname))
                }
            }

            if (!('R2' %in% names(x)) && 'R2' %in% tag){
                x$R2 <- rep(NA, length(x$qname))
            }
            if (!('Q2' %in% names(x)) && 'Q2' %in% tag){
                x$Q2 <- rep(NA, length(x$qname))
            }
            x
        })))

        ## faster CIGAR string parsing with vectorization and data tables
        if (verbose) {
            print(Sys.time() - now)
            print('filling pos2 from cigar')
        }
        if (ignore.indels) {
          cigar <- gsub('[0-9]+D', '', gsub('([0-9]+)I', '\\1M', out$cigar))  ## Remove deletions, turn insertions to matches
          cig = rep(list(), length(cigar))
          if (any(ix <- !is.na(cigar)))
              {
                cig[ix] <- explodeCigarOps(cigar[ix])        # formerly `cig <- splitCigar(cigar)`, splitCigar() now deprecated
              }
            torun=sapply(cig, function(y) any(duplicated((y[[1]][y[[1]]=="M"]))))
            new.cigar <- sapply(cig[torun], function(y) {
                lets <- y[[1]][!duplicated(y[[1]])]
                vals <- y[[2]][!duplicated(y[[1]])]
                vals[lets=="M"] <- sum(y[[2]][y[[1]]=="M"])
                lets <- strsplit(rawToChar(lets), '')[[1]]
                paste(as.vector(t(matrix(c(vals, lets), nrow=length(vals), ncol=length(lets)))), collapse='')
            })
          if (any(torun)){
            out$cigar[torun] <- new.cigar
          }

        }


        cigs <- countCigar(out$cigar)
        # out$pos2 <- out$pos + cigs[, "M"]
        out$pos2 <- out$pos + rowSums(cigs[, c("D", "M")], na.rm=T) - 1

        if (verbose) {
            print(Sys.time() - now)
            print('fixing seqdata')
        }
        out$qwidth = nchar(out$seq)
        unm = is.na(out$pos)
        if (any(unm)){

            out$pos[unm] = 1
            out$pos2[unm] = 0
            out$strand[unm] = '*'
        }
        gr.fields = c('rname', 'strand', 'pos', 'pos2');
        out = as.data.table(out)
        vals = out[, setdiff(names(out), gr.fields), with=FALSE]

        if (!as.data.table) {
            out <- GRanges(out$rname, IRanges(out$pos, pmax(0, out$pos2-1)), strand = out$strand, seqlengths = seqlengths(intervals))
            values(out) <- vals;
        } else {
            out <- data.table(seqnames=out$rname, start=out$pos, end= pmax(out$pos2-1, 0), strand=out$strand)
            val <- data.table(vals)
            out <- cbind(out, val)
        }
    } else {
        if (!as.data.table){
            return(GRanges(seqlengths = seqlengths(intervals)))
        } else{
            return(data.table())
        }
    }

    if (verbose)
    {
        if (as.data.table){
            cat(sprintf('Extracted %s reads\n', nrow(out)))
        } else{
            cat(sprintf('Extracted %s reads\n', length(out)))
        }
        print(Sys.time() - now)
    }

    if (pairs.grl){

        if (verbose){
            cat('Pairing reads\n')
        }
        out = get.pairs.grl(out, pairs.grl.split = pairs.grl.split)
        if (verbose){
            cat('done\n')
            print(paste0('Total time to complete: ', Sys.time() - now))
        }
        if (pairs.grl.split && !as.data.table){
            names(out) = NULL;
            values(out)$col = 'gray';
            values(out)$border = 'gray';
        }
    }

    return(out)
}




#' @name bam.cov.gr
#' @title Get coverage as GRanges from BAM on custom set of GRanges
#' @description
#'
#' Gets coverage from BAM in supplied GRanges using 'countBam()', returning GRanges with coverage counts in
#' each of the provided GRanges (different from 'bamUtils::bam.cov()') specified as the 
#' columns $file, $records, and $nucleotides in the values field
#'
#' Basically a wrapper for 'Rsamtools::countBam()' with some standard settings for 'Rsamtools::ScanBamParams()'
#'
#' @param bam string Input BAM file. Advisable to make the input BAM a BamFile instance instead of a plain string, so that the index does not have to be reloaded.
#' @param bai string Input BAM index file
#' @param intervals GRanges of intervals to retrieve
#' @param all boolean Flag to read in all of BAM as a GRanges via `si2gr(seqinfo())` (default = FALSE)
#' @param isPaired boolean Flag indicates whether unpaired (FALSE), paired (TRUE), or any (NA) read should be returned. See documentation for Rsamtools::scanBamFlag(). (default = NA)
#' @param isProperPair boolean Flag indicates whether improperly paired (FALSE), properly paired (TRUE), or any (NA) read should be returned. A properly paired read is defined by the alignment algorithm and might, e.g., represent reads aligning to identical reference sequences and with a specified distance. See documentation for Rsamtools::scanBamFlag(). (default = NA)
#' @param isUnmappedQuery boolean Flag indicates whether unmapped (TRUE), mapped (FALSE), or any (NA) read should be returned. See documentation for Rsamtools::scanBamFlag(). (default = NA)
#' @param hasUnmappedMate boolean Flag indicates whether reads with mapped (FALSE), unmapped (TRUE), or any (NA) mate should be returned. See documentation for Rsamtools::scanBamFlag(). (default = NA)
#' @param isNotPassingQualityControls boolean Flag indicates whether reads passing quality controls (FALSE), reads not passing quality controls (TRUE), or any (NA) read should be returned. See documentation for Rsamtools::scanBamFlag(). (default = FALSE)
#' @param isDuplicate boolean Flag indicates that un-duplicated (FALSE), duplicated (TRUE), or any (NA) reads should be returned. 'Duplicated' reads may represent PCR or optical duplicates. See documentation for Rsamtools::scanBamFlag(). (default = FALSE)
#' @param mc.cores integer Number of cores in mclapply (default = 1)
#' @param chunksize integer How many intervals to process per core (default = 10)
#' @param verbose boolean "verbose" flag (default = FALSE)
#' @param ... futher arguments passed into Rsamtools::scanBamFlag()
#' @return GRanges parallel to input GRanges, but with metadata filled in.
#' @export
bam.cov.gr = function(bam, bai = NULL, intervals = NULL, all = FALSE, count.all = FALSE, isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = FALSE, 
    hasUnmappedMate = FALSE, isNotPassingQualityControls = FALSE, isDuplicate = FALSE, mc.cores = 1, chunksize = 10, verbose = FALSE, ...)
{

    ## check that the BAM is valid
    check_valid_bam = readChar(gzfile(bam, 'r'), 4)
    if (!identical(check_valid_bam, 'BAM\1')){
        stop("Cannot open BAM. A valid BAM for 'bam_file' must be provided.")
    }

    if (missing(bam) | missing(intervals)){
        stop("Error: arguments 'bam' and 'intervals' are both required for 'bam.cov.gr'. Please see documentation for details.")
    }
    if (!is(intervals, "GRanges")){
        stop("Error: Granges of intervals to retrieve 'intervals' must be in the format 'GRanges'. Please see documentation for details.")
    }

    if (is.character(bam)){
        if (!is.null(bai)){
            bam = BamFile(bam, bai)
        } else{
            if (file.exists(paste(bam, 'bai', sep = '.'))){
                bam = BamFile(bam, paste(bam, 'bai', sep = '.'))
            } else if (file.exists(gsub('.bam$', '.bai', bam))){
                bam = BamFile(bam, paste(bam, 'bai', sep = '.'))
            } else{
                stop('Error: BAM index not found, please find index and specify BAM file argument as valid BamFile object. Please see documentation for details.')
            }
        }
    }


    keep = which(as.character(seqnames(intervals)) %in% seqlevels(bam))
    
    if (length(keep) > 0){
      ix = c(keep[c(seq(1, length(keep), chunksize))], keep[length(keep)] + 1);  ## prevent bam error from improper chromosomes
        chunk.id = unlist(lapply(1:(length(ix) - 1), function(x) rep(x, ix[x+1] - ix[x])))

        gr.chunk = split(intervals[keep], chunk.id[keep]);
        if (count.all){
            flag = scanBamFlag()
        } else{
            flag = scanBamFlag(isPaired = isPaired, isProperPair = isProperPair, isUnmappedQuery = isUnmappedQuery,
                               hasUnmappedMate = hasUnmappedMate, isNotPassingQualityControls = isNotPassingQualityControls,
                               isDuplicate = isDuplicate, ...)
        }

        out = rbindlist(mclapply(1:length(gr.chunk), function(x) {
            if (verbose){
                cat(sprintf('Processing ranges %s to %s of %s, extracting %s bases\n', ix[x], ix[x+1]-1, length(keep), sum(width(gr.chunk[[x]]))))
            }
            as.data.table(countBam(bam, param = ScanBamParam(which = dt2gr(gr2dt(gr.chunk[[x]]), seqlengths = seqlengths(bam)), flag = flag)))
        }, mc.cores = mc.cores));


      gr.tag = paste(as.character(seqnames(intervals)), start(intervals), end(intervals));
      out.tag = paste(out$space, out$start, out$end);
      ix = match(gr.tag, out.tag);
      values(intervals) = cbind(as.data.frame(values(intervals)), out[ix, c('file', 'records', 'nucleotides'), with = FALSE])
    } else{
        values(intervals) = cbind(as.data.frame(values(intervals)), data.frame(file = rep(gsub('.*\\/([^\\/]+)$', '\\1', path(bam)), length(intervals)), records = NA, nucleotides = NA))
    }

    return(intervals)
}





#' @name bam.cov.tile
#' @title Get coverage as GRanges from BAM on genome tiles across seqlengths of genome
#' @description
#'
#' Quick way to get tiled coverage via piping to samtools (e.g. ~10 CPU-hours for 100bp tiles, 5e8 read pairs)
#'
#' Gets coverage for window size "window" [bp], pulling "chunksize" records at a time and incrementing bin
#' corresponding to midpoint or overlaps of corresponding (proper pair) fragment (uses TLEN and POS for positive strand reads that are part of a proper pair)
#'
#' @param bam.file string Input BAM file
#' @param window integer Window size (in bp) (default = 1e2)
#' @param chunksize integer Size of window (default = 1e5)
#' @param min.mapq integer Minimim map quality reads to consider for counts (default = 30)
#' @param verbose boolean Flag to increase vebosity (default = TRUE)
#' @param max.tlen integer Maximum paired-read insert size to consider (default = 1e4)
#' @param st.flag string Samtools flag to filter reads on (default = '-f 0x02 -F 0x10')
#' @param fragments boolean flag (default = FALSE) detremining whether to compute fragment (i.e. proper pair footprint aka insert) density or read density
#' @param do.gc boolean Flag to execute garbage collection via 'gc()' (default = FALSE)
#' @param midpoint boolean Flag if TRUE will only use the fragment midpoint, if FALSE will count all bins that overlap the fragment (default = TRUE)
#' @return GRanges of "window" bp tiles across seqlengths of bam.file with meta data field $counts specifying fragment counts centered (default = TRUE)
#' in the given bin.
#' @export
bam.cov.tile = function(bam.file, window = 1e2, chunksize = 1e5, min.mapq = 30, verbose = TRUE, max.tlen = 1e4, 
                        st.flag = "-f 0x02 -F 0x10", fragments = TRUE, do.gc = FALSE, midpoint = TRUE, bai = NULL)
{

    ## check that the BAM is valid
    check_valid_bam = readChar(gzfile(bam.file, 'r'), 4)
    if (!identical(check_valid_bam, 'BAM\1')){
        stop("Cannot open BAM. A valid BAM for 'bam_file' must be provided.")
    }

    cmd = 'samtools view %s %s -q %s | cut -f "3,4,9"' ## cmd line to grab the rname, pos, and tlen columns

    sl = seqlengths(BamFile(bam.file))

    counts = lapply(sl, function(x) rep(0, ceiling(x/window)))
    numwin = sum(sapply(sl, function(x) ceiling(x/window)))

    cat('Calling', sprintf(cmd, st.flag, bam.file, min.mapq), '\n')
    p = pipe(sprintf(cmd, st.flag, bam.file, min.mapq), open = 'r')

    i = 0
    sl.dt = data.table(chr = names(sl), len = sl)
    counts = sl.dt[, list(start = seq(1, len, window)), by = chr]
    counts = counts[, bin := 1:length(start), by = chr]
    counts[, end := pmin(start + window-1, sl[chr])]
    counts[, count := 0]
    counts[, rowid := 1:length(count)]
    setkeyv(counts, c("chr", "bin")) ## now we can quickly populate the right entries
    totreads = 0

    st = Sys.time()
    if (verbose){
        cat('Starting fragment count on', bam.file, 'with bin size', window, 'and min mapQ', min.mapq, 'and   insert size limit', max.tlen, 'with midpoint set to', midpoint, '\n')
    }

    while (length(chunk <- readLines(p, n = chunksize)) > 0)
    {
        i = i+1

        if (fragments){

            chunk = fread(paste(chunk, collapse = "\n"), header = F)[abs(V3) <= max.tlen, ]  ## only take midpoints

            if (midpoint){
                ## only take midpoints
                chunk[, bin := 1 + floor((V2 + V3/2)/window)] ## use midpoint of template to index the correct bin

            } else {
                ## enumerate all bins containing fragment i.e. where fragments overlap multiple bins  (slightly slower)
                if (verbose){
                    cat('Counting all overlapping bins!\n')
                }
                chunk[, ":="(bin1 = 1 + floor((V2)/window), bin2 = 1 + floor((V2+V3)/window))]
                chunk = chunk[, list(V1, bin = bin1:bin2), by = list(ix = 1:length(V1))]
            }
        } else{   
            ## just count reads
            cat('Just counting reads\n')
            chunk = fread(paste(chunk, collapse = "\n"), header = F)
            chunk[, bin := 1 + floor((V2)/window)]
        }

        tabs = chunk[, list(newcount = length(V1)), by = list(chr = as.character(V1), bin)] ## tabulate reads to bins data.table style
        counts[tabs, count := count + newcount] ## populate latest bins in master data.table

        ## should be no memory issues here since we preallocate the data table .. but they still appear
        if (do.gc){
            print('GC!!')
            print(gc())
        }

        ## report timing
        if (verbose){

            cat('bam.cov.tile.st ', bam.file, 'chunk', i, 'num fragments processed', i*chunksize, '\n')
            timeelapsed = as.numeric(difftime(Sys.time(), st, units = 'hours'))
            meancov = i * chunksize / counts[tabs[nrow(tabs),], ]$rowid  ## estimate comes from total reads and "latest" bin filled
            totreads = meancov * numwin
            tottime = totreads*timeelapsed/(i*chunksize)
            rate = i*chunksize / timeelapsed / 3600
            cat('mean cov:', round(meancov,1), 'per bin, estimated tot fragments:', round(totreads/1e6,2), 'million fragments, processing', rate,
                'fragments/second\ntime elapsed:', round(timeelapsed,2), 'hours, estimated time remaining:', round(tottime - timeelapsed,2), 'hours', ', estimated total time', round(tottime,2), 'hours\n')
        }
    }

    gr = GRanges(counts$chr, IRanges(counts$start, counts$end), count = counts$count, seqinfo = Seqinfo(names(sl), sl))
    if (verbose){
        cat("Finished computing coverage, and making GRanges\n")
    }
    close(p)


    return(gr)
}




#' @name get.mate.gr
#' @title returns GRanges corresponding to mates of reads
#' @description
#'
#' Inputs GRanges or data.frame/data.table of reads. Outputs GRanges corresponding to mates of reads.
#'
#' @param reads GRanges or data.table/data.frame Input reads
#' @return \code{GRanges} corresponding to mates of reads
#' @name get.mate.gr
#' @export
get.mate.gr = function(reads)
{
    if (inherits(reads, 'GRanges')) {
        mpos = values(reads)$mpos
        mrnm = as.vector(values(reads)$mrnm)
        mapq = values(reads)$MQ
        bad.chr = !(mrnm %in% seqlevels(reads)); ## these are reads mapping to chromosomes that are not in the current "genome"
        mrnm[bad.chr] = as.character(seqnames(reads)[bad.chr]) # we set mates with "bad" chromosomes to have 0 width and same seqnames (i.e. as if unmapped)
    } else if (inherits(reads, 'data.table')) {
        mpos <- reads$mpos
        mrnm <- reads$mrnm
        mapq = reads$MQ
        bad.chr <- !(mrnm %in% c(seq(22), 'X', 'Y', 'M'))
        mrnm[bad.chr] <- reads$seqnames[bad.chr]
    }

    if (inherits(reads, 'GappedAlignments')){
        mwidth = qwidth(reads)
    } else{
        mwidth = reads$qwidth
        mwidth[is.na(mwidth)] = 0
    }

    mwidth[is.na(mpos)] = 0
    mwidth[bad.chr] = 0;  # we set mates with "bad" chromosomes to have 0 width
    mpos[is.na(mpos)] = 1;

    if (inherits(reads, 'GappedAlignments')){
        GRanges(mrnm, IRanges(mpos, width = mwidth), strand = c('+', '-')[1+bamflag(reads)[, 'isMateMinusStrand']], seqlengths = seqlengths(reads), qname = values(reads)$qname, mapq = mapq)
    } else if (inherits(reads, 'GRanges')){
        GRanges(mrnm, IRanges(mpos, width = mwidth), strand = c('+', '-')[1+bamflag(reads$flag)[, 'isMateMinusStrand']], seqlengths = seqlengths(reads), qname = values(reads)$qname, mapq = mapq)
    } else if (inherits(reads, 'data.table')){
        ab=data.table(seqnames=mrnm, start=mpos, end=mpos + mwidth - 1, strand=c('+','-')[1+bamflag(reads$flag)[,'isMateMinusStrand']], qname=reads$qname, mapq = mapq)
    }
}





#' @name get.pairs.grl
#' @title Get coverage as GRanges from BAM on custom set of GRanges
#' @description
#'
#' Takes reads object and returns GRangesList with each read and its mate (if exists)
#'
#' @param reads GRanges holding reads
#' @param pairs.grl.split boolean Flag to return as GRangesList if TRUE (default = TRUE)
#' @param verbose boolean verbose flag (default = FALSE)
#' @export
get.pairs.grl = function(reads, pairs.grl.split = TRUE, verbose = FALSE)
{
    isdt = inherits(reads, 'data.table')

    bad.col = c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", 
                "isCircular", "genome", "start", "end", "width", "element")

    if (verbose){
        cat('deduping\n')
    }

    if (is(reads, 'GappedAlignmentPairs') | is(reads, 'GAlignmentPairs')){
        reads = unlist(reads)
    }

    if (inherits(reads, 'GRanges')) {
        d = duplicated(values(reads)$qname)  ## duplicates are already paired up
        qpair.ix = values(reads)$qname %in% unique(values(reads)$qname[d])
    } else if (isdt){
        d = duplicated(reads$qname)
        qpair.ix = reads$qname %in% unique(reads$qname[d])
    }

    if (!inherits(reads, 'GenomicRanges') && !inherits(reads, 'data.table'))
    {
        if (verbose){
            cat('converting to GRanges\n')
        }
        r.gr = granges(reads)
    } else if (!isdt){
        r.gr = reads[, c()]
    } else{
        r.gr = reads
    }

    if (verbose){
        cat('grbinding\n')
    }

    m.gr = get.mate.gr(reads[!qpair.ix]);

    if (inherits(reads, 'GRanges')) {
        m.val = values(m.gr)
        values(m.gr) = NULL;
        r.gr = c(r.gr, m.gr);
        mcols(r.gr) = rrbind(mcols(reads)[, setdiff(colnames(values(reads)), bad.col), drop = FALSE], m.val)
    } else if (isdt) {
        m.gr = m.gr[, setdiff(colnames(reads), colnames(m.gr)) := NA, with = FALSE]
        r.gr = rbind(reads, m.gr, use.names = TRUE)
        setkey(r.gr, qname)
    }

    if (pairs.grl.split && !isdt) {
        if (verbose){
            cat('splitting\n')
        }
        return(split(r.gr, as.character(r.gr$qname)))
    } else {
        return(r.gr)
    }
}





#' @name count.clips
#' @title Return data.frame with fields of "right" soft clips and "left" soft clips
#' @description
#'
#' Takes GRanges or GappedAlignments object and uses cigar field (or takes character vector of cigar strings)
#' and returns data.frame with fields (for character input)
#' $right.clips number of "right" soft clips (e.g. cigar 89M12S)
#' $left.clips number of "left" soft clips (e.g. cigar 12S89M), 
#' or appends these fields to the reads object
#'
#' @param reads GenomicRanges or GappedAlignments or data.frame or data.table holding the reads
#' @return GRanges with 'right.clips' and 'left.clips' columns added
#' @export
count.clips = function(reads)
{
    if (!inherits(reads, 'GRanges') & !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame') & !inherits(reads, 'data.table')){
        stop('Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.')
    }

    if (length(reads) == 0){
        return(reads)
    }
    if (inherits(reads, 'GRanges') | inherits(reads, 'GappedAlignments')){
        cigar = values(reads)$cigar
    } else{
        cigar = reads
    }

    if (!inherits(cigar, 'character') & !inherits(cigar, 'factor')){
        stop('Input must be GRanges, GappedAlignments, or character vector. Please see documentation for details.')
    }

    out = data.frame(left.clips = rep(0, length(cigar)), right.clips = rep(0, length(cigar)));

    re.left = '^(\\d+)S.*'
    re.right = '.*[A-Z](\\d+)S$'

    lclip.ix = grep(re.left, cigar)
    rclip.ix = grep(re.right, cigar)

    if (length(lclip.ix) > 0){
        left.clips = gsub(re.left, '\\1', cigar[lclip.ix])
        out$left.clips[lclip.ix] = as.numeric(left.clips)
    }

    if (length(rclip.ix) > 0){
        right.clips = gsub(re.right, '\\1', cigar[rclip.ix])
        out$right.clips[rclip.ix] = as.numeric(right.clips)
    }

    if (inherits(reads, 'GRanges') | inherits(reads, 'GappedAlignments')){
        values(reads)$right.clips = out$right.clips
        values(reads)$left.clips = out$left.clips
        out = reads
    }

    return(out)
}





## alpha() used in varbase
alpha = function(col, alpha)
{
    col.rgb = col2rgb(col)
    out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
    names(out) = names(col)
    return(out)
}





#' @name varbase
#' @title Returns variant bases and ranges from GRanges or GappedAlignments input
#' @description
#'
#' Takes GRanges or GappedAlignments object "reads" and uses cigar, MD, seq fields
#' to return variant bases and ranges
#'
#' Returns GRangesList (of same length as input) of variant base positions with character vector 
#' $varbase field populated with variant bases for each GRanges item in grl[[k]], 
#' with the following handling for insertions, deletions, and substitution GRange's:
#'
#' Substitutions: nchar(gr$varbase) = width(gr) of the corresponding var
#' Insertions: nchar(gr$varbase)>=1, width(gr) ==0
#' Deletions: gr$varbase = '', width(gr)>=1
#'
#' Each GRanges also has $type flag which shows the cigar string code for the event i.e.
#' S = soft clip --> varbase represents clipped bases
#' I = insertion --> varbase represents inserted bases
#' D = deletion --> varbase is empty
#' X = mismatch --> varbase represents mismatched bases
#'
#' @param reads GenomicRanges or GRangesList or GappedAlignments or data.frame/data.table reads to extract variants from
#' @param soft boolean Flag to include soft-clipped matches (default = TRUE)
#' @param verbose boolean verbose flag (default = TRUE)
#' @name varbase
#' @export
varbase = function(reads, soft = TRUE, verbose = TRUE)
{
    nreads = length(reads)
    if (inherits(reads, 'GRangesList')){
        was.grl = TRUE
        r.id = grl.unlist(reads)$grl.ix
        reads = unlist(reads)
    } else if (inherits(reads, 'data.frame')){
        r.id = 1:nrow(reads)
        nreads = nrow(reads)
        was.grl = FALSE
    } else{
        r.id = 1:length(reads)
        was.grl = FALSE
    }

    if (!inherits(reads, 'GRanges') & !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame')  & !inherits(reads, 'data.table')){
        stop('Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.')
    }

    if (is.data.frame(reads)){

        sl = NULL
        sn =  reads$seqnames
        cigar = as.character(reads$cigar)
        seq = as.character(reads$seq)
        str = reads$strand

        if (!is.null(reads$MD)){
            md = as.character(reads$MD)
        } else{
            md = rep(NA, length(cigar))
        }
    } else{

        sl = seqlengths(reads)
        sn =  seqnames(reads)
        cigar = as.character(values(reads)$cigar)
        seq = as.character(values(reads)$seq)
        str = as.character(strand(reads))

        if (!is.null(values(reads)$MD)){
            md = as.character(values(reads)$MD)
        } else{
            md = rep(NA, length(cigar))
        }
    }

    if (!inherits(cigar, 'character')){
        stop('Error: Input must be GRanges with seq, cigar, and MD fields populated or GappedAlignments object. Please see documentation for details.')
    }

    ix = which(!is.na(cigar))

    if (length(ix)==0){
        return(rep(GRangesList(GRanges()), nreads))
    }

    cigar = cigar[ix]
    seq = seq[ix]
    md = md[ix]
    str = str[ix]

    if (is.data.frame(reads)){
        r.start = reads$start[ix]
        r.end = reads$end[ix]
    } else{
        r.start = start(reads)[ix]
        r.end = end(reads)[ix];
    }

    flip = str == '-'

    if (!is.null(seq)){
        nix = sapply(seq, function(x) all(is.na(x)))
        if (any(nix)){
            seq[nix] = ''
        }
        seq = strsplit(seq, '')
    }

    cigar.vals = explodeCigarOps(cigar)
    cigar.lens = explodeCigarOpLengths(cigar)

    clip.left = sapply(cigar.vals, function(x) x[1] %in% c('H', 'S'))
    clip.right = sapply(cigar.vals, function(x) x[length(x)] %in% c('H', 'S'))

    if (any(clip.left)){
        r.start[clip.left] = r.start[clip.left]-sapply(cigar.lens[which(clip.left)], function(x) x[1])
    }

    if (any(clip.right)){
        r.end[clip.right] = r.end[clip.right]+sapply(cigar.lens[which(clip.right)], function(x) x[length(x)])
    }

    md.vals = strsplit(gsub('([a-zA-Z])', '|\\1|', gsub('\\^[a-zA-Z]+', '|', md)), '\\|')


    ## ranges of different cigar elements relative to query ie read-centric coordinates
    starts.seq = lapply(1:length(cigar.lens), function(i){

        x = c(0, cigar.lens[[i]])
        x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))+1] = 0  ## deletions have 0 width on query
        cumsum(x[1:(length(x)-1)])+1
    })

    ends.seq = lapply(1:length(cigar.lens), function(i){

        x = cigar.lens[[i]];
        x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))] = 0
        cumsum(x)
    })

    ## ranges of different cigar elements relative to reference coordinatse
    starts.ref = lapply(1:length(cigar.lens), function(i){

        x = c(0, cigar.lens[[i]]);
        x[which(cigar.vals[[i]] %in% c('I'))+1] = 0 ## insertions have 0 width on reference / subject
        cumsum(x[1:(length(x)-1)]) + r.start[i]
    })

    ends.ref = lapply(1:length(cigar.lens), function(i){

        x = cigar.lens[[i]];
        x[which(cigar.vals[[i]] %in% c('I'))] = 0
        cumsum(x) + r.start[i] - 1
    })

    ## now using MD find coordinates of mismatched bases (using starts and ends of M regions)
    ## find coords of subs on genome
    tmp.pos = lapply(1:length(md.vals), function(i){

        x = md.vals[[i]]

        nix = grepl('[a-zA-Z]', x);

        if (!any(nix)){
            return(c())
        }

        p = rep(0, length(x))
        p[!nix] = as.numeric(x[!nix])
        p[nix] = 1
        s.pos.m = cumsum(p)[nix] ## position of subs in read
        mix = cigar.vals[[i]]=='M'
        m.st = cumsum(c(1, ends.seq[[i]][mix]-starts.seq[[i]][mix]+1))
        m.st.g = starts.ref[[i]][mix]
        m.st.r = starts.seq[[i]][mix]
        s.match = rep(NA, length(s.pos.m)) ## indices of matching "M" cigar element for each sub
        
        for (ii in 1:length(s.pos.m)){
            j = 0;
            done = FALSE
            for (j in 0:(length(m.st)-1)){
                if (s.pos.m[ii] < m.st[j+1]){
                    break
                }
            }
            s.match[ii] = j
        }

        s.pos.g = m.st.g[s.match] + s.pos.m-m.st[s.match]
        s.pos.r = m.st.r[s.match] + s.pos.m-m.st[s.match]

        return(rbind(s.pos.g, s.pos.r))
    })
    subs.pos = lapply(tmp.pos, function(x) x[1,])
    subs.rpos = lapply(tmp.pos, function(x) x[2,])
    if (is.null(seq)){
        subs.base = lapply(lapply(md.vals, grep, pattern = '[ATGCNatcgn]', value = T), function(x) rep('X', length(x))) ## replace with N
    } else{
        subs.base = lapply(1:length(seq), function(x) ifelse(is.na(seq[[x]][subs.rpos[[x]]]), 'X', seq[[x]][subs.rpos[[x]]]))
    }      

    ## make sure MD and cigar are consistent
    ## (for some reason - sometimes soft clipped mismatches are included in MD leading to a longer MD string)
    ## also some MD are NA
    mlen.cigar = sapply(1:length(ends.seq), function(x) {mix = cigar.vals[[x]]=='M'; sum(ends.seq[[x]][mix]-starts.seq[[x]][mix]+1)})
                                        #    mlen.md = sapply(md.vals, function(x) {ix = grepl('[ATGCN]', x); sum(as.numeric(x[!ix])) + sum(nchar(x[ix]))})
    mlen.md = sapply(md.vals, function(x) {ix = grepl('[a-zA-Z]', x); sum(as.numeric(x[!ix])) + sum(nchar(x[ix]))})
    good.md = which(!is.na(md))
    ##  good.md = which(mlen.md == mlen.cigar & !is.na(md))

    if (any(na = is.na(md))){

        warning('Warning: MD field absent from one or more input reads')
        good.md = which(!na)
    }
    ## now make GRanges of subs
    if (length(good.md) > 0){

        if (any(mlen.md[good.md] != mlen.cigar[good.md])){
            warning('Warning: The lengths of some MD strings do not match the number of M positions on the corresponding CIGAR string: some variants may not be correctly mapped to the genome')
        }

        iix.md = unlist(lapply(good.md, function(x) rep(x, length(subs.pos[[x]]))))
        tmp = unlist(subs.pos[good.md])

        if (!is.null(tmp)){
            subs.gr = GRanges(sn[ix][iix.md], IRanges(tmp, tmp), strand = '*')
            values(subs.gr)$varbase = unlist(subs.base[good.md])
            values(subs.gr)$varlen = 1
            values(subs.gr)$type = 'X'
            values(subs.gr)$iix = ix[iix.md]
        } else{
            subs.gr = GRanges()
        }
    } else{
        subs.gr = GRanges()
    }

    iix = unlist(lapply(1:length(cigar.vals), function(x) rep(x, length(cigar.vals[[x]]))))
    cigar.vals = unlist(cigar.vals)
    cigar.lens = unlist(cigar.lens)
    starts.seq = unlist(starts.seq)
    ends.seq = unlist(ends.seq)
    starts.ref = unlist(starts.ref)
    ends.ref = unlist(ends.ref)

    ## pick up other variants (including soft clipped and indel)
    is.var = cigar.vals != 'M'
    iix = iix[is.var]
    cigar.vals = cigar.vals[is.var]
    cigar.lens = cigar.lens[is.var]
    starts.ref = starts.ref[is.var]
    ends.ref = ends.ref[is.var]
    starts.seq = starts.seq[is.var]
    ends.seq = ends.seq[is.var]
    str <- str[iix] 

    if (length(cigar.vals)>0){
        
        var.seq = lapply(1:length(cigar.vals), function(i){                                 
            
            if (ends.seq[i]<starts.seq[i]){
                return('') # deletion
            } else{
                                
                if (length(seq[[iix[i]]])==0){
                    rep('N', ends.seq[i]-starts.seq[i]+1)
                } else{
                    seq[[iix[i]]][starts.seq[i]:ends.seq[i]] #insertion
                }
            }
        })

        other.gr = GRanges(sn[ix][iix], IRanges(starts.ref, ends.ref), strand = str, seqlengths = sl)        
        values(other.gr)$varbase = sapply(var.seq, paste, collapse = '')
        values(other.gr)$varlen = cigar.lens
        values(other.gr)$type = cigar.vals
        values(other.gr)$iix = ix[iix];
        out.gr = sort(c(subs.gr, other.gr))
    } else{
        out.gr = subs.gr
    }
    ## add default colors to out.gr
    VAR.COL = c('XA' = 'green', 'XG' = 'brown', 'XC' = 'blue', 'XT' = 'red', 'D' = 'white',  'I'= 'purple', 'N' = alpha('gray', 0.2), 'XX' = 'black', 'S' = alpha('pink', 0.9))

    col.sig = as.character(out.gr$type)
    xix = out.gr$type == 'X'
    col.sig[xix] = paste(col.sig[xix], out.gr$varbase[xix], sep = '')
    out.gr$col = VAR.COL[col.sig]
    out.gr$border = out.gr$col

    if (!soft){
        if (any(soft.ix <<- out.gr$type %in% c('H', 'S'))){
            out.gr = out.gr[-which(soft.ix)]
        }
    }

    out.grl = rep(GRangesList(GRanges()), nreads)
    out.iix = r.id[values(out.gr)$iix]
    values(out.gr)$iix = NULL

    tmp.grl = GenomicRanges::split(out.gr, out.iix)
    out.grl[as.numeric(names(tmp.grl))] = tmp.grl
    values(out.grl)$qname[r.id] = reads$qname

    return(out.grl)
}




#' @name splice.cigar
#' @title Get coverage as GRanges from BAM on custom set of GRanges
#' @description
#'
#' Takes GRanges or GappedAlignments object "reads" and parses cigar fields
#' to return GRanges or GRangesList corresponding to spliced alignments on the genome, which 
#' correspond to portions of the cigar
#'
#' i.e. each outputted GRanges/GRangesList element contains the granges corresponding to all non-N portions of cigar string
#'
#' If GRangesList provided as input (e.g. paired reads) then all of the spliced ranges resulting from each
#' input GRangesList element will be put into the corresponding output GRangesList element
#'
#' NOTE: does not update MD tag
#'
#' If use.D = TRUE, then will treat "D" flags (deletion) in addition to "N" flags as indicative of deletion event.
#'
#' @param reads GenomicRanges or GappedAlignments or data.frame input reads
#' @param verbose boolean verbose flag (default = TRUE)
#' @param fast boolean Flag to use 'GenomicAlignments::cigarRangesAlongReferenceSpace()' to translate CIGAR to GRanges (default = TRUE)
#' @param use.D boolean Treats "D" tags as deletions, along with "N" tags (default = TRUE)
#' @param rem.soft boolean Pick up splice 'S', soft-clipped (default = TRUE)
#' @param get.seq boolean Get InDels (default = TRUE)
#' @param return.grl boolean Return as GRangesList (default = TRUE)
#' @export
splice.cigar = function(reads, verbose = TRUE, fast = TRUE, use.D = TRUE, rem.soft = TRUE, get.seq = FALSE, return.grl = TRUE)
{
    if (!inherits(reads, 'GRanges') & !inherits(reads, 'GRangesList')& !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame') & !inherits(reads, 'data.table')){
        stop('Error: Reads must be either GRanges, GRangesList, GappedAlignments, or data.table object. Please see documentation for details.')
    }

    if (inherits(reads, 'GRangesList')){
      reads = unlist(reads)
    }

    nreads = length(reads)

    if (nreads==0){
        if (return.grl){
            return(GRangesList())
        } else{
            return(GRanges)
        }
    }

    if (inherits(reads, 'GRangesList')){
        was.grl = TRUE
        r.id = as.data.frame(reads)$element
        reads = unlist(reads)
    } else{
        r.id = 1:length(reads)
        was.grl = FALSE
    }

    if (is.data.frame(reads)){
        sl = NULL
        sn =  reads$seqnames
        cigar = as.character(reads$cigar)
        seq = as.character(reads$seq)
        str = reads$strand

        if (!is.null(reads$MD)){
            md = as.character(reads$MD)
        } else{
            md = rep(NA, length(cigar))
        }
    } else{
        sl = seqlengths(reads)
        sn =  seqnames(reads)
        cigar = as.character(values(reads)$cigar)
        seq = as.character(values(reads)$seq)
        str = as.character(strand(reads))

        if (!is.null(values(reads)$MD)){
            md = as.character(values(reads)$MD)
        } else{
            md = rep(NA, length(cigar))
        }
    }

    if (is.null(values(reads)$cigar) | is.null(values(reads)$seq)){
        stop('Error: Reads must have cigar and seq fields specified. Please see documentation for details.')
    }

    if (!inherits(cigar, 'character')){
        stop('Error: Input must be GRanges with seq, cigar, and MD fields populated or GappedAlignments object. Please see documentation for details.')
    }


    ix = which(!is.na(reads$cigar))

    if (length(ix)==0){
        if (return.grl){
            return(rep(GRangesList(GRanges()), nreads))
        } else{
            return(GRanges())
        }       
    }

    if (fast){
        ir = cigarRangesAlongReferenceSpace(reads[ix]$cigar, N.regions.removed = FALSE, with.ops = TRUE, reduce.ranges = FALSE)
        irul = unlist(ir)
        out.gr = GRanges(rep(seqnames(reads)[ix], elementNROWS(ir)), IRanges::shift(IRanges(irul), rep(start(reads)[ix]-1, elementNROWS(ir))),
                         strand = rep(strand(reads)[ix], elementNROWS(ir)), seqlengths = seqlengths(reads))
        out.gr$type = names(irul)
        out.gr$rid = ix[rep(1:length(ir), elementNROWS(ir))]
        out.gr$riid = unlist(lapply(elementNROWS(ir), function(x) 1:x))
        out.gr$fid = r.id[out.gr$rid]
        out.gr$qname = reads$qname[out.gr$rid]

        if (return.grl){
            out.grl = rep(GRangesList(GRanges()), nreads)
            tmp.grl = split(out.gr, out.gr$fid)
            out.grl[as.numeric(names(tmp.grl))] = tmp.grl
            return(out.grl)
        } else{
            return(out.gr)
        }
    } else{

        cigar = cigar[ix]
        str = str[ix]

        if (is.data.frame(reads)){
            r.start = reads$start[ix]
            r.end = reads$end[ix]
        } else{
            r.start = start(reads)[ix]
            r.end = end(reads)[ix];
        }

        flip = str == '-'
        cigar.vals = lapply(strsplit(cigar, '\\d+'), function(x) x[2:length(x)])
        cigar.lens = lapply(strsplit(cigar, '[A-Z]'), as.numeric)

        clip.left = sapply(cigar.vals, function(x) x[1] == 'S')
        clip.right = sapply(cigar.vals, function(x) x[length(x)] == 'S')

                                        # ranges of different cigar elements relative to query ie read-centric coordinates
        starts.seq = lapply(1:length(cigar.lens), function(i){
            x = c(0, cigar.lens[[i]])
            x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))+1] = 0  ## deletions have 0 width on query
            cumsum(x[1:(length(x)-1)])+1
        })

        ends.seq = lapply(1:length(cigar.lens), function(i){
            x = cigar.lens[[i]];
            x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))] = 0
            cumsum(x)
        })

        ## ranges of different cigar elements relative to reference coordinatse
        starts.ref = lapply(1:length(cigar.lens), function(i){
            x = c(0, cigar.lens[[i]]);
            x[which(cigar.vals[[i]] %in% c('I'))+1] = 0 ## insertions have 0 width on reference / subject
            cumsum(x[1:(length(x)-1)]) + r.start[i]
        })

        ends.ref = lapply(1:length(cigar.lens), function(i){
            x = cigar.lens[[i]];
            x[which(cigar.vals[[i]] %in% c('I'))] = 0
            cumsum(x) + r.start[i] - 1
        })

        iix = unlist(lapply(1:length(cigar.vals), function(x) rep(x, length(cigar.vals[[x]]))))
        cigar.vals = unlist(cigar.vals)
        cigar.lens = unlist(cigar.lens)
        starts.seq = unlist(starts.seq)
        ends.seq = unlist(ends.seq)
        starts.ref = unlist(starts.ref)
        ends.ref = unlist(ends.ref)

        ## pick up splice (including soft clipped and indel)
        splice.char = 'N'

        if (use.D){
            splice.char = c(splice.char, 'D')
        }

        if (rem.soft){
            splice.char = c(splice.char, 'S')
        }

        is.splice = !(cigar.vals %in% splice.char)
        iix = iix[is.splice]
        cigar.vals = cigar.vals[is.splice]
        cigar.lens = cigar.lens[is.splice]
        starts.ref = starts.ref[is.splice]
        ends.ref = ends.ref[is.splice]
        starts.seq = starts.seq[is.splice]
        ends.seq = ends.seq[is.splice]
        str <- str[iix] #

        out.gr = GRanges()
        if (length(cigar.vals)>0){
            other.gr = GRanges(sn[ix][iix], IRanges(starts.ref, ends.ref), strand = str, seqlengths = sl)

            if (get.seq){
                var.seq = lapply(1:length(cigar.vals),
                    function(i){
                        if (ends.seq[i]<starts.seq[i]){
                            return('') ## deletion
                            } else{
                            seq[[iix[i]]][starts.seq[i]:ends.seq[i]] ## insertion
                        }
                    })
                values(other.gr)$seq = sapply(var.seq, paste, collapse = '')
            }

            values(other.gr)$type = cigar.vals
            values(other.gr)$iix = ix[iix];

            out.gr = c(other.gr, out.gr)
        }

        out.gr$rid = out.gr$iix
        out.iix = r.id[out.gr$rid]
        values(out.gr)$iix = NULL
        out.gr$fid = out.iix
        out.gr$qname = reads$qname[out.gr$rid]

        if (return.grl){
            out.grl = rep(GRangesList(GRanges()), nreads)
            tmp.grl = split(out.gr, out.iix)
            out.grl[as.numeric(names(tmp.grl))] = tmp.grl
            return(out.grl)
        } else{
            return(out.gr)
        }
    }
}




#' @name bamflag
#' @title Returns matrix of bits from BAM flags
#' @description
#'
#' Shortcut function: assumes reads are GappedAlignments with flag variable or actual integers representing BAM flag
#'
#' @param reads GenomicRanges or 'GappedAlignments' or data.table holding the reads
#' @return matrix of bits from BAM flags
#' @export
bamflag = function(reads)
{
    if (inherits(reads, 'GappedAlignments') | inherits(reads, 'data.frame') | inherits(reads, 'GRanges')){
        bf = reads$flag
    } else{
        bf = reads
    }

    out = matrix(as.numeric(intToBits(bf)), byrow = T, ncol = 32)[, 1:12, drop = FALSE]
    colnames(out) = c('isPaired', 'isProperPair', 'isUnmappedQuery', 'hasUnmappedMate', 'isMinusStrand', 'isMateMinusStrand', 'isFirstMateRead', 'isSecondMateRead', 'isNotPrimaryRead', 'isNotPassingQualityControls', 'isDuplicate', 'isSupplementary')

    return(out)
}




#' @name bamtag
#' @title Outputs a tag to identify duplicate reads in GRanges input
#' @description
#'
#' Outputs a tag that cats 'qname', first vs first second mate +/- secondary alignment +/- gr.string
#' to give an identifier for determine duplicates in a read pile
#'
#' @param reads GenomicRanges or GappedAlignments or data.frame holding the reads
#' @param secondary boolean including secondary alignment(s) (default = FALSE)
#' @param gr.string boolean input reads into gr.string() (default = FALSE)
#' @export
bamtag = function(reads, secondary = FALSE, gr.string = FALSE)
{
    if (!inherits(reads, 'GRanges') & !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame') & !inherits(reads, 'data.table')){
        stop('Error: Reads must be either GRanges, GRangesList, or GappedAlignments object. Please see documentation for details.')
    }

    grs = sec = NULL
    if (secondary){
        sec = bamflag(reads$flag)[, 'isNotPrimaryRead']
    }

    if (gr.string){
        grs = gr.string(reads, mb = FALSE)
    }

    return(paste(reads$qname, ifelse(bamflag(reads$flag)[, 'isFirstMateRead'], '1', '2'), grs, sec, sep = '_'))
}




#' @name countCigar
#' @title Count bases in cigar string
#' @description
#'
#' Counts the total number of bases, per cigar, that fall into D, I, M, S categories.
#' countCigar makes no distinction between, for instance 1S2M2S, 2S2M1S, or 3S2M
#'
#' @param cigar character vector of cigar strings
#' @return matrix of dimensions (4-column, length(cigar)) with the total counts for each type
#' @export
countCigar = function(cigar){

    cigar.vals = unlist(strsplit(cigar, "\\d+"))
    cigar.lens = strsplit(cigar, "[A-Z]")
    lens = nchar(gsub('\\d+', '', cigar))
    lens[is.na(cigar)] = 1

    cigar.lens = as.numeric(unlist(cigar.lens))
    cigar.vals = cigar.vals[cigar.vals != ""]
    repr       = rep(seq_along(cigar), lens)
    dt         = data.table(val=cigar.vals, lens=cigar.lens, group=repr, key="val")
    
    smr.d = dt["D",][, sum(lens), by=group]
    smr.i = dt["I",][, sum(lens), by=group]
    smr.m = dt["M",][, sum(lens), by=group]
    smr.s = dt["S",][, sum(lens), by=group]

    out = matrix(nrow=length(cigar), ncol=4, 0)
    out[smr.d$group, 1] = smr.d$V1
    out[smr.i$group, 2] = smr.i$V1
    out[smr.m$group, 3] = smr.m$V1
    out[smr.s$group, 4] = smr.s$V1
    colnames(out) = c('D','I','M','S')

    return(out)
}




#' @name is.paired.end
#' @title Check if BAM file is paired end by using 0x1 flag
#' @description
#'
#' Check if BAM file is paired-end by using 0x1 flag, 
#' pipes to 'samtools' via command line
#'
#' @param bams vector of input BAMs
#' @return boolean returns TRUE if BAM file is paired-end, returns FALSE if BAM not paired-end
#' @export
is.paired.end = function(bams)
{
    out = sapply(bams, function(x){
        if (is.na(x)){
            return(NA)
        }
        if (!file.exists(x)){
            return(NA)
        }
        out = FALSE                
        p = pipe(sprintf('samtools view -h  %s | head -n 100 | samtools view -f 0x1 - | wc -l', x))
        ln = as.numeric(readLines(p))   
        out = ln > 0
        close(p)
        return(out)                
    })   

    return(out)
}




#' @name chunk
#' @title  chunk
#' @description
#'
#' Internal function takes same input as seq (from, to, by, length.out) and outputs a 2 column matrix of indices
#' corresponding to "chunks"
#'
#' @param from integer Where to begin sequence
#' @param to integer To end sequence (default = NULL)
#' @param by integer Interval to space sequence (default = 1)
#' @param length.out integer Number of desired chunks, i.e. nrows of output matrix (default = NULL)
#' @return 2-column matrix of indices, each row representing a chunk
#' @author Marcin Imielinski
chunk = function(from, to = NULL, by = 1, length.out = NULL)
{
    if (is.null(to)){
        to = from;
        from = 1;
    }

    if (is.null(length.out)){
        tmp = c(seq(from = from, to = to, by = by), to + 1)
    } else{
        tmp = c(seq(from = from, to = to, length.out = length.out), to + 1)
    }

    out = floor(cbind(tmp[-length(tmp)], tmp[-1]-1))

    return(out)
}




#' @name varcount
#' @title Wrapper around applyPileups
#' @description 
#'
#' Takes in vector of bam paths or GRanges corresponding to sites / territories to query,
#' and outputs a list with fields:
#'
#' $counts = 3D matrix of base counts
#' (A, C, G, T, N) x sites x bams subject to mapq and baseq thresholds
#'
#  $gr = output ranges corresponding to "sites" columns of output
#'
#'
#' varcount() relies upon varbase() 
#'
#' @param bams character vector of paths to bam files
#' @param gr GRanges of (width=1) sites i.e. intervals at which to compute base coujnts
#' @param min.mapq integer Minimal mapping quality at which to compute bases (default = 0)
#' @param min.baseq integer Minimal base quality at which to compute bases (default = 20)
#' @param max.depth integer Maximum read depth to consider (default = 500)
#' @param indel boolean Flag whether to consider indels (default = FALSE)
#' @param ... other args be passed to read.bam(). Please see documentation for read.bam()
#' @return GRanges annotated with fields $alt.count.t, $ref.count.t, $alt.count.n, $ref.count.n
#' @author Marcin Imielinski
#' @export
varcount = function(bams, gr, min.mapq = 0, min.baseq = 20, max.depth = 500, indel = FALSE, ...)
{
    require(abind)
    require(Rsamtools)

    out = list()

    if (any(width(gr)!=1)){
        gr = gr.start(gr)
    }

    
    if (is.character(bams)){    

        bami = gsub('\\.bam$', '.bai', bams)
        ix = file.exists(bami)
        if (any(!ix)){
            bami[!ix] = paste(bams[!ix], 'bai', sep = '.')
        }
        if (any(!file.exists(bami))){
            stop('Error: one or more BAM file indices missing')
        }
        bams = BamFileList(mapply(function(bam, bai) BamFile(bam, index = bai), bams, bami, SIMPLIFY = FALSE))
    } else if (is(bams, 'BamFile')){
        bams = BamFileList(bams)
    }

    sls = lapply(bams, seqlengths)
    sl0 = sort(sls[[1]])
    if (length(sls)>1)
    {
        for (i in 2:length(sls))
        {
            if (!identical(sl0, sort(sls[[i]])))
            {
                stop('Seqlengths are mismatched between two or more bam files in provided list.  Please check your bam files and ensure they have been defined on identical genomes')
            }
        }
    }
    

    ix = as.logical(as.character(seqnames(gr)) %in% seqlevels(bams))

    if (any(ix)){

        pp = ApplyPileupsParam(which = gr[ix], what = c("seq"), minBaseQuality = min.baseq, minMapQuality = min.mapq, maxDepth = max.depth)
        ## ## xtYao fix: function applyPileups fail at heterogeneous BAM seqlevels
        ## ## do them separately and put back
        ## if (length(bams)==2){
        ##     if (identical(seqlengths(bams[1]), seqlengths(bams[2]))){
        ##         pu = applyPileups(PileupFiles(bams), function(x) x, param = pp)
        ##     } else {
        ##         pu1 = applyPileups(PileupFiles(bams[1]), function(x) x, param = pp)
        ##         pu2 = applyPileups(PileupFiles(bams[2]), function(x) x, param = pp)
        ##         pu = lapply(which(ix),
        ##                     function(ii){
        ##                         out = list()
        ##                         out$seqnames = pu1[[ii]]$seqnames
        ##                         out$pos = pu1[[ii]]$pos
        ##                         if (all(dim(pu1[[ii]]$seq)==dim(pu2[[ii]]$seq))){
        ##                             out$seq = abind(pu1[[ii]]$seq, pu2[[ii]]$seq, along=2)
        ##                             return(out)
        ##                         } else {
        ##                             return(NULL)
        ##                         }
        ##                     })
        ##     }
        ## } else {
        pu = applyPileups(PileupFiles(bams), function(x) x, param = pp)
        ## }
    }

    if (is(bams, 'BamFile') | is(bams, 'BamFileList')){
        bam.paths = Rsamtools::path(bams)
    } else if (is(bams, 'BamFileList')){
        bam.paths = sapply(bams, path)
    } else if (is(bams, 'list')){
        bam.paths = sapply(bams, path)
    } else if (is(bams, 'character')){
        bam.paths = bams
    }


    if (!indel){
        cnames = c('A', 'C', 'G', 'T', 'N')
        out$counts = array(NA, dim = c(length(cnames), length(gr), length(bams)), dimnames = list(cnames, NULL, bam.paths))
        if (any(ix)){
            nna = sapply(pu, function(x) length(x$seq)>0)
            out$counts[,which(ix)[nna],] = aperm(do.call('abind', lapply(pu, function(x){
                x$seq[cnames,,, drop = F]
            })), c(1,3,2))
        }
    } else{
        cnames = unique(unlist(lapply(pu, function(x) rownames(x$seq))))
        cnames = cnames[order(nchar(cnames), cnames)]
        out$counts = array(NA, dim = c(length(cnames), length(gr), length(bams)), dimnames = list(cnames, NULL, bam.paths))
        if (any(ix)){
            nna = sapply(pu, function(x) length(x$seq)>0)
            out$counts[,which(ix)[nna],] = aperm(do.call('abind', lapply(pu, function(x){
                out = array(NA, dim = c(length(cnames), dim(x$seq)[2:3]), dimnames = list(cnames));
                out[rownames(x$seq),, ] = x$seq
                })), c(1,3,2))
            return(out)
            }
    }    

    out$gr = gr

    return(out)    
}






#' @name mafcount
#' @title Wrapper around varcount adapted to tumor and normal "paired" bams
#' @description 
#'
#' Returns base counts for reference and alternative allele for an input tumor (and normal bam) and import MAF as a GRanges specifying substitutions
#'
#' maf is a single width GRanges describing variants and field 'ref' (or 'Reference_Allele'), 'alt' (or 'Tum_Seq_Allele1') specifying reference and alt allele.
#' maf is assumed to have width 1 and strand is ignored.  
#'
#' @param tum.bam string path to tumor sample, input to Bamfile()
#' @param norm.bam optional string path to normal sample, input to Bamfile() (optional) (default = NULL)
#' @param maf GRanges of imported MAF (e.g. output of read.delim or dt2gr(fread(MAF)))
#' @param chunk.size integer Number of variants to extract from bam file at each iteration (default = 100)
#' @param verbose logical Flag whether to print verbose output (default = TRUE)
#' @param mc.cores integer Number of cores in mclapply (default = 1)
#' @param ...  additional pparams to pass to varcount
#' @return GRanges of MAF annotated with fields $alt.count.t, $ref.count.t, $alt.count.n, $ref.count.n
#' @author Marcin Imielinski
#' @export
mafcount = function(tum.bam, norm.bam = NULL, maf, chunk.size = 100, verbose = TRUE, mc.cores = 1, ...)
{


    if (is.character(tum.bam)){
        tum.bam = BamFile(tum.bam)
    }

    ## xtYao: fix here rather than `varcount`
    bams = BamFileList(tum.bam)
        
    if (!is.null(norm.bam)){
            
        if (is.character(norm.bam)){
            norm.bam = BamFile(norm.bam)
        }

        ## prevent incompatible BAM headers
        if (identical(seqlengths(bams), seqlengths(norm.bam))){
            bams = c(bams, BamFileList(norm.bam))
        } else{
            bams2 = BamFileList(norm.bam)
        }
    }
    

    chunks = chunk(1, length(maf), chunk.size)

        
    if (is.null(maf$Tumor_Seq_Allele1)){
        maf$Tumor_Seq_Allele1 = maf$alt
    }
    
   
    if (is.null(maf$Tumor_Seq_Allele1)){
        maf$Tumor_Seq_Allele1 = maf$ALT
    }


    if (is.null(maf$Reference_Allele)){
        maf$Reference_Allele = maf$ref
    }
        
    if (is.null(maf$Reference_Allele)){
        maf$Reference_Allele = maf$REF
    }

                        
    if (is.null(maf$Reference_Allele) | is.null(maf$Tumor_Seq_Allele1)){
        stop('Error: Cannot locate variant columns in input GRanges, please check input to make sure it either has standard VCF ALT / REF columns or MAF file columns specifying alt and ref allele')
    }

    if (!all(is.character(maf$Tumor_Seq_Allele1))){
        maf$Tumor_Seq_Allele1 = sapply(maf$Tumor_Seq_Allele1, function(x) as.character(x)[1])
    }

        
    if (!all(is.character(maf$Reference_Allele))){
        maf$Reference_Allele = as.character(maf$Reference_Allele)
    }

            
    maf$alt.count.t = maf$ref.count.t = NA

    if (!is.null(norm.bam)){
        maf$alt.count.n =  maf$ref.count.n = NA
    }

    if (verbose){
        cat('Initialized\n')
    }

    if (is.data.frame(maf)){
        maf = seg2gr(maf)   ## maf now a GRanges
    }

    tmp = do.call('rbind', mclapply(1:nrow(chunks), function(i){

        if (verbose){
            cat('Starting chunk ', chunks[i, 1], ' to ', chunks[i, 2], '\n')
        }
                
        ix = chunks[i,1]:chunks[i,2]

        if (verbose){
            now = Sys.time()
        }
               
        vc = varcount(bams, maf[ix])

        if (exists("bams2")){
            vc2 = varcount(bams2, maf[ix])
            ## vc$counts = abind(vc$count, vc2$count, along=3)
        }
               
        if (verbose){
            print(Sys.time() - now)
        }
               
        tum.count = vc$counts[,,1]

        if (exists("bams2")){
            norm.count = vc2$counts[,,1]
        }

        if (is.null(dim(tum.count))){
            tum.count = cbind(tum.count)
        }
        
        out = cbind(
            tum.count[cbind(match(maf$Tumor_Seq_Allele1[ix], rownames(tum.count)), 1:length(ix))],
            tum.count[cbind(match(maf$Reference_Allele[ix], rownames(tum.count)), 1:length(ix))]
        )

        if (verbose){
            cat('Num rows:', nrow(out), '\n')
        }
                     
        if (!is.null(norm.bam)){

            ## prevent incompatible BAM headers
            if (identical(seqlengths(bams), seqlengths(norm.bam))){
                norm.count = vc$counts[, , 2]                      
            } else{
                norm.count = vc2$counts[, , 1]
            }
            if (is.null(dim(norm.count))){
                norm.count = cbind(norm.count)
            }

            out = cbind(out, 
                norm.count[cbind(match(maf$Tumor_Seq_Allele1[ix], rownames(norm.count)), 1:length(ix))],
                norm.count[cbind(match(maf$Reference_Allele[ix], rownames(norm.count)), 1:length(ix))]
            )
        }
               
    return(out) 

    }, mc.cores = mc.cores))

    ### write check if only NA's
    if (all(is.na(tmp))){
        return(GRanges()) ### return empty GRanges
    } else{
        maf$alt.count.t = tmp[,1]
        maf$ref.count.t = tmp[,2]
        maf$alt.frac.t = maf$alt.count.t / (maf$alt.count.t + maf$ref.count.t)
        maf$ref.frac.t = 1 - maf$alt.frac.t

        if (!is.null(norm.bam)){
            maf$alt.count.n = tmp[,3]
            maf$ref.count.n = tmp[,4]
            maf$alt.frac.n = maf$alt.count.n / (maf$alt.count.n + maf$ref.count.n)
            maf$ref.frac.n = 1 - maf$alt.frac.n
        }
        return(maf)
    }
}
