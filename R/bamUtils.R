#' @import GenomicRanges
#' @import GenomicAlignments


#' @name read.bam
#' @title read.bam
#' @description 
#' 
#' Read BAM file into GRanges or data.table
#'
#' Wrapper around Rsamtools bam scanning functions,
#' by default, returns GRangesList of read pairs for which <at least one> read lies in the supplied interval
#' @param bam Input bam file. Advisable to make "bam" a BamFile instance instead of a plain string, so that the index does not have to be reloaded.
#' @param bami Input bam index file.
#' @param gr GRanges of intervals to retrieve
#' @param intervals GRanges of intervals to retrieve
#' @param stripstrand Flag to ignore strand information on the query intervals. Default TRUE
#' @param what What fields to pull down from BAM. Default \code{scanBamWhat()}
#' @param unpack.flag Add features corresponding to read flags. Default FALSE
#' @param verbose Increase verbosity
#' @param tag Additional tags to pull down from the BAM (e.g. 'R2')
#' @param isPaired See documentation for \code{scanBamFlag}. Default NA
#' @param isProperPair See documentation for \code{scanBamFlag}. Default NA
#' @param isUnmappedQuery See documentation for \code{scanBamFlag}. Default NA
#' @param hasUnmappedMate See documentation for \code{scanBamFlag}. Default NA
#' @param isNotPassingQualityControls See documentation for \code{scanBamFlag}. Default NA
#' @param isDuplicate See documentation for \code{scanBamFlag}. Default FALSE
#' @param isValidVendorRead See documentation for \code{scanBamFlag}. Default TRUE
#' @param as.grl Return reads as GRangesList. Controls whether \code{get.pairs.grl} does split. Default TRUE
#' @param as.data.table Return reads in the form of a data.table rather than GRanges/GRangesList
#' @param ignore.indels messes with cigar to read BAM with indels removed. Useful for breakpoint mapping on contigs
#' @param size.limit Default 1e6
#' @param ... passed to \code{scanBamFlag}
#' @return Reads in one of GRanges, GRangesList or data.table
#' @export
read.bam = function(bam, intervals = NULL,## GRanges of intervals to retrieve
                    gr = intervals,
                    all = FALSE,
                    bami = NULL,
                    pairs.grl = TRUE, # if TRUE will return GRangesList of read pairs for whom at least one read falls in the supplied interval
                                        #  paired = F, # if TRUE, will used read bam gapped alignment pairs warning: will throw out pairs outside of supplied window
                                        #  gappedAlignment = T, # if false just read alignments using scanbam
                    stripstrand = TRUE,
                    what = scanBamWhat(),
                    unpack.flag = FALSE, # will add features corresponding to read flags
                    verbose = FALSE,
                    tag = NULL,
                    isPaired = NA, ## if these features are NA, then reads satisfying both T and F will be returned
                    isProperPair = NA,
                    isUnmappedQuery = NA,
                    hasUnmappedMate = NA,
                    isNotPassingQualityControls = NA,
                    isDuplicate = F,
                    isValidVendorRead = TRUE,
                    as.grl=TRUE, ## return pairs as grl, rather than GRanges .. controls whether get.pairs.grl does split (t/c rename to pairs.grl.split)
                    as.data.table=FALSE, ## returns reads in the form of a data table rather than GRanges/GRangesList
                    ignore.indels=FALSE, ## messes with cigar to read BAM with indels removed. Useful for breakpoint mapping on contigs
                    size.limit = 1e6,
                    ... # passed to scanBamFlag (
                    )
{
    if (!inherits(bam, 'BamFile'))
    {
        if (is.null(bami))
        {
            if (file.exists(bai <- gsub('.bam$', '.bai', bam)))
                bam = BamFile(bam, bai)
            else if (file.exists(bai <- paste(bam, '.bai', sep = '')))
                bam = BamFile(bam, bai)
            else
                bam = BamFile(bam)
        }
        else
            bam = BamFile(bam, index = bami)
    }


                                        # if intervals unspecified will try to pull down entire bam file (CAREFUL)

    if (length(intervals)==0)
        intervals = NULL

    if (is.null(intervals))
        intervals = gr

    if (is.null(intervals))
    {
        if (all)
            intervals = si2gr(seqinfo(bam))
        else
            stop('Must provide non empty interval list')
    }

    if (class(intervals) == 'data.frame')
        intervals = seg2gr(intervals);

    if (inherits(intervals, 'GRangesList'))
        intervals = unlist(intervals);

    if (stripstrand)
        strand(intervals) = '*'

    intervals = reduce(intervals);

    now = Sys.time();

    if (pairs.grl)
        paired = F

    flag = scanBamFlag(isPaired = isPaired, isProperPair = isProperPair, isUnmappedQuery = isUnmappedQuery,
                       hasUnmappedMate = hasUnmappedMate, isNotPassingQualityControls = isNotPassingQualityControls,
                       isDuplicate = isDuplicate, ...)

    tag = unique(c('MD', 'MQ', tag))
    param = ScanBamParam(which = gr.fix(intervals, bam, drop = T), what = what, flag = flag, tag = tag)

    if (verbose)
        cat('Reading bam file\n')
    if (class(bam) == 'BamFile')
        out <- scanBam(bam, param=param)
    else
        out <- scanBam(bam, index=bami, param=param)
    if (verbose) {
        print(Sys.time() - now)
        print('BAM read. Making into data.frame')
    }

    out <- out[sapply(out, function(x) length(x$qname)>0)]

    if (length(out)>0)
    {
        if (verbose) {
            print(Sys.time() - now)
            print('combining lists')
        }
        out <- as.data.frame(rbindlist(lapply(out, function(x)
        {
            x <- c(x[-match('tag', names(x))], x$tag)

            x <- x[sapply(x, length)>0]
            conv <- which(!(sapply(x, class) %in% c('integer', 'numeric', 'character')))
            x[conv] <- lapply(x[conv], as.character)

            for (t in tag)
                if (!(t %in% names(x)))
                    x[[t]] = rep(NA, length(x$qname))

            if (!('R2' %in% names(x)) && 'R2' %in% tag)
                x$R2 <- rep(NA, length(x$qname))
            if (!('Q2' %in% names(x)) && 'Q2' %in% tag)
                x$Q2 <- rep(NA, length(x$qname))
            x
        })))

        ## faster CIGAR string parsing with vectorization and data tables
        if (verbose) {
            print(Sys.time() - now)
            print('filling pos2 from cigar')
        }
        if (ignore.indels) {
            cigar <- gsub('[0-9]+D', '', gsub('([0-9]+)I', '\\1M', out$cigar))  ## Remove deletions, turn insertions to matches
            cig <- splitCigar(cigar)
            torun=sapply(cig, function(y) any(duplicated((y[[1]][y[[1]]==M]))))
            M <- charToRaw('M')
            new.cigar <- sapply(cig[torun], function(y) {
                lets <- y[[1]][!duplicated(y[[1]])]
                vals <- y[[2]][!duplicated(y[[1]])]
                vals[lets==M] <- sum(y[[2]][y[[1]]==M])
                lets <- strsplit(rawToChar(lets), '')[[1]]
                paste(as.vector(t(matrix(c(vals, lets), nrow=length(vals), ncol=length(lets)))), collapse='')
            })
            out$cigar[torun] <- new.cigar
        }
        cigs <- countCigar(out$cigar)
        out$pos2 <- out$pos + cigs[, "M"]

        if (verbose) {
            print(Sys.time() - now)
            print('fixing seqdata')
        }
        out$qwidth = nchar(out$seq)
        unm = is.na(out$pos)
        if (any(unm))
        {
            out$pos[unm] = 1
            out$pos2[unm] = 0
            out$strand[unm] = '*'
        }
        gr.fields = c('rname', 'strand', 'pos', 'pos2');
        vals = out[, setdiff(names(out), gr.fields)]

        if (!as.data.table) {
            out <- GRanges(out$rname, IRanges(out$pos, pmax(0, out$pos2-1)), strand = out$strand, seqlengths = seqlengths(intervals))
            values(out) <- vals;
        } else {
            out <- data.table(seqnames=out$rname, start=out$pos, end= pmax(out$pos2-1, 0), strand=out$strand)
            val <- data.table(vals)
            out <- cbind(out, val)
        }
                                        #out$uname = paste(out$qname, ifelse(bamflag(out$flag)[, 'isFirstMateRead'], '_r1', '_r2'), sep = '')
    }
    else {
        if (!as.data.table)
            return(GRanges(seqlengths = seqlengths(intervals)))
        else
            return(data.table())
    }

    if (verbose)
    {
        if (as.data.table)
            cat(sprintf('Extracted %s reads\n', nrow(out)))
        else
            cat(sprintf('Extracted %s reads\n', length(out)))
        print(Sys.time() - now)
    }

    if (pairs.grl)
    {
        if (verbose)
            cat('Pairing reads\n')
        out <- get.pairs.grl(out, as.grl=as.grl)
        if (verbose)
        {
            cat('done\n')
            print(Sys.time() - now)
        }
        if (as.grl && !as.data.table) {
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
#' gets coverage from bam in supplied ranges using "countBam", returning gr with coverage counts in
#' each of the provided ranges (different from bam.cov above) specified as $file, $records, and $nucleotides
#' columns in the values field
#' basically a wrapper for countBam with some standard settings for ScanBamParams
#'
#' @param bam Input bam file. Advisable to make "bam" a BamFile instance instead of a plain string, so that the index does not have to be reloaded.
#' @param bami Input bam index file.
#' @param gr GRanges of intervals to retrieve
#' @param verbose Increase verbosity
#' @param isPaired See documentation for \code{scanBamFlag}. Default NA
#' @param isProperPair See documentation for \code{scanBamFlag}. Default NA
#' @param isUnmappedQuery See documentation for \code{scanBamFlag}. Default NA
#' @param hasUnmappedMate See documentation for \code{scanBamFlag}. Default NA
#' @param isNotPassingQualityControls See documentation for \code{scanBamFlag}. Default NA
#' @param isDuplicate See documentation for \code{scanBamFlag}. Default FALSE
#' @param isValidVendorRead See documentation for \code{scanBamFlag}. Default TRUE
#' @param mc.cores Number of cores in \code{mclapply} call.
#' @param chunksize How many intervals to process per core. Default 10.
#' @param ... passed to \code{scanBamFlag}
#' @return GRanges parallel to input GRanges, but with metadata filled in.
#' @export
bam.cov.gr = function(bam, gr, bami = NULL, count.all = FALSE, isPaired = T, isProperPair = T, isUnmappedQuery = F, hasUnmappedMate = F, isNotPassingQualityControls = F, isDuplicate = F, isValidVendorRead = T, mc.cores = 1, chunksize = 10, verbose = F, ...)
{
    if (is.character(bam))
        if (!is.null(bami))
            bam = BamFile(bam, bami)
        else
        {
            if (file.exists(paste(bam, 'bai', sep = '.')))
                bam = BamFile(bam, paste(bam, 'bai', sep = '.'))
            else if (file.exists(gsub('.bam$', '.bai', bam)))
                bam = BamFile(bam, paste(bam, 'bai', sep = '.'))
            else
                stop('BAM index not found, please find index and specify bam file argument as valid BamFile object')
        }

    keep = which(as.character(seqnames(gr)) %in% seqlevels(bam))
    
    if (length(keep)>0)
    {
        ix = c(keep[c(seq(1, length(keep), chunksize))], keep[length(keep)]+1);  ## prevent bam error from improper chromosomes
        chunk.id = unlist(lapply(1:(length(ix)-1), function(x) rep(x, ix[x+1]-ix[x])))

        gr.chunk = split(gr[keep], chunk.id[keep]);
        if (count.all)
            flag = scanBamFlag()
        else
            flag = scanBamFlag(isPaired = isPaired, isProperPair = isProperPair, isUnmappedQuery = isUnmappedQuery,
                               hasUnmappedMate = hasUnmappedMate, isNotPassingQualityControls = isNotPassingQualityControls,
                               isDuplicate = isDuplicate, ...)
        out = rbindlist(mclapply(1:length(gr.chunk),
                                 function(x) {
                                     if (verbose)
                                         cat(sprintf('Processing ranges %s to %s of %s, extracting %s bases\n', ix[x], ix[x+1]-1, length(keep), sum(width(gr.chunk[[x]]))))
                                     as.data.table(countBam(bam, param = ScanBamParam(which = gr.chunk[[x]], flag = flag)))
                                 }, mc.cores = mc.cores));

        gr.tag = paste(as.character(seqnames(gr)), start(gr), end(gr));
        out.tag = paste(out$space, out$start, out$end);
        ix = match(gr.tag, out.tag);
        values(gr) = cbind(as.data.frame(values(gr)), out[ix, c('file', 'records', 'nucleotides'), with = FALSE])
    }
    else
        values(gr) = cbind(as.data.frame(values(gr)), data.frame(file = rep(gsub('.*\\/([^\\/]+)$', '\\1', path(bam)), length(gr)), records = NA, nucleotides = NA))

    return(gr)
}

#' @name bam.cov.tile
#' @title Get coverage as GRanges from BAM on genome tiles across seqlengths of genome
#' @description
#' Quick way to get tiled coverage via piping to samtools (~10 CPU-hours for 100bp tiles, 5e8 read pairs)
#'
#' Gets coverage for window size "window", pulling "chunksize" records at a time and incrementing bin
#' corresponding to midpoint or overlaps of corresponding (proper pair) fragment (uses TLEN and POS for positive strand reads that are part of a proper pair)
#'
#' @param bam.file character scalar input bam file
#' @param window integer scalar window size (in bp)
#' @param chunksize integer scalar, size of window
#' @param min.mapq integer scalar, minimim map quality reads to consider for counts
#' @param verbose dummy
#' @param max.tlen max paired-read insert size to consider
#' @param st.flag samtools flag to filter reads on [Default: -f 0x02 -F 0x10]
#' @param fragments dummy
#' @param region dummy
#' @param do.gc dummy
#' @param midpoint if TRUE will only use the fragment midpoint, if FALSE will count all bins that overlap the fragment
#' @return GRanges of "window" bp tiles across seqlengths of bam.file with meta data field $counts specifying fragment counts centered
#' in the given bin.
#' @export
bam.cov.tile = function(bam.file, window = 1e2, chunksize = 1e5, min.mapq = 30, verbose = TRUE,
                        max.tlen = 1e4, ## max insert size to consider
                        st.flag = "-f 0x02 -F 0x10",
                        fragments = TRUE,
                        region = NULL,
                        do.gc = FALSE,
                        midpoint = TRUE ## if TRUE will only use the fragment midpoint, if FALSE will count all bins that overlap the fragment
                        )
{
    cmd = 'samtools view %s %s -q %s | cut -f "3,4,9"' ## cmd line to grab the rname, pos, and tlen columns

    sl = seqlengths(BamFile(bam.file))

    counts = lapply(sl, function(x) rep(0, ceiling(x/window)))
    numwin = sum(sapply(sl, function(x) ceiling(x/window)))

    if (!is.null(region))
    {
        cat(sprintf('Limiting to region %s\n', region))
        cmd = 'samtools view %s %s -q %s %s | cut -f "3,4,9"' ## cmd line to grab the rname, pos, and tlen columns
        if (!file.exists(paste(bam.file, '.bam', sep = '')))
            if (file.exists(bai.file <- gsub('.bam$', '.bai', bam.file)))
            {
                .TMP.DIR = '~/temp/.samtools'
                system(paste('mkdir -p', TMP.DIR))
                tmp.fn = paste(normalizePath(TMP.DIR), '/tmp', runif(1), sep = '')
                system(sprintf('ln -s %s %s.bam', bam.file, tmp.fn))
                system(sprintf('ln -s %s %s.bam.bai', bai.file, tmp.fn))
            }

        cat('Calling', sprintf(cmd, st.flag, paste(tmp.fn, 'bam', sep = '.'), min.mapq, region), '\n')
        p = pipe(sprintf(cmd, st.flag, paste(tmp.fn, 'bam', sep = '.'), min.mapq, region), open = 'r')
    }
    else
    {
        cat('Calling', sprintf(cmd, st.flag, bam.file, min.mapq), '\n')
        p = pipe(sprintf(cmd, st.flag, bam.file, min.mapq), open = 'r')
    }

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
    if (verbose)
        cat('Starting fragment count on', bam.file, 'with bin size', window, 'and min mapQ', min.mapq, 'and insert size limit', max.tlen, 'with midpoint set to', midpoint, '\n')

    while (length(chunk <- readLines(p, n = chunksize))>0)
    {
        i = i+1

        if (fragments)
        {
            chunk = fread(paste(chunk, collapse = "\n"), header = F)[abs(V3)<=max.tlen, ]
            if (midpoint) ## only take midpionts
                chunk[, bin := 1 + floor((V2 + V3/2)/window)] ## use midpoint of template to index the correct bin
            else ## enumerate all bins containing fragment i.e. where fragments overlap multiple bins  (slightly slower)
            {
                if (verbose)
                    cat('!!!! Counting all overlapping bins !!!\n')
                chunk[, ":="(bin1 = 1 + floor((V2)/window), bin2 = 1 + floor((V2+V3)/window))]
                chunk = chunk[, list(V1, bin = bin1:bin2), by = list(ix = 1:length(V1))]
            }
        }
        else ## just count reads
        {
            cat('counting reads\n')
            chunk = fread(paste(chunk, collapse = "\n"), header = F)
            chunk[, bin := 1 + floor((V2)/window)]
        }

        tabs = chunk[, list(newcount = length(V1)), by = list(chr = as.character(V1), bin)] ## tabulate reads to bins data.table style
        counts[tabs, count := count + newcount] ## populate latest bins in master data.table

        ## should be no memory issues here since we preallocate the data table .. but they still appear
        if (do.gc)
        {
            print('GC!!')
            print(gc())
        }
                                        #          print(tables())

        ## report timing
        if (verbose)
        {
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
    if (verbose)
        cat("Finished computing coverage, and making GRanges\n")
    close(p)

    if (!is.null(region))
        system(sprintf('rm %s.bam %s.bam.bai', tmp.fn, tmp.fn))

    return(gr)
}

#' Compute rpkm counts from counts
#'
#' takes countbam (or bam.cov.gr) output "counts" and computes rpkm by aggregating across "by" variable
#' @param counts GRanges, data.table or data.frame with records, width fields
#' @param by Field to group counts by
#' @note The denominator (ie total reads) is just the sum of counts$records
#' @export
counts2rpkm = function(counts, by)
{
    out = aggregate(1:nrow(counts), by = list(by), FUN = function(x) sum(counts$records[x])/ sum(counts$width[x]/1000));
    out[,2] = out[,2]/sum(counts$records)*1e6;
    names(out) = c('by', 'rpkm');
    return(out);
}

#' Takes reads object and returns grl with each read and its mate (if exists)
#'
#' @param reads \code{GRanges} holding reads
#' @param as.grl Default TRUE. Return as a \code{GRangesList}
#' @param verbose Default FALSE
#' @name get.pairs.grl
#' @export
get.pairs.grl = function(reads, as.grl = TRUE, verbose = F)
{

    isdt <- inherits(reads, 'data.table')

    bad.col = c("seqnames", "ranges", "strand", "seqlevels",
                "seqlengths", "isCircular", "genome", "start", "end", "width", "element")

    if (verbose)
        cat('deduping\n')

    if (is(reads, 'GappedAlignmentPairs'))
        reads = unlist(reads)

    if (inherits(reads, 'GRanges')) {
        d <- duplicated(values(reads)$qname) ## duplicates are already paired up
        qpair.ix <- values(reads)$qname %in% unique(values(reads)$qname[d])
    } else if (isdt) {
        d <- duplicated(reads$qname)
        qpair.ix <- reads$qname %in% unique(reads$qname[d])
    }

    if (!inherits(reads, 'GenomicRanges') && !inherits(reads, 'data.table'))
    {
        if (verbose)
            cat('converting to granges\n')
        r.gr = granges(reads)
    }
    else if (!isdt)
        r.gr = reads[, c()]
    else
        r.gr <- reads

    if (verbose)
        cat('grbinding\n')

    m.gr = get.mate.gr(reads[!qpair.ix]);

    if (inherits(reads, 'GRanges')) {
        m.val <- values(m.gr)
        values(m.gr) = NULL;
        r.gr = c(r.gr, m.gr);
        mcols(r.gr) <- rrbind2(mcols(reads)[, setdiff(colnames(values(reads)), bad.col)], m.val)
    } else if (isdt) {
        m.gr <- m.gr[, setdiff(colnames(reads), colnames(m.gr)) := NA, with=FALSE]
        r.gr <- rbind(reads, m.gr, use.names=TRUE)
        setkey(r.gr, qname)
    }

    if (as.grl && !isdt) {
        if (verbose)
            cat('splitting\n')
        return(split(r.gr, as.character(r.gr$qname)))
    } else {
        return(r.gr)
    }
}

#' count.clips
#'
#' takes gr or gappedalignment object and uses cigar field (or takes character vector of cigar strings)
#' and returns data frame with fields (for character input)
#' $right.clips number of "right" soft clips (eg cigar 89M12S)
#' #left.clips number of "left" soft clips (eg cigar 12S89M)
#' or appends these fields to the reads object
#'
#' @param reads GenomicRanges holding the reads
#' @param hard [Default TRUE] option counts hard clips
#' @name count.clips
#' @export
count.clips = function(reads, hard = FALSE)
{
    if (length(reads) == 0)
        return(reads)
    if (inherits(reads, 'GRanges') | inherits(reads, 'GappedAlignments'))
        cigar = values(reads)$cigar
    else
        cigar = reads;

    if (!inherits(cigar, 'character') & !inherits(cigar, 'factor'))
        stop('Input must be GRanges, GappedAlignments, or character vector')

    out = data.frame(left.clips = rep(0, length(cigar)), right.clips = rep(0, length(cigar)));

    re.left = '^(\\d+)S.*'
    re.right = '.*[A-Z](\\d+)S$'

    lclip.ix = grep(re.left, cigar)
    rclip.ix = grep(re.right, cigar)

    if (length(lclip.ix)>0)
    {
        left.clips = gsub(re.left, '\\1', cigar[lclip.ix])
        out$left.clips[lclip.ix] = as.numeric(left.clips)
    }

    if (length(rclip.ix)>0)
    {
        right.clips = gsub(re.right, '\\1', cigar[rclip.ix])
        out$right.clips[rclip.ix] = as.numeric(right.clips)
    }

    if (inherits(reads, 'GRanges') | inherits(reads, 'GappedAlignments'))
    {
        values(reads)$right.clips = out$right.clips
        values(reads)$left.clips = out$left.clips
        out = reads
    }

    return(out)
}



#' varbase
#'
#' takes gr or gappedalignment object "reads" and uses cigar, MD, seq fields
#' to return variant bases and ranges
#'
#' returns grl (of same length as input) of variant base positions with character vector $varbase field populated with variant bases
#' for each gr item in grl[[k]], with the following handling for insertions, deletions, and substitution gr's:
#'
#' substitutions: nchar(gr$varbase) = width(gr) of the corresponding var
#' insertions: nchar(gr$varbase)>=1, width(gr) ==0
#' deletions: gr$varbase = '', width(gr)>=1
#'
#' Each gr also has $type flag which shows the cigar string code for the event ie
#' S = soft clip --> varbase represents clipped bases
#' I = insertion --> varbase represents inserted bases
#' D = deletion --> varbase is empty
#' X = mismatch --> varbase represents mismatched bases
#' @param reads GenomicRanges to extract variants from
#' @param soft [Default TRUE]
#' @param verbose [Default TRUE]
#' @name varbase
#' @export
varbase = function(reads, soft = TRUE, verbose = TRUE)
{
    nreads = length(reads)
    if (inherits(reads, 'GRangesList'))
    {
        was.grl = TRUE
        r.id = grl.unlist(reads)$grl.ix
        reads = unlist(reads)
    }
    else if (inherits(reads, 'data.frame'))
    {
        r.id = 1:nrow(reads)
        nreads = nrow(reads)
        was.grl = FALSE
    }
    else
    {
        r.id = 1:length(reads)
        was.grl = FALSE
    }

    if (!inherits(reads, 'GRanges') & !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame'))
        stop('Reads must be either GRanges, GRangesList, or GappedAlignment object')
    else if (inherits(reads, 'data.frame'))
    {
#        if (is.null(reads$cigar) | is.null(reads$seq))
#            stop('Reads must have cigar and seq fields specified')
    }
 #   else if (is.null(values(reads)$cigar) | is.null(values(reads)$seq))
 #       stop('Reads must have cigar and seq fields specified')

    if (is.data.frame(reads))
    {
        sl = NULL
        sn =  reads$seqnames
        cigar = as.character(reads$cigar)
        seq = as.character(reads$seq)
        str = reads$strand

        if (!is.null(reads$MD))
            md = as.character(reads$MD)
        else
            md = rep(NA, length(cigar))
    }
    else
    {
        sl = seqlengths(reads)
        sn =  seqnames(reads)
        cigar = as.character(values(reads)$cigar)
        seq = as.character(values(reads)$seq)
        str = as.character(strand(reads))

        if (!is.null(values(reads)$MD))
            md = as.character(values(reads)$MD)
        else
            md = rep(NA, length(cigar))
    }

    if (!inherits(cigar, 'character') & !inherits(cigar, 'character') & !inherits(md, 'character'))
        stop('Input must be GRanges with seq, cigar, and MD fields populated or GappedAlignments object')

    ix = which(!is.na(cigar))

    if (length(ix)==0)
        return(rep(GRangesList(GRanges()), nreads))

    cigar = cigar[ix]
    seq = seq[ix]
    md = md[ix]
    str = str[ix]

    if (is.data.frame(reads))
    {
        r.start = reads$start[ix]
        r.end = reads$end[ix]
    }
    else
    {
        r.start = start(reads)[ix]
        r.end = end(reads)[ix];
    }

    flip = str == '-'

    if (!is.null(seq))
        {
            nix = sapply(seq, function(x) all(is.na(x)))
            if (any(nix))
                seq[nix] = ''
            seq = strsplit(seq, '')
        }

    cigar.vals = explodeCigarOps(cigar)
    cigar.lens = explodeCigarOpLengths(cigar)
    
    ## cigar.vals = lapply(strsplit(cigar, "\\d+"), function(x) x[2:length(x)])
    ## cigar.lens = lapply(strsplit(cigar, "[A-Z]"), as.numeric)

    clip.left = sapply(cigar.vals, function(x) x[1] %in% c('H', 'S'))
    clip.right = sapply(cigar.vals, function(x) x[length(x)] %in% c('H', 'S'))

    if (any(clip.left))
        r.start[clip.left] = r.start[clip.left]-sapply(cigar.lens[which(clip.left)], function(x) x[1])

    if (any(clip.right))
        r.end[clip.right] = r.end[clip.right]+sapply(cigar.lens[which(clip.right)], function(x) x[length(x)])

                                        # split md string into chars after removing "deletion" signatures and also
                                        # any soft clipped base calls (bases followed by a 0)
    ##    md.vals = strsplit(gsub('([ATGCN])', '|\\1|', gsub('\\^[ATGCN]+', '|', md)), '\\|')
    md.vals = strsplit(gsub('([A-Z])', '|\\1|', gsub('\\^[A-Z]+', '|', md)), '\\|')

                                        # ranges of different cigar elements relative to query ie read-centric coordinates
    starts.seq = lapply(1:length(cigar.lens), function(i)
    {
        x = c(0, cigar.lens[[i]])
        x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))+1] = 0  ## deletions have 0 width on query
        cumsum(x[1:(length(x)-1)])+1
    })

    ends.seq = lapply(1:length(cigar.lens), function(i)
    {
        x = cigar.lens[[i]];
        x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))] = 0
        cumsum(x)
    })




                                        # ranges of different cigar elements relative to reference coordinatse
    starts.ref = lapply(1:length(cigar.lens), function(i)
    {
        x = c(0, cigar.lens[[i]]);
        x[which(cigar.vals[[i]] %in% c('I'))+1] = 0 ## insertions have 0 width on reference / subject
        cumsum(x[1:(length(x)-1)]) + r.start[i]
    })

    ends.ref = lapply(1:length(cigar.lens), function(i)
    {
        x = cigar.lens[[i]];
        x[which(cigar.vals[[i]] %in% c('I'))] = 0
        cumsum(x) + r.start[i] - 1
    })

                                        # now using MD find coordinates of mismatched bases (using starts and ends of M regions)

                                        # find coords of subs on genome
    tmp.pos = lapply(1:length(md.vals), function(i)
    {
        x = md.vals[[i]]
                                        #        nix = grepl('[ATGCN]', x);
        nix = grepl('[A-Z]', x);
        if (!any(nix))
            return(c())
        p = rep(0, length(x))
        p[!nix] = as.numeric(x[!nix])
        p[nix] = 1
        s.pos.m = cumsum(p)[nix] ## position of subs in read
        mix = cigar.vals[[i]]=='M'
        m.st = cumsum(c(1, ends.seq[[i]][mix]-starts.seq[[i]][mix]+1))
        m.st.g = starts.ref[[i]][mix]
        m.st.r = starts.seq[[i]][mix]
        s.match = rep(NA, length(s.pos.m)) ## indices of matching "M" cigar element for each sub
        for (ii in 1:length(s.pos.m))
        {
            j = 0;
            done = FALSE
            for (j in 0:(length(m.st)-1))
                if (s.pos.m[ii] < m.st[j+1])
                    break
            s.match[ii] = j
        }

        s.pos.g = m.st.g[s.match] + s.pos.m-m.st[s.match]
        s.pos.r = m.st.r[s.match] + s.pos.m-m.st[s.match]

        return(rbind(s.pos.g, s.pos.r))
    })
    subs.pos = lapply(tmp.pos, function(x) x[1,])
    subs.rpos = lapply(tmp.pos, function(x) x[2,])
    if (is.null(seq))
        subs.base = lapply(lapply(md.vals, grep, pattern = '[ATGCN]', value = T), function(x) rep('X', length(x))) ## replace with N
    else        
        subs.base = lapply(1:length(seq), function(x) ifelse(is.na(seq[[x]][subs.rpos[[x]]]), 'X', seq[[x]][subs.rpos[[x]]]))


                                        # make sure md and cigar are consistent
                                        # (for some reason - sometimes soft clipped mismatches are included in MD leading to a longer MD string)
                                        # also some MD are NA
    mlen.cigar = sapply(1:length(ends.seq), function(x) {mix = cigar.vals[[x]]=='M'; sum(ends.seq[[x]][mix]-starts.seq[[x]][mix]+1)})
                                        #    mlen.md = sapply(md.vals, function(x) {ix = grepl('[ATGCN]', x); sum(as.numeric(x[!ix])) + sum(nchar(x[ix]))})
    mlen.md = sapply(md.vals, function(x) {ix = grepl('[A-Z]', x); sum(as.numeric(x[!ix])) + sum(nchar(x[ix]))})
    good.md = which(!is.na(md))
                                        #  good.md = which(mlen.md == mlen.cigar & !is.na(md))

    if (any(na <- is.na(md)))
    {
        warning('MD field absent from one or more input reads')
        good.md = which(!na)
    }

                                        # now make granges of subs
    if (length(good.md)>0)
    {
        if (any(mlen.md[good.md] != mlen.cigar[good.md]))
            warning('the lengths of some MD strings do not match the number of M positions on the corresponding CIGAR string: some variants may not be correctly mapped to the genome')

        iix.md = unlist(lapply(good.md, function(x) rep(x, length(subs.pos[[x]]))))
        tmp = unlist(subs.pos[good.md])

        if (!is.null(tmp))
        {
            subs.gr = GRanges(sn[ix][iix.md], IRanges(tmp, tmp), strand = '*')
            values(subs.gr)$varbase = unlist(subs.base[good.md])
            values(subs.gr)$varlen = 1
            values(subs.gr)$type = 'X'
            values(subs.gr)$iix = ix[iix.md]
        }
        else
            subs.gr = GRanges()
    }
    else
        subs.gr = GRanges()

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
    str <- str[iix] # JEREMIAH

    if (length(cigar.vals)>0)
        {
        var.seq = lapply(1:length(cigar.vals),
                         function(i)
                             {                                 
                             if (ends.seq[i]<starts.seq[i])
                                 return('') # deletion
                             else
                                 if (length(seq[[iix[i]]])==0)
                                     rep('N', ends.seq[i]-starts.seq[i]+1)
                                 else
                                     seq[[iix[i]]][starts.seq[i]:ends.seq[i]] #insertion
                         })
        other.gr = GRanges(sn[ix][iix], IRanges(starts.ref, ends.ref), strand = str, seqlengths = sl)        
        values(other.gr)$varbase = sapply(var.seq, paste, collapse = '')
        values(other.gr)$varlen = cigar.lens
        values(other.gr)$type = cigar.vals
        values(other.gr)$iix = ix[iix];
        out.gr = sort(c(subs.gr, other.gr))
    }
    else
        out.gr = subs.gr
    
                                        # add default colors to out.gr
    VAR.COL = get.varcol()

    col.sig = as.character(out.gr$type)
    xix = out.gr$type == 'X'
    col.sig[xix] = paste(col.sig[xix], out.gr$varbase[xix], sep = '')
    out.gr$col = VAR.COL[col.sig]
    out.gr$border = out.gr$col

    if (!soft)
        if (any(soft.ix <<- out.gr$type %in% c('H', 'S')))
            out.gr = out.gr[-which(soft.ix)]

    out.grl = rep(GRangesList(GRanges()), nreads)
    out.iix = r.id[values(out.gr)$iix]
    values(out.gr)$iix = NULL

    tmp.grl = GenomicRanges::split(out.gr, out.iix)
    out.grl[as.numeric(names(tmp.grl))] = tmp.grl
    values(out.grl)$qname[r.id] = reads$qname

    return(out.grl)
}

#' splice.cigar
#'
#' takes gr or gappedalignment object "reads" and parses cigar fields
#' to return grl corresponding to spliced alignments on the genome corresponding to
#' portions of the cigar
#'
#' ie each outputted grl item contains the granges corresponding to all non-N portions of cigar string
#'
#' if grl provided as input (e.g. paired reads) then all of the spliced ranges resulting from each
#' input grl item will be put into the corresponding output grl item
#'
#' NOTE: does not update MD tag
#'
#' if use.D = TRUE, then will treat "D" (deletion) in addition to "N" flags as indicative of deletion event.
#' @param reads \code{GRanges} reads
#' @param verbose Default TRUE
#' @param fast Default TRUE
#' @param use.D Default TRUE
#' @param rem.soft Default TRUE
#' @param get.seq Default FALSE
#' @param return.grl Default TRUE
#' @name splice.cigar
#' @export
splice.cigar = function(reads, verbose = TRUE, fast = TRUE, use.D = TRUE, rem.soft = TRUE, get.seq = FALSE, return.grl = TRUE)
{
    nreads = length(reads)

    if (nreads==0)
        if (return.grl)
            return(GRangesList())
        else
            return(GRanges)

    if (inherits(reads, 'GRangesList'))
    {
        was.grl = TRUE
        r.id = as.data.frame(reads)$element
        reads = unlist(reads)
    }
    else
    {
        r.id = 1:length(reads)
        was.grl = FALSE
    }


    if (is.data.frame(reads))
    {
        sl = NULL
        sn =  reads$seqnames
        cigar = as.character(reads$cigar)
        seq = as.character(reads$seq)
        str = reads$strand

        if (!is.null(reads$MD))
            md = as.character(reads$MD)
        else
            md = rep(NA, length(cigar))
    }
    else
    {
        sl = seqlengths(reads)
        sn =  seqnames(reads)
        cigar = as.character(values(reads)$cigar)
        seq = as.character(values(reads)$seq)
        str = as.character(strand(reads))

        if (!is.null(values(reads)$MD))
            md = as.character(values(reads)$MD)
        else
            md = rep(NA, length(cigar))
    }


    if (!inherits(reads, 'GRanges') & !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame'))
        stop('Reads must be either GRanges, GRangesList, or GappedAlignment object')
    else if (is.null(values(reads)$cigar) | is.null(values(reads)$seq))
        stop('Reads must have cigar and seq fields specified')

    if (!inherits(cigar, 'character') & !inherits(cigar, 'character') & !inherits(md, 'character'))
        stop('Input must be GRanges with seq, cigar, and MD fields populated or GappedAlignments object')


    ix = which(!is.na(reads$cigar))

    if (length(ix)==0)
        if (return.grl)
            return(rep(GRangesList(GRanges()), nreads))
        else
            return(GRanges())

    if (fast)
    {
        ir = cigarRangesAlongReferenceSpace(reads[ix]$cigar, N.regions.removed = FALSE, with.ops = TRUE, reduce.ranges = FALSE)
        irul = unlist(ir)
        out.gr = GRanges(rep(seqnames(reads)[ix], elementLengths(ir)), shift(IRanges(irul), rep(start(reads)[ix]-1, elementLengths(ir))),
                         strand = rep(strand(reads)[ix], elementLengths(ir)), seqlengths = seqlengths(reads))
        out.gr$type = names(irul)
        out.gr$rid = ix[rep(1:length(ir), elementLengths(ir))]
        out.gr$riid = unlist(lapply(elementLengths(ir), function(x) 1:x))
        out.gr$fid = r.id[out.gr$rid]
        out.gr$qname = reads$qname[out.gr$rid]

        if (return.grl)
        {
            out.grl = rep(GRangesList(GRanges()), nreads)
            tmp.grl = split(out.gr, out.gr$fid)
            out.grl[as.numeric(names(tmp.grl))] = tmp.grl
            return(out.grl)
        }
        else
            return(out.gr)
    }
    else
    {

        cigar = cigar[ix]
        str = str[ix]

        if (is.data.frame(reads))
        {
            r.start = reads$start[ix]
            r.end = reads$end[ix]
        }
        else
        {
            r.start = start(reads)[ix]
            r.end = end(reads)[ix];
        }

        flip = str == '-'
        cigar.vals = lapply(strsplit(cigar, "\\d+"), function(x) x[2:length(x)])
        cigar.lens = lapply(strsplit(cigar, "[A-Z]"), as.numeric)

        clip.left = sapply(cigar.vals, function(x) x[1] == 'S')
        clip.right = sapply(cigar.vals, function(x) x[length(x)] == 'S')

                                        # ranges of different cigar elements relative to query ie read-centric coordinates
        starts.seq = lapply(1:length(cigar.lens), function(i)
        {
            x = c(0, cigar.lens[[i]])
            x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))+1] = 0  ## deletions have 0 width on query
            cumsum(x[1:(length(x)-1)])+1
        })

        ends.seq = lapply(1:length(cigar.lens), function(i)
        {
            x = cigar.lens[[i]];
            x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))] = 0
            cumsum(x)
        })

                                        # ranges of different cigar elements relative to reference coordinatse
        starts.ref = lapply(1:length(cigar.lens), function(i)
        {
            x = c(0, cigar.lens[[i]]);
            x[which(cigar.vals[[i]] %in% c('I'))+1] = 0 ## insertions have 0 width on reference / subject
            cumsum(x[1:(length(x)-1)]) + r.start[i]
        })

        ends.ref = lapply(1:length(cigar.lens), function(i)
        {
            x = cigar.lens[[i]];
            x[which(cigar.vals[[i]] %in% c('I'))] = 0
            Cumsum(x) + r.start[i] - 1
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

        if (use.D)
            splice.char = c(splice.char, 'D')

        if (rem.soft)
            splice.char = c(splice.char, 'S')

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
        if (length(cigar.vals)>0)
            {
            other.gr = GRanges(sn[ix][iix], IRanges(starts.ref, ends.ref), strand = str, seqlengths = sl)

            if (get.seq)
            {
                var.seq = lapply(1:length(cigar.vals),
                                 function(i)
                                 {
                                     if (ends.seq[i]<starts.seq[i])
                                         return('') # deletion
                                     else
                                         seq[[iix[i]]][starts.seq[i]:ends.seq[i]] #insertion
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

        if (return.grl)
        {
            out.grl = rep(GRangesList(GRanges()), nreads)
            tmp.grl = split(out.gr, out.iix)
            out.grl[as.numeric(names(tmp.grl))] = tmp.grl
            return(out.grl)
        }
        else
            return(out.gr)
    }
}

#' bamflag
#'
#' shortcut .. assumes reads are GappedAlignments with flag variable or actual integers representing bam flag
#' @param reads GenomicRanges holding the reads
#' @name bamflag
#' @export
bamflag = function(reads)
{
    if (inherits(reads, 'GappedAlignments') | inherits(reads, 'data.frame') | inherits(reads, 'GRanges'))
        bf = reads$flag
    else
        bf = reads

    out = matrix(as.numeric(intToBits(bf)), byrow = T, ncol = 32)[, 1:12, drop = FALSE]
    colnames(out) = c('isPaired', 'isProperPair', 'isUnmappedQuery', 'hasUnmappedMate', 'isMinusStrand', 'isMateMinusStrand', 'isFirstMateRead', 'isSecondMateRead', 'isNotPrimaryRead', 'isNotPassingQualityControls', 'isDuplicate', 'isSupplementary')

    return(out)
                                        #    if (inherits(reads, 'GappedAlignments'))
                                        #      return(bamFlagAsBitMatrix(values(reads)$flag))
                                        #    else
                                        #      return(bamFlagAsBitMatrix(reads))
}


#' bamtag
#'
#' outputs a tag that cats qname, first vs first second mate +/- secondary alignment +/- gr.string
#' to give an identifier for determine duplicates in a read pile
#' @param reads GenomicRanges holding the reads
#' @name bamflag
#' @export
bamtag = function(reads, secondary = F, gr.string = F)
{
    grs = sec = NULL
    if (secondary)
        sec = bamflag(read$flag[, 'isNotPrimaryRead'])

    if (gr.string)
        grs = gr.string(reads, mb  = F)

    return(paste(reads$qname, ifelse(bamflag(reads$flag)[, 'isFirstMateRead'], '1', '2'), grs, sec, sep = '_'))
}


#' Count bases in cigar string
#'
#' Counts the total number of bases, per cigar, that fall into D, I, M, S categories.
#' countCigar makes no distinction between, for instance 1S2M2S, 2S2M1S, or 3S2M
#' @param cigar character vector of cigar strings
#' @return a 4-column, length(cigar)-row matrix with the total counts for each type
#' @export
countCigar <- function(cigar) {

    cigar.vals <- unlist(strsplit(cigar, "\\d+"))
    cigar.lens <- strsplit(cigar, "[A-Z]")
    lens <- nchar(gsub('\\d+', '', cigar))
    lens[is.na(cigar)] <- 1

    cigar.lens <- as.numeric(unlist(cigar.lens))
    cigar.vals <- cigar.vals[cigar.vals != ""]
    repr       <- rep(seq_along(cigar), lens)
    dt         <- data.table(val=cigar.vals, lens=cigar.lens, group=repr, key="val")
    
    smr.d      <- dt["D",][, sum(lens), by=group]
    smr.i      <- dt["I",][, sum(lens), by=group]
    smr.m      <- dt["M",][, sum(lens), by=group]
    smr.s      <- dt["S",][, sum(lens), by=group]

    out <- matrix(nrow=length(cigar), ncol=4, 0)
    out[smr.d$group,1] <- smr.d$V1
    out[smr.i$group,2] <- smr.i$V1
    out[smr.m$group,3] <- smr.m$V1
    out[smr.s$group,4] <- smr.s$V1
    colnames(out) <- c('D','I','M','S')

    return(out)
}



########################
#' get.var.col
#'
#' simple function storing default
#' variant color scheme
#' @name get.var.col
#' @export
########################
get.varcol = function()
  {
    VAR.COL = c('XA' = 'green', 'XG' = 'brown', 'XC' = 'blue', 'XT' = 'red', 'D' = 'white', 
    'I'= 'purple', 'N' = alpha('gray', 0.2), 'XX' = 'black', 'S' = alpha('pink', 0.9))
    return(VAR.COL)
  }



#' is.paired.end
#'
#' @description
#'
#' Check if bam file is paired end by using 0x1 flag
#' 
#' @name is.paired.end
#' @export
is.paired.end = function(bams)
    {
        out = sapply(bams, function(x)            
            {
                if (is.na(x))
                    return(NA)
                if (!file.exists(x))
                    return(NA)
                out = FALSE                
                p = pipe(sprintf('samtools view -h  %s | head -n 100 | samtools view -f 0x1 - | wc -l', x))
                ln = as.numeric(readLines(p))
                out = ln>0
                close(p)
                return(out)                
            })     
        return(out)
    }

alpha = function(col, alpha)
{
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}
