##
##
## $$$$$$$\                                   $$\                        $$\                 $$\   $$\   $$\     $$\ $$\
## $$  __$$\                                  $$ |                       $$ |                $$ |  $$ |  $$ |    \__|$$ |
## $$ |  $$ | $$$$$$$\ $$$$$$\  $$$$$$\$$$$\$$$$$$\   $$$$$$\   $$$$$$\  $$ | $$$$$$$\       $$ |  $$ |$$$$$$\   $$\ $$ |
## $$$$$$$  |$$  _____|\____$$\ $$  _$$  _$$\_$$  _| $$  __$$\ $$  __$$\ $$ |$$  _____|      $$ |  $$ |\_$$  _|  $$ |$$ |
## $$  __$$< \$$$$$$\  $$$$$$$ |$$ / $$ / $$ |$$ |   $$ /  $$ |$$ /  $$ |$$ |\$$$$$$\        $$ |  $$ |  $$ |    $$ |$$ |
## $$ |  $$ | \____$$\$$  __$$ |$$ | $$ | $$ |$$ |$$\$$ |  $$ |$$ |  $$ |$$ | \____$$\       $$ |  $$ |  $$ |$$\ $$ |$$ |
## $$ |  $$ |$$$$$$$  \$$$$$$$ |$$ | $$ | $$ |\$$$$  \$$$$$$  |\$$$$$$  |$$ |$$$$$$$  |      \$$$$$$  |  \$$$$  |$$ |$$ |
## \__|  \__|\_______/ \_______|\__| \__| \__| \____/ \______/  \______/ \__|\_______/        \______/    \____/ \__|\__|
##
##
## Rsamtools util
##
## wrapper functions around Rsamtools and rtracklayer
## to help extract read info from bam file, mutations pileups, and coverage from wig / bigwig files
##
##


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
            intervals = seqinfo2gr(seqinfo(bam))
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

    keep = which(seqnames(gr) %in% seqlevels(bam))

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
