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
