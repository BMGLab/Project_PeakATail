#!/usr/bin/env Rscript
## Convert scAPAtrap 0.2.0 output to the benchmark's standard PAS BED6.
##
## Usage:
##   Rscript convert_scapatrap.R <outdir_with_rda> <workdir_with_.peaks> <label> <genome> [beddir]
##
## Emits into [beddir], defaulting to <outdir_with_rda>:
##   pas.bed      PRIMARY -- quantified/filtered peak set (scAPAtrapData.rda peaks.meta)
##   pas_all.bed  VARIANT -- every called peak (peaks.saf), pre-quantification
##
## Geometry is scAPAtrap's OWN PAS definition, not our reinterpretation:
## scAPAtrap:::.computePAcoord() sets coord = end ('+') / start ('-'), and that
## `coord` is what scAPAtrap itself feeds to .findIntersectBypeak() when matching
## peaks to polyA tails.  See the header written into each BED for the full proof.

suppressMessages({library(scAPAtrap); library(Matrix)})
invisible(Sys.setlocale("LC_ALL", "C"))

args    <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(4L, 5L))
OUTDIR  <- args[1]; WORKDIR <- args[2]; LABEL <- args[3]; GENOME <- args[4]
## BED destination: the benchmark's convention is <dataset>/<tool>[/<sample>]/pas.bed,
## which for pbmc_10k_v3 is the scapatrap/ dir, not the scAPAtrap outputDir (trap/).
BEDDIR  <- if (length(args) == 5L) args[5] else OUTDIR
dir.create(BEDDIR, showWarnings = FALSE, recursive = TRUE)

rdaFile <- file.path(OUTDIR, "scAPAtrapData.rda")
safFile <- file.path(OUTDIR, "peaks.saf")
fwdPk   <- file.path(WORKDIR, "input.UniqSorted.dedup.forward.bam.peaks")
revPk   <- file.path(WORKDIR, "input.UniqSorted.dedup.reverse.bam.peaks")
for (f in c(rdaFile, safFile, fwdPk, revPk)) stopifnot(file.exists(f))

## ---------------------------------------------------------------- helpers ----
## BED is 0-based half-open; SAF/GRanges coords are 1-based inclusive.
## A 1-based point P becomes the 1-bp interval [P-1, P).
to_point <- function(coord) data.frame(start = coord - 1L, end = coord)

## derfinder::fullCoverage() renames the requested Ensembl contigs to UCSC style
## (`names(result) <- extendedMapSeqlevels(chrs, ...)`), so scAPAtrap emits
## 'chr1' even though the BAM contigs were '1'.  Strip it back: the benchmark's
## reference (atlas/TES/chrom.sizes) is bare-Ensembl.  Lossless here because the
## runs restricted chrs to the autosomes + X + Y (no chrM/MT involved).
strip_chr <- function(x) sub("^chr", "", as.character(x))

## Write '#'-commented header + C-locale-sorted BED body.
write_bed <- function(df, path, header) {
    stopifnot(all(df$end == df$start + 1L), all(df$start >= 0L))
    tmp <- paste0(path, ".body.tmp")
    write.table(df[, c("chrom","start","end","name","score","strand")], tmp,
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    writeLines(header, path)
    rc <- system(sprintf("LC_ALL=C sort -k1,1 -k2,2n %s >> %s", shQuote(tmp), shQuote(path)))
    if (rc != 0L) stop("sort failed for ", path)
    unlink(tmp)
    invisible(nrow(df))
}

score_cap <- function(x) pmin(1000L, as.integer(round(pmax(0, x))))

## -------------------------------------------------- PRIMARY: pas.bed ---------
e <- new.env(); load(rdaFile, envir = e)
o <- get("scAPAtrapData", envir = e); rm(e)
meta <- o$peaks.meta
stopifnot(identical(meta$peakID, rownames(o$peaks.count)))

## Re-derive coord from start/end and assert it equals the tool's own column --
## this is the executable proof that our reduction == .computePAcoord().
chk <- ifelse(meta$strand == "+", meta$end, meta$start)
stopifnot(identical(as.numeric(chk), as.numeric(meta$coord)))
stopifnot(all(meta$strand %in% c("+", "-")))

umis <- Matrix::rowSums(o$peaks.count)
n_cells <- ncol(o$peaks.count)
rm(o); invisible(gc())

pt <- to_point(meta$coord)
prim <- data.frame(chrom  = strip_chr(meta$chr),
                   start  = as.integer(pt$start),
                   end    = as.integer(pt$end),
                   name   = meta$peakID,
                   score  = score_cap(umis),
                   strand = meta$strand,
                   stringsAsFactors = FALSE)

dup_key <- paste(prim$chrom, prim$start, prim$strand, sep = ":")
n_dup   <- sum(duplicated(dup_key))

## -------------------------------------------------- VARIANT: pas_all.bed -----
saf <- read.delim(safFile, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("peakID","chr","start","end","strand"))
fw  <- read.delim(fwdPk, header = TRUE, stringsAsFactors = FALSE)
rv  <- read.delim(revPk, header = TRUE, stringsAsFactors = FALSE)
pk  <- rbind(fw, rv)                       # generateSAF(): rbind(forward, reverse)
stopifnot(nrow(pk) == nrow(saf))
## generateSAF() assigns PeakID = paste0("peak_", 1:n) over exactly this row
## order, so `value` (derfinder regionMatrix mean bp coverage) maps positionally.
stopifnot(identical(as.character(pk$chr),    as.character(saf$chr)),
          identical(as.numeric(pk$start),    as.numeric(saf$start)),
          identical(as.numeric(pk$end),      as.numeric(saf$end)),
          identical(as.character(pk$strand), as.character(saf$strand)),
          identical(saf$peakID, paste0("peak_", seq_len(nrow(saf)))))

saf_coord <- ifelse(saf$strand == "+", saf$end, saf$start)
pt2 <- to_point(saf_coord)
allb <- data.frame(chrom  = strip_chr(saf$chr),
                   start  = as.integer(pt2$start),
                   end    = as.integer(pt2$end),
                   name   = saf$peakID,
                   score  = score_cap(pk$value),
                   strand = saf$strand,
                   stringsAsFactors = FALSE)
n_dup_all <- sum(duplicated(paste(allb$chrom, allb$start, allb$strand, sep = ":")))

## Every primary peak must be a subset of the SAF peak set, unchanged.
stopifnot(all(prim$name %in% allb$name))
m <- match(prim$name, allb$name)
stopifnot(identical(prim$start[order(prim$name)], allb$start[m][order(prim$name)]))

## ------------------------------------------------------------- headers -------
ver <- as.character(packageVersion("scAPAtrap"))
geom <- function(kind, nrows, extra) c(
sprintf("# scAPAtrap %s PAS point BED6 -- %s, %s (%s)", ver, LABEL, kind, GENOME),
"# BED6, 0-based half-open, 1-bp intervals. Point = strand-aware 3'-most base of the peak.",
"#",
"# GEOMETRY PROOF (installed scAPAtrap 0.2.0 in conda env bench_scapatrap,",
"# sources dumped from the lazy-load DB via print()/getFromNamespace()):",
"#",
"# 1) COORDINATE SYSTEM: 1-based, fully closed.  findPeaks() takes regions from",
"#    derfinder::regionMatrix() as a GRanges (1-based inclusive) and returns",
"#    as.data.frame(DefPeak)[, c(chr,strand,start,end,value)].  The wide-peak",
"#    split branch stays in the same frame: rawstart <- start(regions[j]) - 1,",
"#    bumphunter::regionFinder() runs on pos = 1..width, then",
"#    split.peaks$start/end <- split.peaks$start/end + rawstart, restoring",
"#    1-based genomic coordinates.  generateSAF() only rbind()s forward+reverse,",
"#    adds PeakID, and write.table()s columns PeakID/chr/start/end/strand with",
"#    no header -- coordinates pass through untouched (SAF is 1-based inclusive).",
"#",
"# 2) WHICH BASE IS THE PAS -- scAPAtrap's OWN definition, not our inference.",
"#    scAPAtrap:::.computePAcoord(peak):",
"#        peak$coord <- 0",
"#        peak$coord[peak$strand == '+'] <- peak$end  [peak$strand == '+']",
"#        peak$coord[peak$strand == '-'] <- peak$start[peak$strand == '-']",
"#    generatescExpMa() calls it on every peak and stores the result as the",
"#    `coord` column of peaks.meta in scAPAtrapData.rda.  That same `coord` is",
"#    the ONLY coordinate scAPAtrap uses to decide whether a peak carries a",
"#    polyA tail: .findIntersectBypeak() builds GRanges(start=coord, end=coord)",
"#    for both peaks and tails and findOverlaps() them.  So `coord` IS the",
"#    tool's polyA site.  Verified executably: ifelse(strand=='+', end, start)",
"#    reproduces peaks.meta$coord for all rows (stopifnot in the converter).",
"#",
"# 3) BED CONVERSION: 1-based inclusive point P -> [P-1, P).",
"#      strand '+': P = end   -> chrom  end-1    end    +",
"#      strand '-': P = start -> chrom  start-1  start  -",
"#    Identical reduction to PeakATail/polyApipe/Sierra and idempotent under",
"#    score_tool.py make_point() ('+': $3-1,$3 ; '-': $2,$2+1).",
"#",
"# 4) CHROMOSOME NAMES: the BAM contigs were bare Ensembl ('1'..'22','X','Y'),",
"#    but derfinder::fullCoverage() renames them via",
"#      names(result) <- extendedMapSeqlevels(chrs, ...)   [UCSC style]",
"#    so scAPAtrap's .peaks/.saf/.rda all carry 'chr1'.  (findPeaks() has a",
"#    paste0('Chr', .) / gsub('^Chr','', .) round-trip, but it never fires here:",
"#    grep('^[0-9]', names(fullCov)) is empty once derfinder has prefixed 'chr'.)",
"#    Stripped back to bare Ensembl for the benchmark reference; lossless --",
"#    the runs restricted chrs to autosomes + X + Y, so no chrM/MT ambiguity.",
"#",
"# 5) NO POLYA-TAIL-ANCHORED COORDINATE EXISTS FOR THIS RUN.  TRAP.PARAMS()",
"#    defaults tails.search = 'no' and neither driver overrode it, so the",
"#    scAPAtrap() wrapper took the `else { tailsfile = NULL }` branch: findTails/",
"#    findTailsByPeaks never ran, generatescExpMa(tails = NULL) skipped the",
"#    .findIntersectBypeak() filter, and no *.tails file was written.  There is",
"#    therefore no pas_tail.bed variant -- these peaks are coverage-only, exactly",
"#    like the coverage callers scAPAtrap advertises itself against.",
"#",
sprintf("# SCORE COLUMN: %s", extra),
sprintf("# n_sites=%d", nrows))

hdr_prim <- c(geom("PRIMARY (quantified peak set)", nrow(prim),
   sprintf("min(1000, total UMIs over the %d filtered cells).", n_cells)),
sprintf("# SOURCE: %s  (peaks.meta + peaks.count)", rdaFile),
"# PEAK SET: generatescExpMa() -> reducePeaks(min.cells=10, min.count=10) applied",
"#   to umi_tools-counted peak x cell UMIs restricted to the whitelisted barcodes;",
"#   this is scAPAtrap's final deliverable (the peaks it actually quantifies).",
sprintf("# duplicate (chrom,pos,strand) points retained as separate rows: %d", n_dup),
"#")

hdr_all <- c(geom("VARIANT (all called peaks, pre-quantification)", nrow(allb),
   "min(1000, round(derfinder regionMatrix mean bp coverage `value`))."),
sprintf("# SOURCE: %s + the two .peaks tables in %s", safFile, WORKDIR),
"# PEAK SET: every peak findPeaksByStrand() called, BEFORE the min.cells=10 /",
"#   min.count=10 UMI filter.  Sensitivity variant only -- pas.bed is PRIMARY.",
sprintf("# duplicate (chrom,pos,strand) points retained as separate rows: %d", n_dup_all),
"#")

n1 <- write_bed(prim, file.path(BEDDIR, "pas.bed"),     hdr_prim)
n2 <- write_bed(allb, file.path(BEDDIR, "pas_all.bed"), hdr_all)

cat(sprintf("[%s] pas.bed=%d (cells=%d, dup_pts=%d)  pas_all.bed=%d (dup_pts=%d)\n",
            LABEL, n1, n_cells, n_dup, n2, n_dup_all))
