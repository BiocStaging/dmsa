# ============================================================================
# GENE MODELS: exons, introns, UTRs and strand, on real genomic coordinates
#
# The locus panel used to draw a gene as a featureless grey bar, and - when a
# probe had no coordinate - to space probes evenly along it. Both are lies of a
# particular kind: they look like a figure that says where something is while
# saying nothing of the sort. A CpG in the first exon and a CpG 40 kb into an
# intron are different claims about regulation, and DMSA's whole premise is
# that WHERE a probe sits determines which way its methylation points. A panel
# that hides structure hides the reasoning.
#
# So the gene model is now data, fetched from a real annotation, with the same
# rule the polarity table follows: every row states where it came from.
#
# FOUR SOURCES, ONE SHAPE. Whichever you use, you get the same table back:
#
#   ensembl  Ensembl REST, over the network. No install, no annotation package,
#            works for any species Ensembl carries. Cached to disk.
#   gff      A GTF or GFF3 file you already have. Same parser as `ensembl`,
#            because that is the format REST returns.
#   txdb     A TxDb object (e.g. TxDb.Hsapiens.UCSC.hg38.knownGene) plus an
#            OrgDb to turn a symbol into a gene id.
#   ensdb    An EnsDb object (e.g. EnsDb.Hsapiens.v86), which resolves symbols
#            itself.
#
# The Bioconductor sources are Suggests, gated behind requireNamespace(), and
# the package neither imports nor requires them: `ensembl` and `gff` need
# nothing beyond base R. That is deliberate. A figure this central should not
# be reachable only by users who can install a gigabyte of annotation.
#
# WHAT "CANONICAL" MEANS HERE. A gene has many transcripts and they disagree
# about where exons start. Drawing all of them is honest but unreadable at
# gene scale; drawing one is readable but a choice. The default draws the
# annotation's own canonical transcript and SAYS SO in the panel. `which =
# "all"` stacks every transcript when the isoform structure is the point.
# ============================================================================

.GM_COLS <- c("gene", "gene_id", "chr", "gene_start", "gene_end", "strand",
              "transcript", "transcript_biotype", "canonical", "feature",
              "start", "end", "exon_rank", "source", "genome")

## The feature vocabulary is closed. "exon" means an exon whose coding status
## could not be resolved - not an error, but the panel draws it at uniform
## height rather than inventing a UTR boundary it does not know.
.GM_FEATURES <- c("cds", "utr5", "utr3", "exon")

.gm_empty <- function() {
  d <- data.frame(matrix(character(0), nrow = 0, ncol = length(.GM_COLS)),
                  stringsAsFactors = FALSE)
  names(d) <- .GM_COLS
  for (v in c("gene_start", "gene_end", "start", "end", "exon_rank"))
    d[[v]] <- numeric(0)
  d$strand <- character(0)
  d$canonical <- logical(0)
  d
}

## ---------------------------------------------------------------------------
## GFF3 / GTF parsing, base R only.
##
## The two formats differ only in how the ninth column is written:
##   GFF3   key=value;key=value
##   GTF    key "value"; key "value";
## so one attribute reader handles both and the caller never has to care.
## ---------------------------------------------------------------------------

.gm_attr <- function(x, key) {
  x <- as.character(x)
  ## GFF3 first: key=value up to the next ';'
  pat3 <- paste0("(^|;)\\s*", key, "=([^;]*)")
  m <- regmatches(x, regexpr(pat3, x, perl = TRUE))
  out <- rep(NA_character_, length(x))
  hit <- vapply(regexpr(pat3, x, perl = TRUE), function(z) z > 0, TRUE)
  if (any(hit)) out[hit] <- sub(paste0("^.*?", key, "="), "", m, perl = TRUE)
  ## then GTF: key "value"
  miss <- is.na(out)
  if (any(miss)) {
    pat2 <- paste0("(^|;)\\s*", key, "\\s+\"([^\"]*)\"")
    mm <- regmatches(x[miss], regexpr(pat2, x[miss], perl = TRUE))
    h2 <- vapply(regexpr(pat2, x[miss], perl = TRUE), function(z) z > 0, TRUE)
    v <- rep(NA_character_, sum(miss))
    if (any(h2)) v[h2] <- gsub(paste0("^.*?", key, "\\s+\"|\"$"), "", mm,
                               perl = TRUE)
    out[miss] <- v
  }
  ## Ensembl prefixes ids with their type: "transcript:ENST...", "gene:ENSG..."
  sub("^(transcript|gene|CDS):", "", out)
}

.gm_read_gff <- function(path_or_text, text = FALSE) {
  con <- if (text) textConnection(path_or_text) else path_or_text
  d <- utils::read.delim(con, header = FALSE, comment.char = "#",
                         quote = "", stringsAsFactors = FALSE,
                         colClasses = "character")
  if (text) close(con)
  if (!ncol(d)) return(NULL)
  if (ncol(d) < 9)
    stop("this does not look like GFF3/GTF: ", ncol(d), " columns, expected 9",
         call. = FALSE)
  names(d)[seq_len(9)] <- c("chr", "src", "type", "start", "end", "score", "strand",
                     "phase", "attr")
  d$start <- as.numeric(d$start); d$end <- as.numeric(d$end)
  d[is.finite(d$start) & is.finite(d$end), , drop = FALSE]
}

## ---------------------------------------------------------------------------
## CDS/UTR resolution.
##
## An annotation gives exons and, separately, coding segments. What a reader
## needs is each exon split into its coding and non-coding parts, with the
## non-coding part labelled 5' or 3' - which depends on strand, and is the step
## everyone gets wrong once. On the minus strand the piece with the HIGHER
## coordinate is the 5' UTR.
## ---------------------------------------------------------------------------

.gm_split_cds <- function(ex, cds, strand) {
  if (!nrow(ex)) return(ex[0, , drop = FALSE])
  if (is.null(cds) || !nrow(cds)) {
    ex$feature <- "exon"                 # coding status unknown - say so
    return(ex)
  }
  cs <- min(cds$start); ce <- max(cds$end)
  out <- list()
  for (i in seq_len(nrow(ex))) {
    e <- ex[i, , drop = FALSE]
    s <- e$start; en <- e$end
    ## the coding slice of this exon
    c1 <- max(s, cs); c2 <- min(en, ce)
    if (c2 >= c1) {
      r <- e; r$start <- c1; r$end <- c2; r$feature <- "cds"; out[[length(out) + 1L]] <- r
    }
    ## everything left of the CDS
    if (s < cs) {
      r <- e; r$start <- s; r$end <- min(en, cs - 1)
      r$feature <- if (strand == "-") "utr3" else "utr5"
      if (r$end >= r$start) out[[length(out) + 1L]] <- r
    }
    ## everything right of the CDS
    if (en > ce) {
      r <- e; r$start <- max(s, ce + 1); r$end <- en
      r$feature <- if (strand == "-") "utr5" else "utr3"
      if (r$end >= r$start) out[[length(out) + 1L]] <- r
    }
  }
  if (!length(out)) { ex$feature <- "exon"; return(ex) }
  do.call(rbind, out)
}

.gm_assemble <- function(g, ex, cds, gene, source, genome) {
  if (is.null(ex) || !nrow(ex)) return(.gm_empty())
  ex$transcript[is.na(ex$transcript)] <- "unknown"
  txs <- unique(ex$transcript)
  rows <- list()
  for (tx in txs) {
    e <- ex[ex$transcript == tx, , drop = FALSE]
    e <- e[order(e$start), , drop = FALSE]
    st <- unique(e$strand); st <- st[st %in% c("+", "-")][1]
    if (is.na(st)) st <- "+"
    cc <- if (!is.null(cds)) cds[cds$transcript %in% tx, , drop = FALSE] else NULL
    f <- .gm_split_cds(e, cc, st)
    f$strand <- st
    rows[[length(rows) + 1L]] <- f
  }
  d <- do.call(rbind, rows)
  out <- .gm_empty()[rep(NA_integer_, nrow(d)), , drop = FALSE]
  rownames(out) <- NULL
  out$gene <- gene
  out$gene_id <- if (!is.null(g$gene_id)) g$gene_id else NA_character_
  out$chr <- .gm_chr(d$chr)
  out$gene_start <- if (!is.null(g$start)) g$start else min(d$start)
  out$gene_end <- if (!is.null(g$end)) g$end else max(d$end)
  out$strand <- d$strand
  out$transcript <- d$transcript
  out$transcript_biotype <- if (!is.null(d$biotype)) d$biotype else NA_character_
  out$canonical <- d$transcript %in%
    (if (!is.null(g$canonical) && !all(is.na(g$canonical))) g$canonical
     else .gm_pick_canonical(d))
  out$feature <- d$feature
  out$start <- d$start
  out$end <- d$end
  out$exon_rank <- if (!is.null(d$rank)) as.numeric(d$rank) else NA_real_
  out$source <- source
  out$genome <- genome
  out[order(!out$canonical, out$transcript, out$start), , drop = FALSE]
}

## Without an explicit canonical flag, the longest coding span is the least
## arbitrary stand-in - and the panel labels it as a fallback, never as the
## annotation's own choice.
.gm_pick_canonical <- function(d) {
  span <- tapply(seq_len(nrow(d)), d$transcript, function(i)
    sum(d$end[i] - d$start[i] + 1))
  names(span)[which.max(span)]
}

.gm_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

## ---------------------------------------------------------------------------
## Sources
## ---------------------------------------------------------------------------

.gm_from_gff <- function(path, gene, genome = NA, text = FALSE) {
  d <- .gm_read_gff(path, text = text)
  if (is.null(d) || !nrow(d)) return(.gm_empty())
  d$transcript <- .gm_attr(d$attr, "Parent")
  na_tx <- is.na(d$transcript)
  if (any(na_tx))
    d$transcript[na_tx] <- .gm_attr(d$attr[na_tx], "transcript_id")
  d$rank <- suppressWarnings(as.numeric(.gm_attr(d$attr, "rank")))
  d$biotype <- .gm_attr(d$attr, "biotype")
  nb <- is.na(d$biotype)
  if (any(nb)) d$biotype[nb] <- .gm_attr(d$attr[nb], "transcript_biotype")
  ## a gene row, if the file carries one, fixes the full span
  gi <- which(tolower(d$type) == "gene")
  g <- list()
  if (length(gi)) {
    g$start <- min(d$start[gi]); g$end <- max(d$end[gi])
    g$gene_id <- .gm_attr(d$attr[gi[1]], "gene_id")
  }
  ex <- d[tolower(d$type) == "exon", , drop = FALSE]
  cds <- d[tolower(d$type) %in% c("cds"), , drop = FALSE]
  .gm_assemble(g, ex, cds, gene, "gff", genome)
}

## One flat field out of a one-level JSON object, with no JSON dependency.
## Only used on Ensembl lookup responses, whose relevant fields are scalars.
.gm_json_field <- function(txt, key) {
  m <- regmatches(txt, regexec(
    sprintf('"%s"[[:space:]]*:[[:space:]]*"?([^",}]+)"?', key), txt))[[1]]
  if (length(m) < 2L) NA_character_ else m[2]
}

## Ensembl REST. Two calls: the gene's span and canonical transcript, then the
## features in that span. GFF3 rather than JSON so no JSON parser is needed -
## and so the same parser covers a user's own GTF.
.gm_from_ensembl <- function(gene, species = "homo_sapiens", genome = "hg38",
                             host = "https://rest.ensembl.org", quiet = FALSE) {
  ## host = NULL means "do not go to the network": used by tests, and by anyone
  ## who wants a hard guarantee that a figure was built offline.
  if (is.null(host) || !nzchar(host)) return(.gm_empty())
  rd <- function(u) {
    con <- url(u, open = "rb")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "\n")
  }
  ## lookup/symbol serves json/xml/jsonp ONLY - requesting text/x-gff3 here
  ## was an unconditional HTTP 400 from Ensembl ("HTTP response code said
  ## error" on every machine, online or not; found 2026-08-29 when
  ## gene_models = "auto" first exercised this path live). The four fields
  ## needed are flat in the JSON, so they are read with anchored regexes
  ## and no JSON dependency; the feature call below stays GFF3, which
  ## overlap/region does support.
  loc <- tryCatch(
    rd(sprintf("%s/lookup/symbol/%s/%s?content-type=application/json", host,
               utils::URLencode(species), utils::URLencode(gene))),
    error = function(e) NULL)
  if (is.null(loc) || !nzchar(loc)) {
    if (!quiet)
      message("could not reach Ensembl for ", gene,
              " - the gene model will be omitted rather than guessed")
    return(.gm_empty())
  }
  .jv <- function(key) .gm_json_field(loc, key)
  chr <- .jv("seq_region_name")
  gs <- suppressWarnings(as.numeric(.jv("start")))
  ge <- suppressWarnings(as.numeric(.jv("end")))
  if (is.na(chr) || !is.finite(gs) || !is.finite(ge)) {
    if (!quiet)
      message("Ensembl's answer for ", gene, " carried no genomic span - ",
              "the gene model will be omitted rather than guessed")
    return(.gm_empty())
  }
  feat <- tryCatch(
    rd(sprintf("%s/overlap/region/%s/%s:%d-%d?feature=exon;feature=cds;content-type=text/x-gff3",
               host, utils::URLencode(species), chr, gs, ge)),
    error = function(e) NULL)
  if (is.null(feat) || !nzchar(feat)) return(.gm_empty())
  d <- .gm_read_gff(feat, text = TRUE)
  if (is.null(d) || !nrow(d)) return(.gm_empty())
  d$transcript <- .gm_attr(d$attr, "Parent")
  d$rank <- suppressWarnings(as.numeric(.gm_attr(d$attr, "rank")))
  d$biotype <- .gm_attr(d$attr, "biotype")
  ## The region query returns every gene overlapping the span, so keep only
  ## transcripts whose exons fall inside this gene's own bounds. Without this a
  ## neighbouring gene's exons are drawn as if they belonged here.
  keep <- d$start >= gs & d$end <= ge
  d <- d[keep, , drop = FALSE]
  g <- list(start = gs, end = ge,
            gene_id = .jv("id"),
            canonical = .jv("canonical_transcript"))
  if (!is.na(g$canonical)) g$canonical <- sub("\\.\\d+$", "", g$canonical)
  ex <- d[tolower(d$type) == "exon", , drop = FALSE]
  cds <- d[tolower(d$type) == "cds", , drop = FALSE]
  out <- .gm_assemble(g, ex, cds, gene, "ensembl", genome)
  out$chr <- .gm_chr(chr)
  out
}

.gm_from_txdb <- function(txdb, gene, orgdb = NULL, genome = NA) {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE))
    stop("source = 'txdb' needs the GenomicFeatures package:\n",
         "  BiocManager::install(\"GenomicFeatures\")", call. = FALSE)
  gid <- gene
  if (!is.null(orgdb)) {
    if (!requireNamespace("AnnotationDbi", quietly = TRUE))
      stop("mapping a symbol through an OrgDb needs AnnotationDbi", call. = FALSE)
    gid <- tryCatch(AnnotationDbi::mapIds(orgdb, keys = gene, column = "ENTREZID",
                                          keytype = "SYMBOL"),
                    error = function(e) gene)
    gid <- unname(gid[!is.na(gid)])
    if (!length(gid))
      stop("could not map the symbol '", gene, "' to a gene id in the OrgDb",
           call. = FALSE)
  }
  ex <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  cd <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  txg <- suppressWarnings(GenomicFeatures::transcriptsBy(txdb, by = "gene"))
  if (!as.character(gid[1]) %in% names(txg))
    stop("gene id '", gid[1], "' is not in this TxDb", call. = FALSE)
  txn <- as.character(txg[[as.character(gid[1])]]$tx_name)
  .gm_from_granges_list(ex[names(ex) %in% txn], cd[names(cd) %in% txn],
                        gene, gid[1], "txdb", genome)
}

.gm_from_ensdb <- function(ensdb, gene, genome = NA) {
  if (!requireNamespace("ensembldb", quietly = TRUE))
    stop("source = 'ensdb' needs the ensembldb package:\n",
         "  BiocManager::install(\"ensembldb\")", call. = FALSE)
  flt <- AnnotationFilter::GeneNameFilter(gene)
  ex <- ensembldb::exonsBy(ensdb, by = "tx", filter = flt)
  cd <- tryCatch(ensembldb::cdsBy(ensdb, by = "tx", filter = flt),
                 error = function(e) NULL)
  .gm_from_granges_list(ex, cd, gene, NA_character_, "ensdb", genome)
}

.gm_from_granges_list <- function(ex, cds, gene, gene_id, source, genome) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("this source needs GenomicRanges", call. = FALSE)
  flat <- function(gl, type) {
    if (is.null(gl) || !length(gl)) return(NULL)
    d <- as.data.frame(unlist(gl, use.names = TRUE))
    d$transcript <- rep(names(gl), lengths(gl))
    data.frame(chr = as.character(d$seqnames), type = type,
               start = d$start, end = d$end,
               strand = as.character(d$strand),
               transcript = d$transcript,
               rank = if (!is.null(d$exon_rank)) d$exon_rank else NA_real_,
               biotype = NA_character_, stringsAsFactors = FALSE)
  }
  e <- flat(ex, "exon"); c2 <- flat(cds, "cds")
  if (is.null(e)) return(.gm_empty())
  .gm_assemble(list(gene_id = gene_id), e, c2, gene, source, genome)
}

## ---------------------------------------------------------------------------
## The user-facing entry point
## ---------------------------------------------------------------------------

#' A gene's exon, intron and UTR structure on real coordinates
#'
#' Returns one row per drawable feature - coding exon slice, 5' UTR, 3' UTR, or
#' an exon whose coding status the annotation did not resolve - for every
#' transcript of a gene, with strand and genomic start/end. This is what
#' \code{\link{dmsa_plot_locus}} draws the gene from, and it is deliberately a
#' plain data frame so it can be built, inspected, corrected and cached like
#' any other reference layer in this package.
#'
#' @param gene Gene symbol, e.g. \code{"NR3C1"}.
#' @param source One of \code{"auto"}, \code{"ensembl"}, \code{"gff"},
#'   \code{"txdb"}, \code{"ensdb"}, \code{"table"}. \code{"auto"} takes a cached
#'   copy if there is one, then Ensembl.
#' @param file GTF/GFF3 path when \code{source = "gff"}.
#' @param db TxDb or EnsDb object when \code{source} is \code{"txdb"} /
#'   \code{"ensdb"}.
#' @param orgdb OrgDb for symbol lookup with \code{"txdb"} (e.g.
#'   \code{org.Hs.eg.db}).
#' @param table A data frame already in this shape, when \code{source =
#'   "table"}.
#' @param species Ensembl species name. Default \code{"homo_sapiens"}.
#' @param genome Genome label recorded on every row. Default \code{"hg38"} -
#'   the assembly the bundled Alpha coordinates use.
#' @param cache Directory for fetched models, or \code{NULL} to disable.
#'   Default is a \code{dmsa} folder under the session temp directory.
#' @param host Ensembl REST host. \code{NULL} disables the network entirely and
#'   returns an empty model, which is how you guarantee an offline build.
#' @param quiet Suppress progress messages.
#' @return A data frame with the columns in \code{.GM_COLS}, class
#'   \code{dmsa_gene_model}. Zero rows if the model could not be obtained -
#'   never a fabricated one.
#' @seealso \code{\link{dmsa_gene_model_check}}, \code{\link{dmsa_plot_locus}}
#' @examples
#' ## from a table you already have: AVP on the minus strand
#' gm <- dmsa_gene_model("AVP", source = "table", table = data.frame(
#'   gene = "AVP", gene_id = "ENSG00000101200", chr = "chr20",
#'   gene_start = 3082556, gene_end = 3084724, strand = "-",
#'   transcript = "ENST00000380293", transcript_biotype = "protein_coding",
#'   canonical = TRUE, feature = c("utr5", "cds", "cds", "utr3"),
#'   start = c(3084665, 3084555, 3082700, 3082556),
#'   end   = c(3084724, 3084664, 3082802, 3082699),
#'   exon_rank = c(1, 1, 3, 3), source = "gff", genome = "hg38"))
#' gm
#'
#' \donttest{
#' ## Needs network: one lookup plus one region query per gene.
#' gm <- dmsa_gene_model("NR3C1")
#' subset(gm, canonical)[, c("feature", "start", "end", "exon_rank")]
#' }
#' @export
dmsa_gene_model <- function(gene, source = c("auto", "ensembl", "gff", "txdb",
                                             "ensdb", "table"),
                            file = NULL, db = NULL, orgdb = NULL, table = NULL,
                            species = "homo_sapiens", genome = "hg38",
                            cache = NULL, host = "https://rest.ensembl.org",
                            quiet = FALSE) {
  source <- match.arg(source)
  if (!is.character(gene) || length(gene) != 1L || !nzchar(gene))
    stop("gene must be a single non-empty symbol", call. = FALSE)
  if (is.null(cache)) cache <- file.path(tempdir(), "dmsa_gene_models")

  if (source == "table") {
    if (is.null(table)) stop("source = 'table' needs `table`", call. = FALSE)
    out <- as.data.frame(table, stringsAsFactors = FALSE)
    miss <- setdiff(.GM_COLS, names(out))
    if (length(miss))
      stop("the table is missing column(s): ", paste(miss, collapse = ", "),
           call. = FALSE)
    return(.gm_class(out))
  }

  cf <- if (!is.null(cache))
    file.path(cache, paste0(gsub("[^A-Za-z0-9._-]", "_", gene), "__",
                            genome, ".csv")) else NULL
  if (source %in% c("auto") && !is.null(cf) && file.exists(cf)) {
    out <- utils::read.csv(cf, stringsAsFactors = FALSE,
                           colClasses = c(chr = "character",
                                          transcript = "character"))
    if (!quiet) message("gene model for ", gene, " read from cache")
    return(.gm_class(out))
  }

  out <- switch(
    source,
    gff = { if (is.null(file)) stop("source = 'gff' needs `file`", call. = FALSE)
            .gm_from_gff(file, gene, genome) },
    txdb = { if (is.null(db)) stop("source = 'txdb' needs `db`", call. = FALSE)
             .gm_from_txdb(db, gene, orgdb, genome) },
    ensdb = { if (is.null(db)) stop("source = 'ensdb' needs `db`", call. = FALSE)
              .gm_from_ensdb(db, gene, genome) },
    .gm_from_ensembl(gene, species, genome, host = host, quiet = quiet))

  if (nrow(out) && !is.null(cf)) {
    dir.create(dirname(cf), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(out, cf, row.names = FALSE)
  }
  if (!nrow(out) && !quiet)
    message("no gene model for ", gene,
            " - the locus panel will draw the probes on their true ",
            "coordinates without the exon structure, and say so")
  .gm_class(out)
}

.gm_class <- function(x) {
  x$canonical <- as.logical(x$canonical)
  for (v in c("gene_start", "gene_end", "start", "end", "exon_rank"))
    if (v %in% names(x)) x[[v]] <- suppressWarnings(as.numeric(x[[v]]))
  class(x) <- unique(c("dmsa_gene_model", "data.frame"))
  x
}

#' Check a gene model before drawing from it
#'
#' @param gm A \code{dmsa_gene_model}, or any data frame in that shape.
#' @param verbose Print the report.
#' @return Invisibly, a list with \code{ok} and one entry per check.
#' @examples
#' ## AVP sits on the minus strand, so its 5' UTR is the HIGH coordinate
#' gm <- data.frame(
#'   gene = "AVP", gene_id = "ENSG00000101200", chr = "chr20",
#'   gene_start = 3082556, gene_end = 3084724, strand = "-",
#'   transcript = "ENST00000380293", transcript_biotype = "protein_coding",
#'   canonical = TRUE, feature = c("utr5", "cds", "cds", "utr3"),
#'   start = c(3084665, 3084555, 3082700, 3082556),
#'   end   = c(3084724, 3084664, 3082802, 3082699),
#'   exon_rank = c(1, 1, 3, 3), source = "gff", genome = "hg38")
#' dmsa_gene_model_check(gm)$ok
#' ## a feature outside the gene span means a neighbour's exons survived
#' bad <- gm; bad$start[1] <- 3000000; bad$end[1] <- 3000100
#' dmsa_gene_model_check(bad, verbose = FALSE)$ok
#' @export
dmsa_gene_model_check <- function(gm, verbose = TRUE) {
  gm <- as.data.frame(gm, stringsAsFactors = FALSE)
  res <- list(); fail <- character(); warn <- character()
  add <- function(nm, ok, detail = "", level = "fail") {
    res[[nm]] <<- list(ok = ok, detail = detail, level = level)
    if (!ok) { if (level == "warn") warn <<- c(warn, nm) else fail <<- c(fail, nm) }
  }
  miss <- setdiff(.GM_COLS, names(gm))
  add("required columns", !length(miss),
      if (length(miss)) paste("missing:", paste(miss, collapse = ", ")) else "")
  if (length(miss)) {
    if (verbose) .cas_report(res, fail, warn)
    return(invisible(list(ok = FALSE, checks = res)))
  }
  add("has at least one feature", nrow(gm) > 0,
      if (nrow(gm)) paste(nrow(gm), "rows") else "empty model")
  bad <- setdiff(unique(gm$feature), .GM_FEATURES)
  add("feature vocabulary is closed", !length(bad),
      if (length(bad)) paste("unknown:", paste(bad, collapse = ", ")) else "")
  add("start <= end on every row", all(gm$end >= gm$start, na.rm = TRUE),
      paste(sum(gm$end < gm$start, na.rm = TRUE), "reversed"))
  add("one strand per transcript",
      all(tapply(gm$strand, gm$transcript, function(z) length(unique(z))) == 1),
      "")
  add("one chromosome", length(unique(gm$chr)) == 1,
      paste(unique(gm$chr), collapse = ", "))
  ## Features outside the stated gene span mean the wrong gene's exons were
  ## kept from a region query - the failure this is here to catch.
  out_of_span <- sum(gm$start < gm$gene_start | gm$end > gm$gene_end, na.rm = TRUE)
  add("every feature inside the gene span", out_of_span == 0,
      paste(out_of_span, "feature(s) outside"))
  add("a canonical transcript is marked", any(gm$canonical),
      if (any(gm$canonical)) paste(length(unique(gm$transcript[gm$canonical])),
                                   "canonical transcript(s)") else
        "none - the panel will fall back to the longest", level = "warn")
  ## Overlapping CDS and UTR on one transcript means the split went wrong.
  ov <- 0L
  for (tx in unique(gm$transcript)) {
    g <- gm[gm$transcript == tx, , drop = FALSE]
    g <- g[order(g$start), , drop = FALSE]
    if (nrow(g) > 1) ov <- ov + sum(g$start[-1] <= g$end[-nrow(g)])
  }
  add("no overlapping features within a transcript", ov == 0,
      paste(ov, "overlap(s)"), level = "warn")
  if (verbose) {
    cat(sprintf("gene model: %s | %s | %d transcript(s) | %d feature(s) | %s\n",
                gm$gene[1], gm$chr[1], length(unique(gm$transcript)), nrow(gm),
                paste0(unique(gm$source), collapse = "/")))
    .cas_report(res, fail, warn)
  }
  invisible(list(ok = !length(fail), checks = res, warnings = warn))
}

#' @export
print.dmsa_gene_model <- function(x, ...) {
  if (!nrow(x)) { cat("dmsa gene model: empty (no annotation obtained)\n")
    return(invisible(x)) }
  cn <- unique(x$transcript[x$canonical])
  cat(sprintf("dmsa gene model: %s  %s:%s-%s (%s strand)\n", x$gene[1], x$chr[1],
              format(x$gene_start[1], big.mark = ","),
              format(x$gene_end[1], big.mark = ","), x$strand[1]))
  cat(sprintf("  %d transcript(s), %d feature(s), from %s (%s)\n",
              length(unique(x$transcript)), nrow(x), x$source[1], x$genome[1]))
  ft <- table(factor(x$feature, levels = .GM_FEATURES))
  cat("  ", paste(sprintf("%s %d", names(ft)[ft > 0], ft[ft > 0]),
                  collapse = " | "), "\n", sep = "")
  if (length(cn)) cat("  canonical: ", paste(cn, collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' The locus panel drawn by Gviz
#'
#' The same probes and the same gene model as \code{\link{dmsa_plot_locus}},
#' handed to \pkg{Gviz} instead of drawn here. Use this when you want what Gviz
#' gives and this package does not: a real cytoband ideogram, transcript
#' collapsing, and the ability to stack any other \code{GdObject} you already
#' build - CpG islands, ATAC peaks, ChIP tracks, conservation - in the same
#' coordinate frame.
#'
#' Gviz, GenomicRanges and IRanges are Suggests, not dependencies. This function
#' errors with the install line if they are absent; everything else in the
#' package keeps working without them.
#'
#' @param probes Data frame with \code{probe}, \code{pos}, \code{b}, and
#'   optionally \code{d}, \code{p}, \code{chr} - the same frame
#'   \code{dmsa_plot_locus()} takes.
#' @param gene Gene symbol, used for the title and the model lookup.
#' @param gene_model A \code{dmsa_gene_model}, or \code{TRUE} to fetch one.
#' @param chrom Chromosome, e.g. \code{"chr5"}. Taken from the model or the
#'   probes when absent.
#' @param genome Genome label passed to Gviz, default \code{"hg38"}.
#' @param extra A list of further Gviz tracks to stack under the model.
#' @param ideogram Draw the cytoband ideogram. Needs network on first use, as
#'   Gviz fetches the band table from UCSC.
#' @param file Optional PNG path.
#' @param width,height,res PNG geometry.
#' @return Invisibly, the list of tracks plotted.
#' @seealso \code{\link{dmsa_plot_locus}} for the dependency-free panel.
#' @examples
#' \donttest{
#' ## Needs the Gviz package; the gene model comes from a table, offline.
#' if (requireNamespace("Gviz", quietly = TRUE)) {
#'   gm <- dmsa_gene_model("AVP", source = "table", table = data.frame(
#'     gene = "AVP", gene_id = "ENSG00000101200", chr = "chr20",
#'     gene_start = 3082556, gene_end = 3084724, strand = "-",
#'     transcript = "ENST00000380293", transcript_biotype = "protein_coding",
#'     canonical = TRUE, feature = c("utr5", "cds", "cds", "utr3"),
#'     start = c(3084665, 3084555, 3082700, 3082556),
#'     end   = c(3084724, 3084664, 3082802, 3082699),
#'     exon_rank = c(1, 1, 3, 3), source = "gff", genome = "hg38"))
#'   probes <- data.frame(probe = paste0("cg", 1:5), chr = "chr20",
#'                        pos = c(3082600, 3082750, 3083050, 3084600, 3084700),
#'                        b = c(-.04, .02, .03, -.05, -.06), se = rep(.015, 5),
#'                        d = c(-1, 1, 1, -1, -1), p = c(.01, .3, .04, .002, .001))
#'   dmsa_plot_locus_gviz(probes, gene = "AVP", gene_model = gm,
#'                        file = tempfile(fileext = ".png"))
#' }
#' }
#' @export
dmsa_plot_locus_gviz <- function(probes, gene = "", gene_model = NULL,
                                 chrom = NA, genome = "hg38", extra = list(),
                                 ideogram = FALSE, file = NULL, width = 1800,
                                 height = 1100, res = 190) {
  need <- c("Gviz", "GenomicRanges", "IRanges")
  absent <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(absent))
    stop("this engine needs ", paste(absent, collapse = ", "), ":\n",
         "  BiocManager::install(c(",
         paste0("\"", absent, "\"", collapse = ", "), "))\n",
         "dmsa_plot_locus() draws the same gene model with no dependencies.",
         call. = FALSE)

  p <- as.data.frame(probes, stringsAsFactors = FALSE)
  for (v in c("probe", "pos", "b"))
    if (!v %in% names(p)) stop("probes needs a '", v, "' column", call. = FALSE)
  p$pos <- suppressWarnings(as.numeric(p$pos))
  p <- p[is.finite(p$pos) & is.finite(suppressWarnings(as.numeric(p$b))), ,
         drop = FALSE]
  if (!nrow(p)) stop("no probe has both a position and an effect", call. = FALSE)

  if (isTRUE(gene_model))
    gene_model <- tryCatch(dmsa_gene_model(gene, genome = genome, quiet = TRUE),
                           error = function(e) NULL)
  if (is.na(chrom)) {
    chrom <- if (!is.null(gene_model) && nrow(gene_model)) gene_model$chr[1]
    else if (!is.null(p$chr)) .gm_chr(unique(as.character(p$chr))[1]) else NA
  }
  if (is.na(chrom)) stop("cannot tell which chromosome this is; pass `chrom`",
                         call. = FALSE)

  trk <- list()
  if (ideogram)
    trk$ideo <- try(Gviz::IdeogramTrack(genome = genome, chromosome = chrom),
                    silent = TRUE)
  if (inherits(trk$ideo, "try-error")) trk$ideo <- NULL
  trk$axis <- Gviz::GenomeAxisTrack(add53 = TRUE, add35 = TRUE,
                                    labelPos = "below")

  ## The probes, signed the way DMSA reads them: the bar is the effect, so a
  ## reader sees direction and magnitude at the position it happened.
  dat <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(start = p$pos, width = 1),
    effect = as.numeric(p$b))
  trk$data <- Gviz::DataTrack(dat, data = "effect", type = c("h", "p"),
                              name = "effect", baseline = 0,
                              col.baseline = "grey40", cex = .8, lwd = 2)
  if (!is.null(gene_model) && nrow(gene_model)) {
    gr <- .gm_to_granges(gene_model)
    trk$gene <- Gviz::GeneRegionTrack(
      gr, genome = genome, chromosome = chrom,
      name = if (nzchar(gene)) gene else "gene",
      transcriptAnnotation = "transcript", collapseTranscripts = FALSE,
      fill = "grey40", col = "grey30")
  }
  if (length(extra)) trk <- c(trk, extra)
  trk <- trk[!vapply(trk, is.null, TRUE)]

  if (!is.null(file)) {
    grDevices::png(file, width, height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  from <- min(c(p$pos, if (!is.null(gene_model)) gene_model$start))
  to   <- max(c(p$pos, if (!is.null(gene_model)) gene_model$end))
  pad <- max((to - from) * .06, 30)
  Gviz::plotTracks(unname(trk), from = from - pad, to = to + pad,
                   chromosome = chrom,
                   main = if (nzchar(gene)) gene else NULL, cex.main = 1)
  invisible(trk)
}

## GRanges handoff for the Gviz engine. Kept here so the drawing code never
## touches Bioconductor classes directly.
.gm_to_granges <- function(gm) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("this needs GenomicRanges", call. = FALSE)
  keep <- gm$feature %in% c("cds", "utr5", "utr3", "exon")
  g <- gm[keep, , drop = FALSE]
  GenomicRanges::GRanges(
    seqnames = g$chr,
    ranges = IRanges::IRanges(start = g$start, end = g$end),
    strand = g$strand,
    feature = ifelse(g$feature == "cds", "protein_coding",
                     ifelse(g$feature == "exon", "unknown", g$feature)),
    gene = g$gene, exon = paste0(g$transcript, ":", seq_len(nrow(g))),
    transcript = g$transcript, symbol = g$gene)
}
