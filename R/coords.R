
## ---------------------------------------------------------------------------
## OPTIONAL COMPANION PACKAGE
##
## `cpgdirection` supplies the direction calls and a CpG coordinate manifest.
## It is genuinely optional: every path below degrades or names an alternative
## when it is absent, and dmsa declares no dependency on it, because it is
## distributed from GitHub and Zenodo rather than from CRAN or Bioconductor
## (Bioconductor forbids declaring dependencies that its build system cannot
## install). Resolving it by name at call time rather than with `::` keeps that
## true statement true for R CMD check as well as for the user.
## ---------------------------------------------------------------------------
.cpgd <- function(fun = NULL) {
  pkg <- "cpgdirection"
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
  if (is.null(fun)) return(TRUE)
  if (!fun %in% getNamespaceExports(pkg)) return(NULL)
  get(fun, envir = asNamespace(pkg))
}
# ============================================================================
# PROBE COORDINATES WITHOUT BIOCONDUCTOR
#
# The previous route was minfi plus IlluminaHumanMethylationEPICanno, which is
# roughly a gigabyte of annotation packages and a Bioconductor install, to
# obtain two numbers per probe. Four cheaper routes, in the order you should
# try them:
#
#   0. dmsa_probe_coords(probes)      no arguments beyond the probe ids. Uses
#      cpgdirection::cpgd_cpg_positions(), 930k probes, already on your machine
#      because cpgdirection is where the direction calls `d` come from. Nothing
#      to install, nothing to download.
#   1. dmsa_probe_coords(file = ...)  read a plain-text manifest you already
#      have. No download, no Bioconductor, works offline.
#   2. dmsa_probe_coords(url  = ...)  read one from a public plain-text
#      manifest. utils::download.file and utils::read.delim only - both base R.
#   3. NOTHING AT ALL. dmsa_plot_locus() works with no coordinates: it orders
#      the probes by their RefGene region if you have one, otherwise keeps the
#      order you gave, spaces them evenly, and says so on the axis. The
#      colours, the direction split and the effect panel are unaffected,
#      because none of them ever needed a genomic position.
#
# Whichever route, run it ONCE and write the result next to your project. A
# three-column csv of probe, chromosome and position is a project asset, not a
# dependency.
# ============================================================================

#' Probe genomic coordinates, without Bioconductor
#'
#' Called with probe ids alone, reads them from
#' \code{cpgdirection::cpgd_cpg_positions()} - the same package the
#' \code{cpgdirection} calls \code{d} come from, so nothing new is installed.
#' Given \code{file} or \code{url} instead, reads any delimited manifest with
#' base R. No annotation packages either way.
#'
#' @param probes Character vector of probe ids to look up.
#' @param file Path to a local manifest (plain text or \code{.gz}). Takes
#'   precedence over \code{url}.
#' @param url URL of a plain-text manifest to download once. Public Infinium
#'   manifests in tsv form are published by the Zhou lab
#'   (\url{https://zwdzwd.github.io/InfiniumAnnotation}); any file with a probe
#'   column and chromosome and position columns will do.
#' @param cache Directory to keep the downloaded manifest in, so the download
#'   happens once. Defaults to \code{tools::R_user_dir("dmsa", "cache")}.
#' @param probe_col,chr_col,pos_col Column names in the manifest. Defaults
#'   cover the common spellings; if they are not found the function lists the
#'   columns it did find rather than guessing.
#' @return data.frame with \code{probe}, \code{chr}, \code{pos}, restricted to
#'   \code{probes} and in that order. Probes not found come back with
#'   \code{NA}, which \code{dmsa_plot_locus()} handles.
#' @seealso \code{\link{dmsa_plot_locus}}, which needs none of this.
#' @examples
#' ## from a manifest you already have (a three-probe stand-in here)
#' f <- tempfile(fileext = ".tsv")
#' utils::write.table(data.frame(probeID = c("cg00000029", "cg00000108",
#'                                           "cg00000109"),
#'                               CpG_chrm = c("chr16", "chr3", "chr3"),
#'                               CpG_beg = c(53434200L, 37417715L, 172198247L)),
#'                    f, sep = "\t", quote = FALSE, row.names = FALSE)
#' dmsa_probe_coords(c("cg00000029", "cg00000109"), file = f)
#'
#' \donttest{
#' ## the full EPIC v2 positions, through cpgdirection when it is installed,
#' ## or a published manifest downloaded once into a cache
#' ## co <- dmsa_probe_coords(my_probes)
#' ## co <- dmsa_probe_coords(my_probes, url = "<manifest URL>")
#' }
#' @export
dmsa_probe_coords <- function(probes, file = NULL, url = NULL, cache = NULL,
                              probe_col = c("probeID", "Probe_ID", "probe",
                                            "IlmnID", "Name", "cpg", "cpg_id"),
                              chr_col = c("CpG_chrm", "chrm", "chr", "CHR",
                                          "seqnames", "chromosome"),
                              pos_col = c("CpG_beg", "start", "pos", "MAPINFO",
                                          "position", "CpG_end")) {
  probes <- as.character(probes)

  ## ---- route 0: cpgdirection, already installed ---------------------------
  if (is.null(file) && is.null(url)) {
    if (is.null(.cpgd()))
      stop("no coordinate source. Either install cpgdirection (which you ",
           "already need for the direction calls) from ",
           "https://github.com/teindor/cpgdirection or ",
           "https://doi.org/10.5281/zenodo.22024185, or pass file= / url= ",
           "with a manifest. Or skip coordinates entirely: dmsa_plot_locus() ",
           "works without them.", call. = FALSE)
    .pos <- .cpgd("cpgd_cpg_positions")
    if (is.null(.pos))
      stop("this cpgdirection has no cpgd_cpg_positions(); upgrade it, or ",
           "pass file= / url=.", call. = FALSE)
    man <- as.data.frame(
      .pos(), stringsAsFactors = FALSE)
    return(.coords_match(man, probes, probe_col, chr_col, pos_col,
                         "cpgdirection"))
  }

  ## ---- routes 1 and 2: a plain-text manifest ------------------------------
  if (is.null(file)) {
    if (is.null(cache))
      cache <- tryCatch(tools::R_user_dir("dmsa", "cache"),
                        error = function(e) tempdir())
    dir.create(cache, showWarnings = FALSE, recursive = TRUE)
    file <- file.path(cache, basename(url))
    if (!file.exists(file)) {
      message("downloading manifest once to ", file)
      utils::download.file(url, file, mode = "wb", quiet = FALSE)
    } else message("using cached manifest ", file)
  }
  if (!file.exists(file)) stop("manifest not found: ", file, call. = FALSE)
  con <- if (grepl("\\.gz$", file)) gzfile(file) else file
  man <- utils::read.delim(con, stringsAsFactors = FALSE, check.names = FALSE,
                           comment.char = "#")
  .coords_match(man, probes, probe_col, chr_col, pos_col, basename(file))
}

## Match probe/chr/pos columns in `man` and pull out `probes`, in that order.
.coords_match <- function(man, probes, probe_col, chr_col, pos_col, what) {
  pick <- function(cands, lab) {
    hit <- intersect(cands, names(man))
    if (!length(hit))
      stop("no ", lab, " column found in ", what, ". It has: ",
           paste(utils::head(names(man), 25), collapse = ", "),
           ". Pass the right name explicitly.", call. = FALSE)
    hit[1]
  }
  pc <- pick(probe_col, "probe"); cc <- pick(chr_col, "chromosome")
  po <- pick(pos_col, "position")
  i <- match(probes, as.character(man[[pc]]))
  .pos <- suppressWarnings(as.numeric(man[[po]][i]))
  ## E10 (PI-approved): the Zhou-lab *_beg columns follow the BED convention
  ## (0-based), while everything drawn beside them (Ensembl gene models,
  ## MAPINFO) is 1-based. Convert, so a probe never sits 1 bp left of its
  ## own exon in a locus figure. A generic `start` column is NOT shifted:
  ## R-side exports (GRanges and friends) use 1-based starts.
  if (grepl("_beg$", po)) .pos <- .pos + 1
  out <- data.frame(probe = probes,
                    chr = sub("^chr", "", as.character(man[[cc]])[i]),
                    pos = .pos,
                    stringsAsFactors = FALSE)
  n_ok <- sum(is.finite(out$pos))
  absent <- if (n_ok < length(probes))
    paste0(" (", length(probes) - n_ok, " absent - dmsa_plot_locus() ",
           "then spaces every probe evenly rather than dropping them)")
  else ""
  message("matched ", n_ok, " of ", length(probes), " probes from ", what,
          absent)
  out
}

#' Save a coordinate table so the lookup never has to be repeated
#'
#' @param coords Result of \code{dmsa_probe_coords()}.
#' @param path csv path.
#' @return \code{path}, invisibly.
#' @examples
#' ## run the lookup once, then keep the csv with the project
#' co <- data.frame(probe = c("cg00000029", "cg00000108", "cg00000109"),
#'                  chr = c("16", "3", "3"),
#'                  pos = c(53434200, 37417716, 171916037))
#'
#' f <- dmsa_write_coords(co, tempfile(fileext = ".csv"))
#' utils::read.csv(f)
#' @export
dmsa_write_coords <- function(coords, path = "probe_coords.csv") {
  utils::write.csv(coords, path, row.names = FALSE)
  message("wrote ", path, " - keep it with the project; the lookup is no ",
          "longer needed")
  invisible(path)
}
