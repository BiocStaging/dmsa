## ===========================================================================
## THE BUNDLED DIRECTION MAPS
##
## DMSA cannot align anything without knowing, per probe, whether higher
## methylation raises or lowers the expression of the gene it belongs to.
## Requiring every user to obtain that map from a separate package before
## they can run anything is the single largest barrier between "I have a
## normalised matrix" and "I have a DMSA result", so usable maps ship here.
##
## WHAT IS BUNDLED, AND WHAT IS NOT. Two layers of cpgdirection ship: blood
## and nasal epithelium, each restricted to the rows that carry an actual
## call. About two in three blood rows and half the nasal rows of the full
## tables are ABSTAIN - the model declined to commit - and an abstention is not
## a direction, so it cannot be aligned to and is not shipped. That restriction
## keeps the two maps at about 1.7 MB each with no loss of anything DMSA could
## have used. inst/scripts/make_direction_maps.R builds them.
##
## WHY TWO LAYERS AND NOT ONE. Direction is tissue-specific, and the array
## studies that will use this package are not one tissue. Whole blood is the
## externally validated layer and the right proxy for an immune-dominated
## sample: adult saliva by passive drool runs about 82% immune, so blood calls
## are the correct calls there. Neonatal and infant buccal swabs are the
## opposite - .89 to 1.00 epithelial - and applying immune-derived calls to a
## near-pure epithelial sample would align every probe against the wrong
## tissue while reporting full coverage. Shipping blood alone would have made
## that silent. Nasal epithelium is the epithelial proxy.
##
## The brain bridge, the SMR causal layer, the remaining tissues and the full
## 3.3M-pair table stay in cpgdirection, which remains the complete resource.
##
## VERSIONING. Annotation drift silently changes enrichment results, and a
## direction map is an annotation. A DMSA finding that moved because the map
## moved would be the same failure this package exists to name. So each map
## carries a version, and every alignment records which version produced it.
## ===========================================================================

.dm_cache <- new.env(parent = emptyenv())

.dm_tissues <- c("blood", "epithelium")

.direction_map <- function(tissue = "blood") {
  tissue <- match.arg(tissue, .dm_tissues)
  key <- paste0("map_", tissue)
  if (!is.null(.dm_cache[[key]])) return(.dm_cache[[key]])
  f <- system.file("extdata", sprintf("direction_%s_hg19.rds", tissue),
                   package = "dmsa")
  if (!nzchar(f))
    stop("the bundled ", tissue, " direction map is missing from the ",
         "installed package; reinstall dmsa, or pass your own direction ",
         "table", call. = FALSE)
  m <- readRDS(f)
  .dm_cache[[key]] <- m
  m
}

#' The direction calls DMSA ships with
#'
#' Returns per-probe methylation-to-expression direction calls from the maps
#' bundled with dmsa, so an alignment can be built without obtaining a
#' direction resource separately.
#'
#' Two layers of \code{cpgdirection} are bundled, each restricted to the
#' probe-gene pairs that carry an actual call. About two in three blood rows
#' and half the nasal rows of the full tables abstain, and an abstention cannot
#' be aligned to, so those rows are not shipped - which is why the reduction
#' costs nothing DMSA could have used.
#'
#' \describe{
#'   \item{\code{"blood"}}{Whole blood: 720,029 pairs over 314,314 probes and
#'     6,811 genes. The externally validated layer, and the right calls for an
#'     immune-dominated sample - adult saliva by passive drool runs about 82%
#'     immune.}
#'   \item{\code{"epithelium"}}{Nasal epithelium: 625,937 pairs over 270,449
#'     probes and 4,574 genes. The right calls for an epithelium-dominated
#'     sample - neonatal and infant buccal swabs run .89 to 1.00 epithelial,
#'     where blood calls would be the wrong tissue.}
#' }
#'
#' Choose the layer by what the sample is made of, not by what it is called.
#' The failure mode this guards against is quiet: an epithelial sample scored
#' against immune-derived calls reports full coverage while every direction is
#' drawn from the wrong tissue.
#'
#' The brain bridge, the SMR causal layer, the remaining tissues and the full
#' 3.3M-pair table live in \code{cpgdirection}
#' (\doi{10.5281/zenodo.22024185}), which remains the complete resource.
#'
#' Map version 2026.11c is built from the cpgdirectionData 0.99.7 lookups, in
#' which every gene name resolves to one Ensembl gene, at a strand-aware hg19
#' TSS, and the tissue models are retrained on harmonized pools. It replaces
#' 2026.11 (the 0.99.5 lookups), whose training pools held wrong-end
#' distances for Framingham's minus-strand genes and counted one gene's
#' several names as several genes; 2026.11 in turn replaced 2026.08, whose
#' lookups paired 84,117 blood and 15,313 epithelium rows with a gene on
#' another chromosome.
#'
#' A probe may map to more than one gene - a median of two, up to seven - and
#' every mapped pair is returned, because which gene a probe is read against
#' is a modelling decision rather than a lookup.
#'
#' @param probes Character vector of CpG identifiers. \code{NULL} (default)
#'   returns the whole map.
#' @param genes Optional character vector restricting the result. Either one
#'   gene per entry of \code{probes}, selecting those specific pairs, or a
#'   shorter vector naming genes to keep.
#' @param tissue Which bundled layer: \code{"blood"} (default) or
#'   \code{"epithelium"}. See Details for how to choose.
#' @return A \code{data.frame} with \code{probe}, \code{gene}, \code{d}
#'   (\code{-1} or \code{+1}), \code{p_plus} (the probability the direction is
#'   \code{+1}) and \code{tier} (evidence grade). The map's provenance and
#'   version are attached as the \code{"dmsa_direction_map"} attribute.
#' @examples
#' # What the bundled blood map knows about one gene's probes.
#' avp <- dmsa_directions(genes = "AVP")
#' nrow(avp)
#' head(avp, 3)
#'
#' # Coverage of a probe set you care about.
#' probes <- c("cg00000029", "cg00000165", "not_a_real_probe")
#' got <- dmsa_directions(probes)
#' got
#'
#' # The layers are not nested: each calls pairs the other declines to call,
#' # so coverage of the same probe set differs by tissue.
#' dmsa_directions(probes, tissue = "epithelium")
#'
#' # Provenance travels with the map.
#' attr(dmsa_directions(), "dmsa_direction_map")[c("tissue", "version")]
#' @export
dmsa_directions <- function(probes = NULL, genes = NULL, tissue = "blood") {
  m <- .direction_map(tissue)
  meta <- attr(m, "dmsa_direction_map")

  if (!is.null(probes)) {
    probes <- as.character(probes)
    keep <- m$probe %in% probes
    m <- m[keep, , drop = FALSE]
  }
  if (!is.null(genes)) {
    genes <- as.character(genes)
    if (!is.null(probes) && length(genes) == length(probes)) {
      want <- paste(probes, genes, sep = "\r")
      m <- m[paste(as.character(m$probe), as.character(m$gene),
                   sep = "\r") %in% want, , drop = FALSE]
    } else {
      m <- m[as.character(m$gene) %in% genes, , drop = FALSE]
    }
  }
  m$probe <- as.character(m$probe)
  m$gene  <- as.character(m$gene)
  m$tier  <- as.character(m$tier)
  rownames(m) <- NULL
  attr(m, "dmsa_direction_map") <- meta
  m
}
