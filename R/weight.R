# ============================================================================
# RELIABILITY (CENTRALITY) WEIGHTS - the default aligned weighting.
#
# DMSA ships two engines. Both align each probe to its expression consequence;
# they differ only in how much each aligned probe counts toward its unit's
# score:
#
#   "reliability" (DEFAULT) weights probe j by its item-REST correlation - how
#       much it tracks the rest of its own unit - so a probe that reliably
#       reports the unit's shared axis counts more than a noisy one. This is
#       classical reliability weighting (corrected item-total correlation), and
#       in simulation it GAINS power on units whose probe quality is
#       heterogeneous while LOSING nothing when it is uniform (unit weights are
#       already near-optimal there). It is the default because it weakly
#       dominates the flat engine.
#
#   "flat" weights every usable aligned probe equally (the original engine).
#       Retained as an option and as the exact reproducibility reference.
#
# THE CALIBRATION GUARANTEE. The weights are a function of the methylation
# matrix ALONE - never the outcome - so under Freedman-Lane permutation of the
# response residuals the weights are frozen. This is permutation-invariant when
# the resampled response is external to methylation, and simulations cover that
# use case. When methylation CHANGE itself is the outcome, the weights are
# outcome-adaptive; dmsa_change() therefore refuses non-flat weighting on
# that arm until dedicated null calibration is available.
#
# THREE GUARDS keep the weights honest:
#   (1) a unit with < 3 usable direction-called probes falls back to flat -
#       the item-rest correlation is too unstable to estimate;
#   (2) a unit whose aligned probes carry no coherent axis - the leading
#       eigenvalue of their correlation matrix is below `w_floor` (default 1.5;
#       independent probes give eigenvalues near 1, so this asks for a shared
#       factor stronger than sampling noise) - falls back to flat; there is
#       nothing to be central to, and reweighting would chase noise;
#   (3) a probe anti-correlated with its unit's axis clamps to weight 0, and if
#       every probe clamps to 0 the unit falls back to flat.
# Surviving weights are renormalised to mean 1 within the unit so the aligned
# score keeps the same scale as the flat engine.
# ============================================================================

## null-coalescing helper (internal; used across the report layer)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Outcome-free reliability weights for aligned probes
#'
#' Item-rest correlation of each aligned probe with the rest of its unit,
#' computed from the (standardised, already winsorised) methylation matrix only.
#' Used internally by \code{\link{dmsa_triangulate}}, \code{\link{dmsa_scores}}
#' and \code{\link{dmsa_change}} when \code{weighting = "reliability"}.
#'
#' @param Z Numeric matrix, samples x probes, standardised methylation (the same
#'   matrix the lenses see; outcome-free).
#' @param units Character/label vector, one per column of \code{Z}.
#' @param mlt Numeric per-probe aligned multiplier (0 for unusable probes).
#' @param w_floor Minimum leading eigenvalue of the aligned probes' correlation
#'   matrix for weighting to engage (independent probes give eigenvalues near 1;
#'   the default 1.5 asks for a shared factor stronger than sampling noise).
#'   Below it the unit is weighted flat.
#' @param weighting \code{"reliability"} (default) or \code{"flat"} (returns all
#'   ones).
#' @return Numeric vector of per-probe weights, length \code{ncol(Z)}, mean 1
#'   within each weighted unit and 1 everywhere a guard fired.
#' @examples
#' set.seed(2)
#' n <- 120
#' f <- rnorm(n)
#' ## three probes report the gene's shared axis well, three barely at all
#' lam <- c(.8, .8, .8, .1, .1, .1)
#' Z <- scale(sapply(lam, function(l) l * f + sqrt(1 - l^2) * rnorm(n)))
#'
#' ## weights come from methylation alone, so the permutation null is untouched
#' round(dmsa_relweights(Z, units = rep("NR3C1", 6), mlt = rep(1, 6)), 2)
#' @export
dmsa_relweights <- function(Z, units, mlt, w_floor = 1.5,
                            weighting = c("reliability", "flat")) {
  weighting <- match.arg(weighting)
  P <- ncol(Z); relw <- rep(1, P)
  if (weighting == "flat") return(relw)
  for (u in unique(units)) {
    j <- which(units == u & is.finite(mlt) & mlt != 0)
    if (length(j) < 3L) next                              # guard (1): small K
    A <- sweep(Z[, j, drop = FALSE], 2, sign(mlt[j]), "*")  # orient to +axis
    C <- suppressWarnings(stats::cor(A))
    if (any(!is.finite(C))) next
    ev1 <- eigen(C, symmetric = TRUE, only.values = TRUE)$values[1]
    if (!is.finite(ev1) || ev1 < w_floor) next             # guard (2): no axis
    r <- vapply(seq_along(j), function(k) {
      rr <- suppressWarnings(stats::cor(A[, k], rowMeans(A[, -k, drop = FALSE])))
      if (is.finite(rr)) rr else 0 }, numeric(1))
    w <- pmax(0, r)                                        # guard (3): clamp
    if (sum(w > 0) < 2L) next                              # need >=2 for an axis
    relw[j] <- w / mean(w[w > 0])                          # mean 1 within unit
  }
  relw
}
## Shared resolution of row-grouping arguments (block, ri_group, id).
##
## Every permutation engine here takes "labels, one per row". Passing a COLUMN
## NAME instead is the natural mistake - and it did not error: a length-1
## string recycles into one block holding every row, the only within-block
## permutation is the identity, every permuted statistic equals the observed
## one, and p = 1 for eternity. Silent total power loss. So names are now
## resolved, vectors validated, and the degenerate one-block case refused.
.dmsa_rows <- function(x, data, what) {
  if (is.null(x)) return(NULL)
  n <- nrow(data)
  if (is.character(x) && length(x) != n && all(x %in% names(data))) {
    x <- if (length(x) == 1L) data[[x]]
         else interaction(data[, x, drop = FALSE], drop = TRUE)
  } else if (length(x) != n) {
    stop(what, " must be a column name of `data` or one label per row ",
         "(got length ", length(x), " for ", n, " rows",
         if (is.character(x) && length(x) <= 3L)
           paste0("; \"", paste(x, collapse = "\", \""),
                  "\" is not a column of data") else "", ")", call. = FALSE)
  }
  x
}

.dmsa_check_block <- function(blk, what = "block") {
  if (is.null(blk)) return(invisible(NULL))
  u <- unique(blk)
  if (length(u) == 1L && length(blk) > 1L)
    stop("every row falls in ONE ", what, " - no permutation is possible ",
         "and every p-value would be 1. This usually means a column name ",
         "was recycled or the wrong column was chosen.", call. = FALSE)
  invisible(NULL)
}

## E9 (2026-08-29, PI-approved): guard the positional join between a matrix's
## columns and an alignment's rows. The engines apply multiplier k to column k
## BY POSITION; when both sides carry names, that assumption is checked on the
## canonical CpG id (the "cg..." prefix), because matrix columns are often the
## panel's suffixed design names while the alignment carries bare probe ids.
## Outcomes: "ok" (orders agree, proceed), a reorder permutation (same unique
## set, different order - fixable unambiguously), a hard error (overlapping
## but different/misordered sets - the silent-scramble case this exists to
## kill), or "positional" (no overlap even canonically: different naming
## conventions; the caller's order is trusted, exactly as before this guard).
.dmsa_align_order <- function(m_names, al_probes, what = "alignment") {
  if (is.null(m_names) || is.null(al_probes) ||
      length(m_names) != length(al_probes)) return("positional")
  canon <- function(x) {
    x <- as.character(x)
    ifelse(grepl("^cg[0-9]+", x), sub("^(cg[0-9]+).*$", "\\1", x), x)
  }
  cm <- canon(m_names); ca <- canon(al_probes)
  if (identical(cm, ca)) return("ok")
  if (!length(intersect(cm, ca))) return("positional")
  if (setequal(cm, ca) && !anyDuplicated(cm) && !anyDuplicated(ca))
    return(match(cm, ca))
  bad <- utils::head(unique(c(setdiff(cm, ca), setdiff(ca, cm),
                              cm[cm != ca])), 5)
  stop("the ", what, "'s probes and the matrix columns overlap but do not ",
       "line up (first differences: ", paste(bad, collapse = ", "), ").\n",
       "Multipliers are applied BY POSITION, so proceeding would put ",
       "direction calls on the wrong probes. Reorder one side so they match ",
       "(dmsa_align() returns rows in the order of the probes you pass it).",
       call. = FALSE)
}
