## ===========================================================================
## THE P-DANCE TEST
##
## A finding is a claim about a unit. If the unit's p-value moves when the
## analyst changes which OTHER probes are in the analysis - probes that carry
## none of the unit's evidence - then the number was never a property of the
## finding alone, and no two laboratories will reproduce it, because no two
## laboratories load the same probe set. Detection-p QC, artifact masks and
## the 450K/EPIC transition put cross-cohort divergence at a few per cent to
## a few tens of per cent by construction.
##
## The test perturbs the analysis set in one of two directions while holding
## the focal unit's own probes fixed, and reports what each engine's p-value
## does. It is engine-agnostic on purpose: an engine here is any function
## from a set of probes to a p-value, so a competitor is tested by exactly
## the procedure that tests DMSA.
## ===========================================================================

#' The p-dance test: does a finding's p-value survive a change of analysis set?
#'
#' Perturbs the analysis set around a fixed finding and records what happens to
#' each engine's p-value. The focal unit's own probes are never touched, so any
#' movement is attributable to the analysis set rather than to the evidence.
#'
#' Two forms, both reported in the same way:
#'
#' \describe{
#'   \item{\code{"dropout"}}{Removes a random fraction of the non-focal probes.
#'     The omitted probes carry no focal evidence, so a finding that depends on
#'     them was never about its own unit.}
#'   \item{\code{"addition"}}{Adds probes drawn from \code{pool}, which sits
#'     outside the analysis set. This is the direction a competitive engine is
#'     most sensitive to, because added probes enter its background.}
#' }
#'
#' The two forms move competitive and whole-set engines in opposite directions,
#' which is why running only one of them understates the problem.
#'
#' @param engines Named list of functions. Each takes one argument - the
#'   members of a perturbed analysis set, in whatever form \code{set} is given
#'   - and returns a single p-value. Anything expressible this way can be
#'   tested, which is the point: DMSA and its competitors are put through an
#'   identical procedure.
#' @param set The full analysis set: a character vector of probe identifiers or
#'   an integer vector of column indices.
#' @param focal The focal unit's own members, a subset of \code{set}. Never
#'   perturbed. These carry the evidence the finding is about.
#' @param pool For \code{form = "addition"}, the reservoir of probes to add,
#'   disjoint from \code{set}. Ignored for dropout.
#' @param form \code{"dropout"} (default) or \code{"addition"}.
#' @param grid Perturbation sizes. For dropout, fractions of the non-focal
#'   probes to remove; defaults to \code{c(.01, .02, .03, .05, .10, .20, .30,
#'   .50)}, which spans the observed range of cross-cohort probe-set divergence
#'   and then some. For addition, counts of probes to add; defaults to
#'   \code{c(1, 2, 5, 10, 25, 50)}.
#' @param R Replicates per grid point. Default 100. The reported quantity is a
#'   spread, so this is the resolution of the answer.
#' @param alpha Significance threshold used only to report whether a replicate's
#'   verdict differs from the unperturbed one. Default 0.05.
#' @param seed Optional integer for reproducibility. The caller's RNG state is
#'   restored on exit.
#' @return A \code{data.frame} with one row per engine, grid point and
#'   replicate: \code{engine}, \code{form}, \code{size} (the fraction dropped or
#'   the count added), \code{rep}, \code{p}, and \code{flip} - whether that
#'   replicate's verdict at \code{alpha} differs from the unperturbed verdict.
#'   The unperturbed run is included as \code{size = 0}, \code{rep = 0}. The
#'   baseline p-value per engine is attached as the \code{"baseline"} attribute.
#' @examples
#' set.seed(1)
#' n <- 120L; K <- 40L
#' f <- rnorm(n)
#' M <- sapply(seq_len(K), function(j) if (j <= 6) 0.7 * f + rnorm(n) else rnorm(n))
#' colnames(M) <- sprintf("cg%04d", seq_len(K))
#' y <- 0.6 * f + rnorm(n)
#' focal <- colnames(M)[1:6]
#'
#' # Two engines on identical evidence. The competitive one scores the focal
#' # unit against whatever else happens to be loaded; the unit-only one uses
#' # the focal probes and nothing else.
#' competitive <- function(u) {
#'   rest <- setdiff(u, focal)
#'   if (length(rest) < 2) return(NA_real_)
#'   t.test(abs(cor(M[, focal], y)), abs(cor(M[, rest], y)))$p.value
#' }
#' unit_only <- function(u) {
#'   summary(lm(rowMeans(M[, focal]) ~ y))$coefficients["y", 4]
#' }
#'
#' pd <- dmsa_pdance(list(competitive = competitive, unit_only = unit_only),
#'                   set = colnames(M), focal = focal,
#'                   grid = c(0.1, 0.3, 0.5), R = 20, seed = 1)
#'
#' # Both start from the same evidence, and it is never touched.
#' round(attr(pd, "baseline"), 4)
#'
#' # How far each engine's p travels when only the OTHER probes change.
#' round(tapply(pd$p, pd$engine, function(x) diff(range(x, na.rm = TRUE))), 4)
#'
#' # And how often the verdict at .05 reverses. The unit-only engine cannot
#' # move, because nothing it reads has changed.
#' tapply(pd$flip, pd$engine, sum)
#' @export
dmsa_pdance <- function(engines, set, focal, pool = NULL,
                        form = c("dropout", "addition"),
                        grid = NULL, R = 100L, alpha = 0.05, seed = NULL) {
  form <- match.arg(form)

  if (!is.list(engines) || !length(engines))
    stop("`engines` must be a non-empty list of functions", call. = FALSE)
  if (is.null(names(engines)) || any(!nzchar(names(engines))))
    stop("every element of `engines` must be named; the names label the ",
         "output", call. = FALSE)
  if (!all(vapply(engines, is.function, logical(1))))
    stop("every element of `engines` must be a function of one argument ",
         "returning a p-value", call. = FALSE)

  if (!length(set)) stop("`set` is empty", call. = FALSE)
  if (!length(focal)) stop("`focal` is empty", call. = FALSE)
  miss <- setdiff(focal, set)
  if (length(miss))
    stop("`focal` must be a subset of `set`; ", length(miss),
         " member(s) are not in `set`, e.g. ",
         paste(utils::head(miss, 3), collapse = ", "), call. = FALSE)

  free <- setdiff(set, focal)
  if (!length(free))
    stop("`set` contains nothing outside `focal`, so there is nothing to ",
         "perturb. The p-dance test asks what happens when the OTHER probes ",
         "change.", call. = FALSE)

  if (identical(form, "addition")) {
    if (is.null(pool) || !length(pool))
      stop("`pool` is required for the addition form: it is the reservoir of ",
           "probes outside `set` that get added", call. = FALSE)
    overlap <- intersect(pool, set)
    if (length(overlap))
      stop("`pool` must be disjoint from `set`; ", length(overlap),
           " member(s) appear in both", call. = FALSE)
  }

  if (is.null(grid))
    grid <- if (identical(form, "dropout"))
      c(.01, .02, .03, .05, .10, .20, .30, .50) else c(1, 2, 5, 10, 25, 50)
  grid <- sort(unique(grid[grid > 0]))
  if (!length(grid)) stop("`grid` has no positive values", call. = FALSE)
  if (identical(form, "dropout") && any(grid >= 1))
    stop("for the dropout form `grid` holds fractions in (0, 1)", call. = FALSE)
  if (identical(form, "addition") && any(grid > length(pool)))
    stop("`grid` asks for more added probes than `pool` holds (",
         length(pool), ")", call. = FALSE)

  R <- as.integer(R)
  if (is.na(R) || R < 1L) stop("`R` must be a positive integer", call. = FALSE)

  if (!is.null(seed)) {
    ## restore the caller's RNG state on exit: a permutation seed is for
    ## reproducing THIS result, not for silently reseeding the user's session.
    ## withr sets the seed and, on exit, restores the caller's state - or
    ## removes it when there was none (E10); the draws are those of set.seed().
    withr::local_seed(seed)
  }

  run <- function(u) vapply(engines, function(e) {
    v <- tryCatch(as.numeric(e(u))[1], error = function(err) NA_real_)
    if (length(v)) v else NA_real_
  }, numeric(1))

  base_p <- run(set)
  base_sig <- !is.na(base_p) & base_p < alpha

  out <- list(data.frame(engine = names(engines), form = form, size = 0,
                         rep = 0L, p = unname(base_p), flip = FALSE,
                         stringsAsFactors = FALSE))

  for (g in grid) for (r in seq_len(R)) {
    u <- if (identical(form, "dropout")) {
      keep <- round((1 - g) * length(free))
      ## sample the ELEMENTS, never sample(n, k): with a single remaining
      ## integer column index k, sample(k, 1) draws from 1:k instead
      c(focal, if (keep > 0) free[sample.int(length(free), keep)] else free[0])
    } else {
      c(set, pool[sample.int(length(pool), g)])
    }
    p <- run(u)
    sig <- !is.na(p) & p < alpha
    out[[length(out) + 1L]] <-
      data.frame(engine = names(engines), form = form, size = g, rep = r,
                 p = unname(p), flip = sig != base_sig,
                 stringsAsFactors = FALSE)
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "baseline") <- base_p
  attr(res, "alpha") <- alpha
  res
}
