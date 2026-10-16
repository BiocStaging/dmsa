# ============================================================================
# DMSA 2.0: ONE SCORE, ONE MODEL - at every level
#
# WHAT WAS ALREADY HERE. mDMSA (R/moderate.R) already collapsed an aligned set
# into a subject-level tone score with dmsa_score(), and already tested
# score x moderator with Freedman-Lane block permutation via dmsa_moderate().
# The idea that the aligned pooling can be applied to the DATA rather than to
# the coefficients is therefore not new to this file.
#
# WHAT THIS FILE ADDS.
#   dmsa_scores()          the same aligned projection PLUS the three
#                          alternatives every competing method implies -
#                          unweighted mean z, the regional first principal
#                          component, and sign-blind mean |z| - so the value of
#                          alignment is measurable instead of assumed. Also
#                          winsorisation, which won the positive-control
#                          estimator grid.
#   dmsa_model()           one model for ANY formula and ANY term, not just
#                          score x moderator: main effects, interactions,
#                          three-way terms, mutual adjustment, whatever the
#                          question is. This is what makes the upper level a
#                          single ordinary regression.
#   dmsa_model_omnibus()   ACAT across the projections. Valid under arbitrary
#                          dependence, so hedging costs no multiplicity.
#
# THE LEVEL STRUCTURE IS UNTOUCHED, and that is the point. dmsa_scores() takes
# whatever probe set you hand it, so a system score is level 1, a gene score is
# level 2, a probe is level 3. At every level the upper question is one
# regression on one number per person. A flat system score dilutes a one-gene
# signal by construction - AVP is 8 of oxytocin's 52 direction-called CpGs - so
# routing the question to the right level is not a workaround, it IS the method.
# ============================================================================

#' Aligned set scores, and the alternatives worth comparing against
#'
#' The aligned projection is \code{dmsa_score()}'s, restated for a matrix input
#' and with optional winsorisation. The other three are what competing methods
#' implicitly use, returned so that the comparison can be run rather than
#' argued.
#'
#' @param M Numeric matrix, samples x probes, on the analysis scale (M-values).
#'   Columns are standardised internally.
#' @param alignment A \code{dmsa_align()} result for the same probes in the same
#'   order, or a numeric vector of multipliers.
#' @param method \code{"expected"} uses \code{2*p_s_plus - 1}; \code{"fixed"}
#'   uses \code{s}.
#' @param winsor Winsorise each probe at this many MADs first; \code{NULL}
#'   disables. 3 is the default because it won the positive-control estimator
#'   grid - outliers otherwise dominate the composite.
#' @param flavours Which projections to return. \code{aligned} is DMSA's;
#' @param weighting Character. Probe weighting engine within a unit: \code{"combined"} (default) fuses the flat and reliability statistics on one shared permutation stream, \code{"flat"} weights every usable aligned probe equally, \code{"reliability"} weights each probe by its item-rest correlation with the rest of its unit. Weights are computed from methylation alone, so the permutation null is unaffected.
#' @param w_floor Numeric. Lower bound applied to reliability weights before normalisation, so a single poorly behaved probe cannot be driven to zero influence. Ignored when \code{weighting = "flat"}.
#'   \code{mean} is the unweighted mean z that "averaging the gene set" means;
#'   \code{pc1} is the regional principal component (DMRPC-style), oriented by
#'   the alignment so its sign stays interpretable; \code{blind} is mean
#'   \code{|z|}, the sign-blind control.
#' @return data.frame of standardised scores, one column per flavour.
#' @seealso \code{\link{dmsa_score}} for the single aligned vector on a
#'   data.frame, \code{\link{dmsa_model}} for the upper-level model.
#' @examples
#' \dontrun{
#' al <- dmsa_align(dir_table, genes = gene, level = "gene")
#' S  <- dmsa_scores(M[, probes], al)
#' dat$oxytocin <- S$aligned
#' }
#' @export
dmsa_scores <- function(M, alignment, method = c("expected", "fixed"),
                        winsor = 3,
                        flavours = c("aligned", "mean", "pc1", "blind"),
                        weighting = c("combined", "reliability", "flat"), w_floor = 1.5) {
  method <- match.arg(method); weighting <- match.arg(weighting)
  M <- as.matrix(M); mode(M) <- "numeric"
  m <- if (is.numeric(alignment)) alignment else {
    al <- as.data.frame(alignment)
    ## E9: multipliers are applied by position, so when both sides carry
    ## probe names the positional assumption is CHECKED (and an unambiguous
    ## scramble of the same probes is repaired) instead of trusted
    if (!is.null(al$probe) && nrow(al) == ncol(M)) {
      .ord <- .dmsa_align_order(colnames(M), al$probe)
      if (is.numeric(.ord)) al <- al[.ord, , drop = FALSE]
    }
    mm <- if (method == "fixed") as.numeric(al$s) else
      2 * as.numeric(al$p_s_plus) - 1
    us <- if (!is.null(al$usable)) as.logical(al$usable) else rep(TRUE, length(mm))
    mm[!is.finite(mm) | !us] <- 0
    mm
  }
  if (length(m) != ncol(M))
    stop("alignment has ", length(m), " entries but M has ", ncol(M),
         " columns; they must be in the same order", call. = FALSE)
  if (!is.null(winsor)) M <- apply(M, 2, function(y) {
    md <- stats::median(y, na.rm = TRUE); s <- stats::mad(y, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(y)
    pmin(pmax(y, md - winsor * s), md + winsor * s) })
  Z <- scale(M)
  keep <- apply(Z, 2, function(z) all(is.finite(z)))
  Z <- Z[, keep, drop = FALSE]; m <- m[keep]
  if (!ncol(Z)) stop("no usable probe columns", call. = FALSE)
  if (all(m == 0) && "aligned" %in% flavours)
    stop("every multiplier is zero: this set carries no usable direction call, ",
         "so an ALIGNED score is undefined. Ask for ",
         "flavours = c('mean','pc1','blind'), or supply directions.",
         call. = FALSE)
  ## reliability weighting (outcome-free): scale the aligned multiplier by each
  ## probe's item-rest centrality; the whole probe set is one unit here. A
  ## single score vector cannot be "both engines", so "combined" uses the
  ## reliability score (the sharper point estimate); combined-engine hedging
  ## lives at the inference level (dmsa_triangulate / dmsa_change).
  wt_eff <- if (weighting == "combined") "reliability" else weighting
  mw <- m * dmsa_relweights(Z, rep(1L, ncol(Z)), m, w_floor = w_floor,
                            weighting = wt_eff)
  al <- as.numeric(Z %*% mw) / sum(abs(mw))
  out <- list()
  if ("aligned" %in% flavours) out$aligned <- al
  if ("mean"    %in% flavours) out$mean    <- rowMeans(Z)
  if ("blind"   %in% flavours) out$blind   <- rowMeans(abs(Z))
  if ("pc1"     %in% flavours) out$pc1 <- tryCatch({
    p <- stats::prcomp(Z, center = TRUE, scale. = FALSE); s <- p$x[, 1]
    ## E10 (PI-approved): with all multipliers zero, `al` is 0/0 = NaN and
    ## cor(s, al) threw inside this tryCatch, silently turning pc1 into an
    ## all-NA column - in exactly the situation the aligned-score error
    ## message steers users to pc1. Orient against `al` only when an aligned
    ## reference exists; otherwise return the (arbitrarily signed) PC1 as is.
    if (all(is.finite(al))) {
      cr <- suppressWarnings(stats::cor(s, al))
      if (is.finite(cr) && cr < 0) s <- -s
    }
    s }, error = function(e) rep(NA_real_, nrow(Z)))
  as.data.frame(lapply(out[flavours[flavours %in% names(out)]],
                       function(v) as.numeric(scale(v))))
}

#' One model, block-permutation calibrated, for any term
#'
#' The upper level of DMSA 2.0. Fits \code{formula} once, then permutes the
#' reduced-model residuals within permutation blocks of equal size
#' (Freedman-Lane) to obtain a permutation-calibrated p-value for one named
#' coefficient. It is finite-sample exact only when the reduced-model residuals
#' are exchangeable under the tested null; otherwise it is an approximate
#' residual-permutation test, including for an interaction. \code{\link{dmsa_moderate}} is the special case of
#' this for \code{score * moderator} with a declared \code{dmsa_design()};
#' use that when you want the design contract enforced and simple slopes
#' reported, and this when the question is an arbitrary term.
#'
#' @param formula Model formula.
#' @param data data.frame.
#' @param term Name of the coefficient to test, exactly as \code{model.matrix}
#'   names it (e.g. \code{"S:ACE"}).
#' @param block Block labels, one per row of \code{data} - the family or
#'   cluster identifier. Rows are permuted only within blocks of EQUAL SIZE,
#'   because concatenating a random order of unequal-size families is not a
#'   valid relabelling.
#' @param B Permutations.
#' @param seed Optional integer.
#' @param nulls Logical (default \code{FALSE}). \code{TRUE} returns the
#'   vector of permutation null statistics as \code{$null_t} alongside the
#'   test, so a caller can compose joint corrections on one stream. (Joint
#'   calibration across lenses, engines and levels comes from sharing the
#'   \code{seed}, which \code{dmsa_model_omnibus()} does for you.)
#' @return An object of class \code{dmsa_model}.
#' @examples
#' set.seed(1)
#' n <- 80
#' d <- data.frame(S = rnorm(n), ACE = rnorm(n),
#'                 cID = rep(seq_len(n / 2), each = 2))   # couples
#' d$y <- 0.5 * d$S * d$ACE + rnorm(n)
#' ## rows are relabelled only within equal-size families, so the couple
#' ## dependence is preserved under the null
#' dmsa_model(y ~ S * ACE, d, term = "S:ACE", block = d$cID, B = 99, seed = 1)
#' @export
dmsa_model <- function(formula, data, term, block = NULL, B = 1999,
                       seed = NULL, nulls = FALSE) {
  if (!is.null(seed)) {
    ## restore the caller's RNG state on exit: a permutation seed is for
    ## reproducing THIS result, not for silently reseeding the user's session.
    .old_seed <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(.old_seed))
      assign(".Random.seed", .old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv()))
        ## E10: a true restore of "no prior RNG state" - otherwise a fresh
        ## session leaves deterministically seeded mid-stream of `seed`
        rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  block <- .dmsa_rows(block, as.data.frame(data), "block")
  .dmsa_check_block(block)
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  keep <- match(rownames(mf), rownames(as.data.frame(data)))
  X <- stats::model.matrix(formula, mf); y <- stats::model.response(mf)
  fi <- which(colnames(X) == term)
  if (!length(fi))
    stop("term '", term, "' is not a column of the design. Columns are: ",
         paste(utils::head(colnames(X), 20), collapse = ", "), call. = FALSE)
  n <- nrow(X); XtXi <- solve(crossprod(X))
  bh <- XtXi %*% crossprod(X, y); r <- y - X %*% bh
  dfr <- n - ncol(X); se <- sqrt(sum(r^2) / dfr * XtXi[fi, fi])
  tobs <- as.numeric(bh[fi] / se)
  Zo <- X[, -fi, drop = FALSE]; PZ <- Zo %*% solve(crossprod(Zo), t(Zo))
  fit <- PZ %*% y; res <- y - fit
  blk <- if (is.null(block)) seq_len(n) else block[keep]
  rows <- split(seq_len(n), blk)
  strata <- split(seq_along(rows), lengths(rows))
  tn <- numeric(B)
  for (b in seq_len(B)) {
    ix <- integer(n)
    for (st in strata) ix[unlist(rows[st], use.names = FALSE)] <-
      unlist(rows[st[sample.int(length(st))]], use.names = FALSE)
    yp <- fit + res[ix]
    bp <- XtXi %*% crossprod(X, yp); rp <- yp - X %*% bp
    tn[b] <- bp[fi] / sqrt(sum(rp^2) / dfr * XtXi[fi, fi])
  }
  structure(list(term = term, n = n, B = B, n_blocks = length(rows),
    b = as.numeric(bh[fi]), se = se, t = tobs,
    r_partial = tobs / sqrt(tobs^2 + dfr),
    p_perm = (1 + sum(abs(tn) >= abs(tobs))) / (B + 1),
    p_param = 2 * stats::pt(-abs(tobs), dfr),
    null_t = if (nulls) tn else NULL), class = "dmsa_model")
}

#' @export
print.dmsa_model <- function(x, ...) {
  cat("DMSA one-model test of '", x$term, "'\n", sep = "")
  cat(sprintf("  n = %d in %d blocks\n", x$n, x$n_blocks))
  cat(sprintf("  b = %+.4f (SE %.4f)   t = %+.2f   partial r = %+.3f\n",
              x$b, x$se, x$t, x$r_partial))
  cat(sprintf("  p = %.4f  (block permutation, B = %d)   parametric p = %.4g\n",
              x$p_perm, x$B, x$p_param))
  invisible(x)
}

#' Omnibus across the score projections
#'
#' Refits the same model once per projection, substituting each into the column
#' named \code{score_col}, and combines the p-values with the aggregated Cauchy
#' test. ACAT is valid under arbitrary dependence, so hedging across projections
#' costs no multiplicity - unlike a Bonferroni over arms, which is what the
#' old \code{min(1, 2*min(p))} gate amounted to.
#'
#' @param scores data.frame from \code{dmsa_scores()}.
#' @param formula Formula referring to the score by \code{score_col}.
#' @param data,block,B,seed As in \code{dmsa_model()}.
#' @param term Term to test.
#' @param score_col Column name the formula uses for the score.
#' @return list: \code{p_by_flavour}, the omnibus \code{p}, and the fits.
#' @examples
#' set.seed(1)
#' n <- 80
#' M <- matrix(rnorm(n * 6), n, 6)
#' al <- dmsa_align(data.frame(cpg = paste0("cg", 1:6), d = c(1, 1, -1, 1, -1, 1),
#'                             p_plus = rep(.9, 6)),
#'                  genes = rep("OXTR", 6), level = "gene")
#' sc <- dmsa_scores(M, al)
#' d <- data.frame(age = rnorm(n), cID = rep(seq_len(n / 2), each = 2))
#' d$y <- 0.6 * sc$aligned + rnorm(n)
#' ## one p per projection, combined by ACAT - hedging costs no multiplicity
#' o <- dmsa_model_omnibus(sc, y ~ S + age, d, term = "S", block = d$cID,
#'                         B = 99, seed = 1)
#' round(c(o$p_by_flavour, omnibus = o$p), 3)
#' @export
dmsa_model_omnibus <- function(scores, formula, data, term, block = NULL,
                               score_col = "S", B = 1999, seed = NULL) {
  ps <- stats::setNames(rep(NA_real_, ncol(scores)), names(scores))
  fits <- list()
  for (fl in names(scores)) {
    d2 <- as.data.frame(data); d2[[score_col]] <- scores[[fl]]
    f <- tryCatch(dmsa_model(formula, d2, term, block = block, B = B,
                             seed = seed), error = function(e) NULL)
    if (!is.null(f)) { ps[fl] <- f$p_perm; fits[[fl]] <- f }
  }
  p <- ps[is.finite(ps)]
  omni <- if (!length(p)) NA_real_ else {
    q <- pmin(pmax(p, 1e-15), 1 - 1e-15)
    0.5 - atan(mean(tan((0.5 - q) * pi))) / pi
  }
  list(p_by_flavour = ps, p = omni, fits = fits)
}
