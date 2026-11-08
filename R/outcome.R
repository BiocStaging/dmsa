# Non-Gaussian outcomes
# --------------------------------------------------------------------------
# Two regression directions live in this package and they behave differently
# when the psychological variable is not continuous.
#
# SET LEVEL (dmsa_fit, dmsa_cascade): methylation is the OUTCOME,
#   M_j ~ exposure + covariates, and the exposure's type is irrelevant - a
#   binary or dummy-coded predictor is perfectly ordinary in a linear model.
#   For a single focal term the t statistic for M ~ X equals the t statistic for
#   X ~ M (same partial correlation), so nothing about the set-level machinery
#   needs to change. Nothing in this file applies there.
#
# SUBJECT LEVEL (dmsa_moderate, and this file): the psychological variable is
#   the OUTCOME and the aligned tone score is the PREDICTOR. Here the outcome's
#   type matters, and `lm` on a 0/1 outcome is wrong.
#
# PERMUTATION. Freedman-Lane permutes residuals of the reduced OUTCOME model,
# which has no clean analogue outside the Gaussian case. So for non-Gaussian
# families the null is built on the PREDICTOR side instead: residualize the tone
# score (and, for a moderated test, the product term) on the covariates, permute
# those residuals within exchangeable blocks, and refit. The outcome's
# distribution is never touched, so one scheme serves binomial, ordinal,
# multinomial and count outcomes alike. Calibration of this scheme is checked by
# simulation rather than asserted.

## residualize a matrix of predictors on a covariate design
.residualize <- function(P, Z) {
  if (is.null(Z) || !ncol(Z)) return(scale(P, TRUE, FALSE))
  Q <- cbind(1, Z)
  P - Q %*% (solve(crossprod(Q), t(Q) %*% P))
}

#' Test an aligned tone score against a non-Gaussian outcome
#'
#' The subject-level arm of DMSA for outcomes that are not continuous: binary
#' (logistic), ordered-categorical (proportional odds), unordered-categorical
#' (multinomial), or counts (Poisson / negative binomial). Optionally moderated.
#'
#' @param data data.frame with one row per subject.
#' @param outcome Name of the outcome column. For \code{"binomial"} it may be
#'   0/1, logical, or a two-level factor; for \code{"ordinal"} an ordered factor
#'   or an integer-coded severity; for \code{"multinomial"} a factor.
#' @param score Name of the tone-score column (from \code{dmsa_score()}).
#' @param family One of \code{"gaussian"}, \code{"binomial"}, \code{"ordinal"},
#'   \code{"multinomial"}, \code{"poisson"}.
#' @param moderator Optional moderator column; when given, the tested term is
#'   \code{score:moderator}.
#' @param covariates Character vector of covariate columns.
#' @param block Optional column naming the exchangeable unit for permutation
#'   (families, repeated measures).
#' @param center Center the moderator (default TRUE).
#' @param B Permutations. 0 skips them and returns model-based inference only.
#' @param seed Optional integer.
#' @return An object of class \code{dmsa_outcome}.
#' @details
#' \strong{Direction.} A binary, count or ordinal outcome yields ONE coefficient
#' for the tone score, so DMSA's directional claim survives intact - the sign of
#' the log-odds (or log-rate, or latent-scale slope) is the expression-aligned
#' direction. A \strong{multinomial} outcome yields K-1 coefficients and there is
#' no single sign: the test becomes an omnibus "tone predicts category" and the
#' per-category log-odds are reported separately. If the categories are in any
#' sense ordered, \code{family = "ordinal"} keeps the direction and is strictly
#' preferable.
#' @examples
#' set.seed(1)
#' n <- 120
#' d <- data.frame(S = rnorm(n), cv = rnorm(n),
#'                 cid = rep(seq_len(n / 2), each = 2))
#' d$y <- rbinom(n, 1, plogis(0.9 * d$S))
#' ## a binary outcome yields ONE coefficient, so the sign of the log-odds
#' ## is still the expression-aligned direction
#' dmsa_outcome(d, outcome = "y", score = "S", family = "binomial",
#'              covariates = "cv", block = "cid", B = 99, seed = 2)
#' @export
dmsa_outcome <- function(data, outcome, score, family = c("gaussian", "binomial",
                         "ordinal", "multinomial", "poisson"),
                         moderator = NULL, covariates = character(),
                         block = NULL, center = TRUE, B = 999L, seed = NULL) {
  family <- match.arg(family)
  if (!is.null(seed)) {
    ## restore the caller's RNG state on exit: a permutation seed is for
    ## reproducing THIS result, not for silently reseeding the user's session.
    ## withr sets the seed and, on exit, restores the caller's state - or
    ## removes it when there was none (E10); the draws are those of set.seed().
    withr::local_seed(seed)
  }
  data <- as.data.frame(data)
  need <- c(outcome, score, moderator, covariates, block)
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("columns absent from data: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  dat <- data[stats::complete.cases(data[, need, drop = FALSE]), , drop = FALSE]
  n <- nrow(dat)
  if (n < 30) stop("only ", n, " complete cases", call. = FALSE)

  y <- dat[[outcome]]
  if (family == "binomial") {
    if (is.logical(y)) y <- as.integer(y)
    if (is.factor(y)) {
      if (nlevels(y) != 2L) stop("binomial outcome must have exactly 2 levels; ",
                                 "got ", nlevels(y), call. = FALSE)
      y <- as.integer(y) - 1L
    }
    u <- sort(unique(y))
    if (!identical(as.numeric(u), c(0, 1)))
      stop("binomial outcome must be 0/1, logical, or a 2-level factor; got ",
           paste(u, collapse = "/"), call. = FALSE)
    if (min(table(y)) < 10)
      warning("the rarer outcome class has ", min(table(y)),
              " cases; logistic estimates will be unstable")
  }
  if (family == "ordinal") {
    if (!is.factor(y)) y <- factor(y, ordered = TRUE)
    y <- factor(y, levels = levels(y), ordered = TRUE)
    if (nlevels(y) < 3L)
      stop("ordinal outcome needs >= 3 levels; use family = 'binomial'",
           call. = FALSE)
  }
  if (family == "multinomial") {
    y <- factor(y)
    if (nlevels(y) < 3L)
      stop("multinomial outcome needs >= 3 levels; use family = 'binomial'",
           call. = FALSE)
  }
  if (family == "poisson" && (any(y < 0) || any(y != round(y))))
    stop("poisson outcome must be non-negative integers", call. = FALSE)

  S <- as.numeric(dat[[score]])
  W <- if (is.null(moderator)) NULL else as.numeric(dat[[moderator]])
  mod_mean <- if (is.null(W)) NA_real_ else mean(W)
  if (!is.null(W) && center) W <- W - mod_mean
  Z <- if (length(covariates))
    stats::model.matrix(stats::reformulate(covariates), dat)[, -1, drop = FALSE]
  else NULL

  ## predictors: the tested one first
  if (is.null(W)) { P <- cbind(S = S); tested <- "S" }
  else { P <- cbind(SW = S * W, S = S, W = W); tested <- "SW" }

  fitfun <- function(Pm) {
    df <- data.frame(y = y, Pm, check.names = FALSE)
    if (!is.null(Z)) df <- cbind(df, as.data.frame(Z))
    rhs <- paste(setdiff(names(df), "y"), collapse = " + ")
    f <- stats::as.formula(paste("y ~", rhs))
    switch(family,
      gaussian    = stats::lm(f, df),
      binomial    = suppressWarnings(stats::glm(f, df, family = stats::binomial())),
      poisson     = suppressWarnings(stats::glm(f, df, family = stats::poisson())),
      ordinal     = suppressWarnings(MASS::polr(f, df, Hess = TRUE, method = "logistic")),
      multinomial = suppressWarnings(nnet::multinom(f, df, trace = FALSE)))
  }
  ## statistic for the tested term: |z| for one coefficient, LR chi-square for
  ## the multinomial block
  statfun <- function(fit, Pm) {
    if (family == "multinomial") {
      df <- data.frame(y = y, Pm, check.names = FALSE)
      if (!is.null(Z)) df <- cbind(df, as.data.frame(Z))
      drop <- setdiff(names(df), c("y", tested))
      f0 <- stats::as.formula(paste("y ~", paste(drop, collapse = " + ")))
      m0 <- suppressWarnings(nnet::multinom(f0, df, trace = FALSE))
      return(as.numeric(stats::deviance(m0) - stats::deviance(fit)))
    }
    cf <- if (inherits(fit, "polr")) {
      se <- sqrt(diag(stats::vcov(fit)))[tested]
      abs(stats::coef(fit)[tested] / se)
    } else {
      cc <- stats::coef(summary(fit))
      abs(cc[tested, 1] / cc[tested, 2])
    }
    as.numeric(cf)
  }

  fit <- fitfun(P)
  obs <- statfun(fit, P)

  ## reported estimate on the model's own scale
  est <- if (inherits(fit, "multinom")) NA_real_ else
    as.numeric(if (inherits(fit, "polr")) stats::coef(fit)[tested]
               else stats::coef(summary(fit))[tested, 1])
  se  <- if (inherits(fit, "multinom")) NA_real_ else
    as.numeric(if (inherits(fit, "polr")) sqrt(diag(stats::vcov(fit)))[tested]
               else stats::coef(summary(fit))[tested, 2])
  p_model <- if (inherits(fit, "multinom")) {
    stats::pchisq(obs, df = nlevels(y) - 1L, lower.tail = FALSE)
  } else 2 * stats::pnorm(-abs(est / se))

  ## per-category log-odds for the multinomial case, since there is no one sign
  per_category <- NULL
  if (inherits(fit, "multinom")) {
    cf <- stats::coef(fit)
    if (is.matrix(cf) && tested %in% colnames(cf))
      per_category <- data.frame(category = rownames(cf),
                                 log_odds = as.numeric(cf[, tested]),
                                 row.names = NULL)
  }

  ## ---- permutation on the predictor side ----------------------------------
  p_perm <- NA_real_; null_stat <- numeric(0)
  if (B > 0) {
    ## E10b (2026-08-29, PI-approved): FULL Freedman-Lane geometry. The
    ## tested term is residualised on every non-tested regressor of the
    ## fitted model - for the moderated test that means S and W as well as
    ## the covariates, not the covariates alone. Otherwise the observed
    ## statistic's regressor carries the product's collinearity with its
    ## parents while every permuted draw has that geometry broken, and the
    ## null is mis-centred exactly when S x W is collinear with S or W.
    .Zfull <- cbind(Z, P[, setdiff(colnames(P), tested), drop = FALSE])
    Pr <- .residualize(P[, tested, drop = FALSE], .Zfull)
    blk <- if (is.null(block)) NULL else dat[[block]]
    idxs <- if (is.null(blk)) lapply(seq_len(B), function(i) sample.int(n))
            else dmsa_block_index(blk, B)
    null_stat <- vapply(idxs, function(ix) {
      Pm <- P; Pm[, tested] <- Pr[ix, 1]
      statfun(fitfun(Pm), Pm)
    }, numeric(1))
    null_stat <- null_stat[is.finite(null_stat)]
    p_perm <- (1 + sum(null_stat >= obs)) / (length(null_stat) + 1)
  }

  structure(list(call = match.call(), family = family, outcome = outcome,
    score = score, moderator = moderator, tested = tested, n = n,
    n_levels = if (family %in% c("ordinal", "multinomial")) nlevels(y) else 2L,
    estimate = est, se = se, statistic = obs, p_model = p_model,
    p_perm = p_perm, B = length(null_stat), null_stat = null_stat,
    per_category = per_category, centered = center, moderator_mean = mod_mean,
    directional = family != "multinomial", fit = fit),
    class = "dmsa_outcome")
}

#' @export
print.dmsa_outcome <- function(x, ...) {
  cat("DMSA subject-level test (", x$family, ")\n", sep = "")
  cat("  ", x$outcome, " ~ ", x$score,
      if (!is.null(x$moderator)) paste0(" x ", x$moderator), "   n = ", x$n,
      if (x$family %in% c("ordinal", "multinomial"))
        paste0(", ", x$n_levels, " levels"), "\n", sep = "")
  if (x$family == "multinomial") {
    cat(sprintf("  omnibus LR chi-square = %.2f on %d df   p(model) = %.4g\n",
                x$statistic, x$n_levels - 1L, x$p_model))
    cat("  NO single direction: a multinomial outcome has one coefficient per\n",
        "  category. Per-category log-odds for the tested term:\n", sep = "")
    if (!is.null(x$per_category)) print(x$per_category, row.names = FALSE)
    cat("  If the categories are ordered at all, family = 'ordinal' keeps the\n",
        "  directional claim and is preferable.\n", sep = "")
  } else {
    sc <- switch(x$family, binomial = "log-odds", poisson = "log-rate",
                 ordinal = "latent-scale slope", "coefficient")
    cat(sprintf("  %s = %+.4f  se = %.4f  |z| = %.2f  p(model) = %.4g\n",
                sc, x$estimate, x$se, x$statistic, x$p_model))
    if (x$family %in% c("binomial", "poisson"))
      cat(sprintf("  exp(coef) = %.3f\n", exp(x$estimate)))
  }
  if (is.finite(x$p_perm))
    cat(sprintf("  permutation p = %.4g  (B = %d, predictor-side residual permutation)\n",
                x$p_perm, x$B))
  else cat("  permutation not run (B = 0): model-based inference only\n")
  invisible(x)
}
