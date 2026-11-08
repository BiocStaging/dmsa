# mDMSA: moderated Directional Methylation Set Analysis
# --------------------------------------------------------------------------
# Methylation rarely acts on a psychological outcome unconditionally. mDMSA
# collapses the aligned set into one subject-level tone score and tests its
# interaction with a moderator, under the same declared design and the same
# block-permutation discipline as dmsa_fit().

#' Subject-level aligned tone score
#'
#' Weighted mean of within-probe standardised M-values, each weighted by its
#' aligned multiplier (\code{s_j} for fixed-sign, \code{E[s_j]} for
#' expected-sign), then standardised. Positive = the set reads as higher
#' system activation tone.
#'
#' @param data data.frame.
#' @param probes methylation columns, in the row order of \code{alignment}.
#' @param alignment a \code{dmsa_align()} result.
#' @param method \code{"expected"} or \code{"fixed"}.
#' @param beta_input TRUE if \code{probes} hold beta values.
#' @return numeric vector, one value per row of \code{data} (NA where the
#'   probes are missing).
#' @examples
#' set.seed(1)
#' n <- 80
#' d <- rep(c(1, -1), each = 3)      # methylation-to-expression direction
#' f <- rnorm(n)                     # the tone the gene actually tracks
#' probes <- paste0("cg", 1:6)
#' dat <- as.data.frame(sapply(d, function(dj) plogis(dj * .8 * f + rnorm(n))))
#' names(dat) <- probes
#'
#' al <- dmsa_align(data.frame(cpg = probes, d = d, p_plus = ifelse(d > 0, .9, .1)),
#'                  genes = rep("FKBP5", 6), level = "gene")
#' ## one standardised tone per subject, despite the probes opposing each other
#' dat$tone <- dmsa_score(dat, probes, al)
#' round(cor(dat$tone, f), 2)
#' @export
dmsa_score <- function(data, probes, alignment, method = c("expected", "fixed"),
                       beta_input = TRUE) {
  method <- match.arg(method)
  al <- as.data.frame(alignment)
  if (length(probes) != nrow(al))
    stop("probes and alignment must be in the same order", call. = FALSE)
  mult <- if (method == "fixed") al$s else (2 * al$p_s_plus - 1)
  use <- !is.na(mult) & mult != 0
  if (!any(use)) stop("no probe carries a usable aligned sign", call. = FALSE)
  M <- as.matrix(data[, probes[use], drop = FALSE]); mode(M) <- "numeric"
  if (beta_input) M <- dmsa_mvalues(M)
  as.numeric(scale(scale(M) %*% mult[use] / sum(abs(mult[use]))))
}

#' Moderated DMSA (mDMSA)
#'
#' @param data data.frame.
#' @param probes,alignment as in \code{dmsa_fit()}.
#' @param outcome Character. The outcome column (a psychological measure).
#' @param moderator Character. The moderator column.
#' @param design A \code{dmsa_design()}. Its \code{focal} is ignored and
#'   replaced by \code{score * moderator}; its \code{fixed}, \code{random},
#'   \code{exchangeable} and \code{forbid} are used as declared.
#' @param method \code{"expected"} or \code{"fixed"} scoring.
#' @param center Center the moderator (default TRUE). Leaving a moderator
#'   uncentered makes the product term collinear with its parents and the
#'   simple slopes uninterpretable; the returned object records what was done.
#' @param scale_outcome Standardise the outcome (default TRUE).
#' @param B Permutations (Freedman-Lane: residuals of the no-interaction model
#'   are permuted within exchangeable blocks, preserving the dependence).
#' @param engine as in \code{dmsa_fit()}.
#' @param beta_input TRUE if \code{probes} hold beta values.
#' @param check Run the design check (default TRUE).
#' @param seed Optional integer.
#' @return An object of class \code{dmsa_moderate}.
#' @examples
#' set.seed(1)
#' n <- 100; probes <- paste0("cg", 1:6)
#' dat <- data.frame(sex_c = rep(c(-.5, .5), length.out = n),
#'                   chip = rep(1:10, each = 10),
#'                   cID = rep(seq_len(n / 2), each = 2), mod = rnorm(n))
#' dat[probes] <- plogis(matrix(rnorm(n * 6), n, 6))
#' al <- dmsa_align(data.frame(cpg = probes, d = rep(c(1, -1), 3),
#'                             p_plus = rep(c(.9, .1), 3)),
#'                  genes = rep("OXTR", 6), level = "gene")
#' dat$out <- 0.6 * dmsa_score(dat, probes, al) * dat$mod + rnorm(n)
#' des <- dmsa_design("placeholder", "sex_c", random = "chip", exchangeable = "cID")
#' dmsa_moderate(dat, probes, al, outcome = "out", moderator = "mod",
#'               design = des, B = 99, engine = "lm", seed = 3)
#' @export
dmsa_moderate <- function(data, probes, alignment, outcome, moderator, design,
                          method = c("expected", "fixed"),
                          center = TRUE, scale_outcome = TRUE,
                          B = 999L, engine = c("auto", "lm", "lmer"),
                          beta_input = TRUE, check = TRUE, seed = NULL) {
  method <- match.arg(method); engine <- match.arg(engine)
  stopifnot(inherits(design, "dmsa_design"))
  if (!is.null(seed)) {
    ## restore the caller's RNG state on exit: a permutation seed is for
    ## reproducing THIS result, not for silently reseeding the user's session.
    ## withr sets the seed and, on exit, restores the caller's state - or
    ## removes it when there was none (E10); the draws are those of set.seed().
    withr::local_seed(seed)
  }
  data <- as.data.frame(data)

  bad <- intersect(design$forbid, c(outcome, moderator))
  if (length(bad))
    stop("never-list term used as outcome/moderator: ", paste(bad, collapse = ", "),
         call. = FALSE)

  data$.S <- dmsa_score(data, probes, alignment, method = method,
                        beta_input = beta_input)
  mvars <- unique(c(".S", outcome, moderator, design$fixed,
                    design$random_groups, design$exchangeable))
  miss <- setdiff(mvars, names(data))
  if (length(miss)) stop("columns absent from data: ", paste(miss, collapse = ", "),
                         call. = FALSE)

  d2 <- design
  d2$focal <- c(".S:.M", ".S", ".M")
  d2$focal_test <- ".S:.M"; d2$focal_vars <- c(".S", ".M")
  d2$vars <- unique(c(".S", ".M", outcome, design$fixed, design$random_groups,
                      design$exchangeable))
  data$.M <- data[[moderator]]
  mod_mean <- mean(data$.M, na.rm = TRUE)
  if (center) data$.M <- data$.M - mod_mean

  keep <- stats::complete.cases(data[, c(d2$vars, outcome), drop = FALSE])
  dat  <- data[keep, , drop = FALSE]
  n    <- nrow(dat)
  dat$.y <- if (scale_outcome) as.numeric(scale(dat[[outcome]])) else dat[[outcome]]
  if (check) dmsa_check_design(d2, dat)

  sdM <- stats::sd(dat$.M)
  rhs_full <- paste(c(".S*.M", design$fixed), collapse = " + ")
  rhs_red  <- paste(c(".S", ".M", design$fixed), collapse = " + ")
  f_full <- stats::as.formula(paste(".y ~", rhs_full))
  f_red  <- stats::as.formula(paste(".y ~", rhs_red))

  fit_lm <- stats::lm(f_full, dat)
  cl <- stats::coef(summary(fit_lm))
  ix <- if (".S:.M" %in% rownames(cl)) ".S:.M" else ".M:.S"
  b_lm <- cl[ix, 1]; se_lm <- cl[ix, 2]
  V_lm <- stats::vcov(fit_lm)

  ## collinearity of the product term with its parents -- the centering check
  vif_int <- {
    r2 <- summary(stats::lm(stats::as.formula(
      paste(c("I(.S*.M) ~ .S + .M", design$fixed), collapse = " + ")), dat))$r.squared
    1 / max(1 - r2, 1e-12)
  }

  use_lmer <- engine == "lmer" ||
    (engine == "auto" && length(design$random_groups) > 0 &&
       requireNamespace("lme4", quietly = TRUE))
  b_rep <- b_lm; se_rep <- se_lm; p_rep <- cl[ix, 4]; icc <- NULL
  if (use_lmer) {
    fml <- stats::as.formula(paste(".y ~", rhs_full, "+",
             paste(sprintf("(1|%s)", design$random_groups), collapse = " + ")))
    fm <- tryCatch(suppressMessages(suppressWarnings(lme4::lmer(fml, dat))),
                   error = function(e) NULL)
    if (!is.null(fm)) {
      cf <- stats::coef(summary(fm))
      k <- if (ix %in% rownames(cf)) ix else grep(":", rownames(cf), value = TRUE)[1]
      b_rep <- cf[k, 1]; se_rep <- cf[k, 2]
      p_rep <- if (ncol(cf) >= 5) cf[k, 5] else 2 * stats::pnorm(-abs(cf[k, 1] / cf[k, 2]))
      vc <- as.data.frame(lme4::VarCorr(fm))
      icc <- stats::setNames(vc$vcov / sum(vc$vcov), vc$grp)
    }
  }

  ## ---- Freedman-Lane: permute reduced-model residuals within blocks --------
  p_perm <- NA_real_; null_t <- numeric(0)
  if (B > 0) {
    r <- stats::lm(f_red, dat); fv <- stats::fitted(r); ev <- stats::resid(r)
    blk <- if (!is.null(design$exchangeable)) dat[[design$exchangeable]] else seq_len(n)
    idxs <- dmsa_block_index(blk, B)
    t_obs <- b_lm / se_lm
    null_t <- vapply(idxs, function(i2) {
      dd <- dat; dd$.y <- fv + ev[i2]
      cc <- stats::coef(summary(stats::lm(f_full, dd)))
      k <- if (ix %in% rownames(cc)) ix else grep(":", rownames(cc), value = TRUE)[1]
      cc[k, 1] / cc[k, 2]
    }, numeric(1))
    p_perm <- dmsa_perm_pvalue(t_obs, null_t)
  }

  ## ---- simple slopes at +/- 1 SD of the moderator --------------------------
  ## E7 fix (2026-08-29, PI-approved): every ingredient - the slope, the
  ## interaction coefficient, and all three (co)variances - comes from ONE
  ## fitted model: the mixed model when it is the reported engine, the lm
  ## otherwise. The old code spliced b_S and the vcov from lm with b_rep
  ## from lmer, so the printed slope belonged to neither model.
  if (use_lmer && exists("fm") && !is.null(fm) &&
      ".S" %in% names(lme4::fixef(fm))) {
    Vr <- as.matrix(stats::vcov(fm))
    kI <- if (ix %in% rownames(Vr)) ix else
      grep(":", rownames(Vr), value = TRUE)[1]
    b_S  <- unname(lme4::fixef(fm)[".S"])
    b_I  <- unname(lme4::fixef(fm)[kI])
    v_SS <- Vr[".S", ".S"]; v_II <- Vr[kI, kI]; v_SI <- Vr[".S", kI]
    ## Wald/normal reference for the mixed model (its residual df is not
    ## uniquely defined); the lm branch keeps its exact t reference
    tc <- stats::qnorm(.975); .pfun <- function(t) 2 * stats::pnorm(-abs(t))
  } else {
    b_S <- cl[".S", 1]; b_I <- b_lm
    v_SS <- V_lm[".S", ".S"]; v_II <- V_lm[ix, ix]; v_SI <- V_lm[".S", ix]
    tc <- stats::qt(.975, stats::df.residual(fit_lm))
    .pfun <- function(t) 2 * stats::pt(-abs(t), df = stats::df.residual(fit_lm))
  }
  at <- c(-1, 1) * sdM
  ss <- data.frame(
    at_moderator = c("-1 SD", "+1 SD"),
    moderator_value = at + if (center) mod_mean else 0,
    slope = b_S + b_I * at,
    se = sqrt(v_SS + at^2 * v_II + 2 * at * v_SI))
  ss$t <- ss$slope / ss$se
  ss$p <- .pfun(ss$t)

  ## ---- Johnson-Neyman region ----------------------------------------------
  ## the slope is significant where A w^2 + Bq w + Cq <= 0. The parabola
  ## opens DOWN when A < 0 (significant OUTSIDE the roots) and UP when A > 0
  ## (significant BETWEEN the roots - the everyday strong-main-effect /
  ## weak-interaction case, which the old code reported inverted). With no
  ## real roots the sign of A alone decides: everywhere (A < 0) or nowhere.
  A <- v_II * tc^2 - b_I^2
  Bq <- 2 * (v_SI * tc^2 - b_S * b_I)
  Cq <- v_SS * tc^2 - b_S^2
  disc <- Bq^2 - 4 * A * Cq
  obs_range <- range(dat[[moderator]], na.rm = TRUE)
  if (A == 0 || disc < 0) {
    jn <- c(NA_real_, NA_real_)
    jn_region <- if (A < 0 || (A == 0 && Cq <= 0)) "everywhere" else "nowhere"
    in_region <- if (jn_region == "everywhere") 1 else 0
  } else {
    jn <- sort((-Bq + c(-1, 1) * sqrt(disc)) / (2 * A))
    jn_region <- if (A > 0) "between" else "outside"
    jn_raw0 <- jn + if (center) mod_mean else 0
    in_region <- if (jn_region == "between")
      mean(dat[[moderator]] >= jn_raw0[1] & dat[[moderator]] <= jn_raw0[2],
           na.rm = TRUE)
    else
      mean(dat[[moderator]] < jn_raw0[1] | dat[[moderator]] > jn_raw0[2],
           na.rm = TRUE)
  }
  jn_raw <- jn + if (center) mod_mean else 0

  structure(list(
    call = match.call(), design = design, outcome = outcome,
    moderator = moderator, method = method, n = n,
    centered = center, moderator_mean = mod_mean, moderator_sd = sdM,
    b_int = b_rep, se_int = se_rep, t_int = b_rep / se_rep, p_int = p_rep,
    p_perm = p_perm, vif_interaction = vif_int,
    engine = if (use_lmer) "lmer (reported) / lm (null)" else "lm",
    icc = icc, simple_slopes = ss,
    jn_bounds = jn_raw, jn_region = jn_region,
    jn_share = in_region, moderator_range = obs_range,
    null_t = null_t, B = B
  ), class = "dmsa_moderate")
}

#' @export
print.dmsa_moderate <- function(x, ...) {
  cat("mDMSA: ", x$outcome, " ~ tone x ", x$moderator, "\n", sep = "")
  if (!is.null(x$design$label)) cat("  ", x$design$label, "\n", sep = "")
  cat(sprintf("  n = %d   moderator %s (mean %.3f, SD %.3f)\n", x$n,
              if (x$centered) "centered" else "UNCENTERED - simple slopes are at 0",
              x$moderator_mean, x$moderator_sd))
  cat(sprintf("  interaction b = %.4f  se = %.4f  t = %.2f  p = %.4g\n",
              x$b_int, x$se_int, x$t_int, x$p_int))
  if (is.finite(x$p_perm))
    cat(sprintf("  Freedman-Lane permutation p = %.4g (B = %d, block = %s)\n",
                x$p_perm, x$B,
                if (is.null(x$design$exchangeable)) "rows" else x$design$exchangeable))
  cat(sprintf("  product-term VIF = %.2f%s\n", x$vif_interaction,
              if (x$vif_interaction > 5) "  <- collinear; check centering" else ""))
  cat("  engine ", x$engine, "\n", sep = "")
  cat("  simple slopes:\n")
  print(format(x$simple_slopes, digits = 3), row.names = FALSE)
  if (all(is.finite(x$jn_bounds))) {
    cat(sprintf("  Johnson-Neyman: significant %s [%.2f, %.2f]; %.0f%% of the sample",
                x$jn_region %||% "outside",
                x$jn_bounds[1], x$jn_bounds[2], 100 * x$jn_share))
    cat(sprintf("\n    (moderator observed range %.2f to %.2f)\n", x$moderator_range[1],
                x$moderator_range[2]))
  } else if (identical(x$jn_region, "everywhere"))
    cat("  Johnson-Neyman: significant across the whole moderator range\n")
  else cat("  Johnson-Neyman: no region of significance\n")
  invisible(x)
}
