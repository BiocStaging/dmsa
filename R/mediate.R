# ============================================================================
# LATENT DMSA: moderated mediation over latent system tones
#
# The claim shape this file serves: "W moderates the X -> M links and M is
# directly associated with Y", with X, M, W latent SYSTEM TONES measured by
# their MODULE tones as reflective indicators. Observed-composite mediation
# attenuates exactly the interaction paths that carry such theories (Feng et
# al. 2020, Front. Psychol.: ~8% underestimation at reliability ~.85, ~57%
# at ~.6), so the structural model here runs on measurement-error-corrected
# moments instead of on composites.
#
# ARCHITECTURE (registered in SPEC_latent_DMSA_v2):
#   measurement  one-factor model per system on its module tones.
#                Default is the essentially tau-equivalent model - equal
#                loadings, free uniquenesses (Strauss 2026, EPM: congeneric
#                one-factor ML commonly fails below n = 150 with 5
#                indicators; tau-equivalent converges at ~100% there, and
#                Heywood cases must NOT be patched to zero). Ladder:
#                congeneric (sensitivity) -> tau-equivalent (default) ->
#                reliability-corrected single indicator (fallback,
#                Savalei 2019).
#   structure    Bartlett scores + SAM/Croon corrected moments (Rosseel,
#                Burghgraeve, Loh & Schermelleh-Engel 2025, Behav. Res.
#                Methods; Croon 2002; Cox & Kelcey 2021): the moment matrix
#                of (latents, latent products, covariates, outcome) is the
#                factor-score moment matrix minus the measurement-error
#                offsets, with the Bogaert-Loh-Rosseel (2022) alpha/lambda
#                small-sample corrections keeping it positive definite.
#                Bartlett scores from DISJOINT indicator sets make the
#                cross-products unbiased for the latent products; only
#                variances and product-product covariances need offsets.
#   inference    Freedman-Lane permutation per structural path (the
#                measurement fit is untouched by permutation - residual
#                permutation never resamples people, so the indicator
#                covariance matrix is literally unchanged; p-values are
#                conditional on the measurement model). Products are tested
#                by JOINT SIGNIFICANCE, max(p_a, p_b) - the permutation
#                supremum test of Kroehl, Lutz & Wagner (2020), the valid
#                small-sample choice; permuting the product statistic
#                itself is grossly liberal at the composite null (Type I to
#                ~.78: Taylor & MacKinnon 2012). Holm over the declared
#                index family. Estimation uncertainty comes from a
#                nonparametric bootstrap of the ENTIRE pipeline -
#                measurement models refit in every resample (Can & Rosseel
#                2025: second-stage-only resampling understates
#                uncertainty).
#   chip         quasi-demeaning (R/ranef.R machinery) applied to module
#                tones BEFORE the measurement stage and to outcome/
#                covariates before the structural stage. Permutation is
#                available unrestricted or within-chip block-restricted
#                (Winkler et al. 2015); the registered simulation gate
#                selects the scheme.
# ============================================================================

## ---------------------------------------------------------------------------
## measurement stage: one-factor ML on a covariance matrix
## ---------------------------------------------------------------------------

## ML discrepancy for a 1-factor model. S: covariance matrix. Returns fitted
## lambda/theta. model = "tau" constrains one loading for all indicators.
.med_fa1 <- function(S, n, model = c("tau", "congeneric"),
                     theta_floor = 0.005) {
  model <- match.arg(model)
  p <- ncol(S); dS <- diag(S)
  floor_j <- theta_floor * dS
  ev <- eigen(stats::cov2cor(S), symmetric = TRUE)
  l0 <- sqrt(max(ev$values[1] - 1, 0.1) / p) * sqrt(dS)  # rough start
  fml <- function(lam, th) {
    Sig <- tcrossprod(lam) + diag(th, p)
    ch <- tryCatch(chol(Sig), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    ld <- 2 * sum(log(diag(ch)))
    tr <- sum(diag(chol2inv(ch) %*% S))
    v <- ld + tr - determinant(S, logarithm = TRUE)$modulus[1] - p
    if (!is.finite(v)) 1e10 else v
  }
  if (model == "tau") {
    par0 <- c(mean(l0), log(pmax(dS - mean(l0)^2, floor_j)))
    obj <- function(par) {
      lam <- rep(par[1], p); th <- pmax(exp(par[-1]), floor_j)
      fml(lam, th)
    }
    o <- stats::optim(par0, obj, method = "BFGS",
                      control = list(maxit = 500, reltol = 1e-10))
    lam <- rep(o$par[1], p); th <- pmax(exp(o$par[-1]), floor_j)
  } else {
    par0 <- c(l0, log(pmax(dS - l0^2, floor_j)))
    obj <- function(par) {
      lam <- par[seq_len(p)]; th <- pmax(exp(par[-seq_len(p)]), floor_j)
      fml(lam, th)
    }
    o <- stats::optim(par0, obj, method = "BFGS",
                      control = list(maxit = 1000, reltol = 1e-10))
    lam <- o$par[seq_len(p)]; th <- pmax(exp(o$par[-seq_len(p)]), floor_j)
  }
  if (mean(lam) < 0) lam <- -lam            # orient the factor with its tones
  heywood <- any(th <= floor_j * (1 + 1e-6))
  conv <- is.finite(o$value) && o$convergence == 0L
  list(lambda = lam, theta = th, conv = conv, heywood = heywood,
       model = model, value = o$value)
}

## Bartlett mapping for one fitted system. Z: n x p CENTERED indicator
## matrix. Returns the score vector (unbiased: f = eta + u), the error
## variance V = Var(u) = (lambda' Theta^-1 lambda)^-1, omega, determinacy.
.med_bartlett <- function(Z, lambda, theta) {
  a <- lambda / theta
  denom <- sum(lambda * a)                  # lambda' Theta^-1 lambda
  f <- as.numeric(Z %*% a) / denom
  Sig <- tcrossprod(lambda) + diag(theta, length(theta))
  det_r <- sqrt(max(0, min(1, as.numeric(
    crossprod(lambda, solve(Sig, lambda))))))
  omega <- sum(lambda)^2 / (sum(lambda)^2 + sum(theta))
  list(f = f, V = 1 / denom, omega = omega, determinacy = det_r)
}

## The measurement ladder for one system. tones: n x p matrix of module
## tones (already chip-transformed if requested). Returns scores + metadata.
.med_measure <- function(tones, measurement = c("tau", "congeneric", "single"),
                         label = "system", single_rel = NULL,
                         mod_w = NULL) {
  measurement <- match.arg(measurement)
  p <- ncol(tones); n <- nrow(tones)
  Z <- scale(tones, scale = FALSE)          # center only: keep variances
  used <- measurement; note <- character()
  if (p < 3L && measurement != "single") {
    note <- c(note, sprintf("%s: %d indicator(s) < 3; single-indicator fallback",
                            label, p))
    used <- "single"
  }
  if (used != "single") {
    S <- stats::cov(tones)
    fit <- .med_fa1(S, n, model = if (used == "congeneric") "congeneric" else "tau")
    if (used == "congeneric" && (!fit$conv || fit$heywood)) {
      note <- c(note, sprintf(
        "%s: congeneric fit %s; stepping down to tau-equivalent (Heywood cases are never fixed to zero)",
        label, if (!fit$conv) "did not converge" else "hit a Heywood bound"))
      used <- "tau"
      fit <- .med_fa1(S, n, model = "tau")
    }
    ## a zero-loading fit is as inadmissible as a Heywood case: the modules
    ## share no common variance, the Bartlett denominator collapses, and the
    ## error-variance offset explodes (observed on real data: omega = 0,
    ## V ~ 1e16 poisoning the joint moment matrix)
    degenerate <- fit$conv && !fit$heywood && {
      bl0 <- .med_bartlett(Z, fit$lambda, fit$theta)
      bl0$omega < 0.30 || !is.finite(bl0$V) || bl0$V > 1e6
    }
    if (degenerate) {
      note <- c(note, sprintf(
        "%s: %s fit empirically underidentified (omega %.2f < .30, the weak-factor boundary - the modules do not cohere enough to define a common factor; Cox 2022, Rhemtulla et al. 2020); reliability-corrected single indicator",
        label, used, {bl0 <- .med_bartlett(Z, fit$lambda, fit$theta); bl0$omega}))
    } else if (fit$conv && !fit$heywood) {
      bl <- .med_bartlett(Z, fit$lambda, fit$theta)
      ## unit-LATENT-variance conditioning: rescale so Var(eta) = 1 -
      ## coefficients are then per SD of the latent tone, and no system's
      ## scale can dominate the joint moment matrix; V transforms as a
      ## variance (guarded away from zero)
      s2 <- max(stats::var(bl$f) - bl$V, 0.05 * stats::var(bl$f))
      if (is.finite(s2) && s2 > 0) { bl$f <- bl$f / sqrt(s2); bl$V <- bl$V / s2 }
      return(list(f = bl$f, V = bl$V, omega = bl$omega,
                  determinacy = bl$determinacy, lambda = fit$lambda,
                  theta = fit$theta, model = used, p = p, notes = note))
    }
    if (!degenerate) note <- c(note, sprintf(
      "%s: %s fit inadmissible; reliability-corrected single indicator", label, used))
    used <- "single"
  }
  ## single indicator: composite mean tone. Reliability is either FIXED a
  ## priori via single_rel (Savalei 2019: fixed values outperform same-sample
  ## estimates in small samples - the registered PoC choice) or estimated by
  ## Spearman-Brown over the mean inter-tone correlation.
  ## single-indicator composite: unweighted mean of the standardized
  ## module tones by default; mod_w (e.g. probe-count precision weights)
  ## makes it a weighted mean. Within-module probe weighting is upstream
  ## in dmsa_scores() either way.
  comp <- if (is.null(mod_w) || length(mod_w) != p) rowMeans(tones) else
    as.numeric(tones %*% (mod_w / sum(mod_w)))
  if (!is.null(single_rel)) {
    rel <- min(max(as.numeric(single_rel), 0.05), 0.99)
  } else if (p >= 2L) {
    R <- stats::cor(tones); rbar <- mean(R[lower.tri(R)])
    rel <- (p * rbar) / (1 + (p - 1) * rbar)
    rel <- min(max(rel, 0.05), 0.99)
  } else rel <- 0.8                          # declared, not estimated: p = 1
  f <- comp - mean(comp)
  ## same unit-LATENT-variance conditioning as the factor branch:
  ## latent variance = rel * Var(composite), so scaling by its sqrt gives
  ## Var(eta) = 1 and V = (1 - rel) / rel exactly
  s2 <- rel * stats::var(f)
  if (is.finite(s2) && s2 > 0) f <- f / sqrt(s2)
  V <- (1 - rel) / rel
  list(f = f, V = V, omega = rel, determinacy = NA_real_,
       lambda = rep(NA_real_, p), theta = rep(NA_real_, p),
       model = "single", p = p, notes = note)
}

## ---------------------------------------------------------------------------
## corrected moments (SAM/Croon with Bogaert-Loh-Rosseel small-sample guard)
## ---------------------------------------------------------------------------

## Build the corrected moment matrix for columns of D (n x k, centered).
## voffset: named vector of measurement-error variances, one entry per LATENT
## column of D (zero for observed columns). products: list of c(a, b) column
## names for each product column, named by the product column.
## Offsets (Bartlett scores; disjoint indicator sets; centered latents):
##   Var(f_a)            = Sig_aa + V_a                       -> -V_a
##   Var(f_a f_b)        = Var(eta_a eta_b) + Sig_aa V_b
##                         + V_a Sig_bb + V_a V_b             -> -(...)
##   Cov(f_a f_b, f_c f_b) (shared b) = ... + Sig_ac V_b      -> -Sig_ac V_b
##   everything else unbiased.
## alpha-correction (Bogaert et al. 2022): subtract only (1 - alpha/(n-1))
## of the offset. lambda-correction: shrink the offset until the corrected
## matrix is positive definite; the shrink factor is reported.
.med_moments <- function(D, voffset, products, n, alpha = 0) {
  k <- ncol(D); nm <- colnames(D)
  Sobs <- crossprod(D) / (n - 1)
  O <- matrix(0, k, k, dimnames = list(nm, nm))
  V <- stats::setNames(numeric(k), nm)
  V[names(voffset)] <- voffset
  ## first-order latent variances
  for (a in names(voffset)) O[a, a] <- voffset[[a]]
  ## corrected first-order covariances needed inside product offsets
  sig <- function(a, b) Sobs[a, b] - if (a == b) V[[a]] else 0
  for (pn in names(products)) {
    ab <- products[[pn]]; a <- ab[1]; b <- ab[2]
    O[pn, pn] <- sig(a, a) * V[[b]] + V[[a]] * sig(b, b) + V[[a]] * V[[b]]
    for (qn in setdiff(names(products), pn)) {
      cd <- products[[qn]]
      shared <- intersect(ab, cd)
      if (length(shared) == 1L) {
        oth <- c(setdiff(ab, shared), setdiff(cd, shared))
        O[pn, qn] <- O[qn, pn] <- sig(oth[1], oth[2]) * V[[shared]]
      }
    }
  }
  tf <- 1 - alpha / (n - 1)                 # alpha-correction
  shrink <- 1
  repeat {
    Sstar <- Sobs - (tf * shrink) * O
    ev <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)$values
    ok_pd <- min(ev) > 1e-6 * max(ev)
    ## CONDITIONING cap (extended lambda-correction, Bogaert spirit): a
    ## corrected matrix can be PD yet so ill-conditioned that coefficients
    ## explode - observed on real data at n ~ 80-110 with rel .70 product
    ## offsets (18x coefficient inflation, degenerate p = 1 cells). Shrink
    ## the offset until the corrected CORRELATION matrix has condition
    ## number <= 1e3; the correction backs off only where it cannot be
    ## afforded, and the shrink factor is reported.
    ok_cond <- ok_pd && {
      evR <- eigen(stats::cov2cor(Sstar), symmetric = TRUE,
                   only.values = TRUE)$values
      min(evR) > 0 && max(evR) / min(evR) <= 1e3
    }
    if (ok_cond || shrink < 1e-3) break
    shrink <- shrink * 0.9                  # lambda-correction
  }
  list(S = Sstar, Sobs = Sobs, O = O, shrink = shrink * tf,
       lambda_triggered = shrink < 1)
}

## Collinearity diagnostics on a (corrected) moment matrix: VIF per
## predictor and the condition number of the predictor correlation matrix.
## NOTE the corrected metric is the right place to look: measurement-error
## correction DISATTENUATES latent correlations, so collinearity is worse
## here than among the raw scores (e.g. r_obs = .5 between two .70-
## reliability tones is r_latent ~ .71).
.med_collin <- function(S, xn) {
  R <- tryCatch(stats::cov2cor(S[xn, xn, drop = FALSE]),
                error = function(e) NULL)
  if (is.null(R)) return(list(vif = NULL, kappa = NA_real_))
  Ri <- tryCatch(solve(R), error = function(e) NULL)
  vif <- if (is.null(Ri)) stats::setNames(rep(Inf, length(xn)), xn) else
    stats::setNames(diag(Ri), xn)
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  list(vif = vif, kappa = sqrt(max(ev) / max(min(ev), 1e-12)))
}

## Moment-space OLS with t statistics. S: corrected moment matrix.
.med_reg <- function(S, xn, yn, n) {
  Sxx <- S[xn, xn, drop = FALSE]; Sxy <- S[xn, yn]
  Sxxi <- tryCatch(solve(Sxx), error = function(e)
    solve(Sxx + diag(1e-6 * mean(diag(Sxx)), nrow(Sxx))))
  b <- as.numeric(Sxxi %*% Sxy); names(b) <- xn
  rss <- max(as.numeric(S[yn, yn] - crossprod(b, Sxy)), 1e-12)
  kk <- length(xn) + 1L                      # + intercept (data are centered)
  se <- sqrt(rss * diag(Sxxi) / (n - kk))
  list(b = b, se = se, t = b / se, rss = rss, Sxxi = Sxxi)
}

## Freedman-Lane permutation p for one term of a moment-space regression.
## D: centered data matrix (columns include xn and, when y is latent, yn).
## Vy: kept in the signature for the callers; the sandwich pivot (E14) is
## built from the UNCORRECTED response variance and does not use it.
## blocks + block_mode:
##   "within"   residuals permute WITHIN each block (chip scheme, Winkler)
##   "exchange" whole blocks of EQUAL SIZE exchange with each other, row
##              order preserved inside the block (the family/cID scheme -
##              the same relabelling dmsa_model() uses, preserving the
##              within-family dependence under the null)
.med_fl <- function(D, S, xn, yn, term, n, B, blocks = NULL, Vy = 0,
                    block_mode = c("within", "exchange"), Sobs = NULL,
                    pivot = c("hc", "homoscedastic"), strata = NULL) {
  block_mode <- match.arg(block_mode)
  pivot <- match.arg(pivot)
  ## PIVOT (E14, 2026-09-02, PI ruling): the errors-in-variables SANDWICH.
  ## The coefficient is the corrected-moment estimate b = S*xx^-1 Sxy; its
  ## variance, treating the error variances as known, is
  ##   Var(b) = sigma2_obs * S*xx^-1 S_obs S*xx^-1 / (n - k),
  ##   sigma2_obs = Var(y) - 2 b'Sxy + b' S_obs b
  ## (Meijer, Spierdijk & Wansbeek 2020; Devlieger, Mayer & Rosseel 2016),
  ## where S_obs is the UNCORRECTED moment matrix of the regressors. The
  ## previous pivot divided the corrected coefficient by the plain-OLS
  ## residual (Var(y) - b'Sxy), which over-subtracts once b is
  ## disattenuated: level was exact (same pivot on observed and permuted
  ## data) but the null had heavy tails and power on error-prone
  ## main-effect paths dropped (b1 .735 -> .815, a1 .855 -> .945 at n = 200,
  ## rel .70, identical permutation; interaction paths unchanged). Vy is no
  ## longer needed by the pivot: sigma2_obs is built from the uncorrected
  ## response variance by construction. Freedman-Lane geometry is
  ## unchanged: reduced-model residuals from the RAW OLS fit (orthogonal to
  ## the reduced design in-sample), the corrected S*xx and the sandwich
  ## matrix held FIXED across permutations, the identical statistic on
  ## observed and permuted responses.
  ## PIVOT (E15, 2026-09-03, PI ruling): the HETEROSCEDASTICITY-CONSISTENT
  ## errors-in-variables sandwich is the default. The corrected estimator
  ## solves  sum_i x_i (y_i - x_i'b)/(n-1) + Ofs b = 0  with
  ## Ofs = S_obs,xx - S*xx (the moment offsets), so its score for person i is
  ##   s_i = x_i e_i + Ofs b,   Var(b) = S*xx^-1 [sum_i s_i s_i'] S*xx^-1/(n-1)^2
  ## (HC0 form; Buonaccorsi 2010 ch. 5), and the studentised Freedman-Lane
  ## of DiCiccio & Romano (2017) recomputes the meat on every permuted
  ## response. With no measurement error and homoscedastic disturbances it
  ## is the ordinary HC0 statistic. GATE_REPORT_hc: under a heteroscedastic
  ## mediator disturbance the interaction test moved from .087-.097 to .053;
  ## no cost under homoscedasticity; -6 points of interaction power when the
  ## disturbance IS heteroscedastic. pivot = "homoscedastic" reproduces the
  ## E14 statistic (0.99.11-0.99.12) exactly.
  Sraw <- if (is.null(Sobs)) S else Sobs
  Sxx <- S[xn, xn, drop = FALSE]
  Sxxi <- tryCatch(solve(Sxx), error = function(e)
    solve(Sxx + diag(1e-6 * mean(diag(Sxx)), nrow(Sxx))))
  So <- Sraw[xn, xn, drop = FALSE]
  A <- Sxxi %*% So %*% Sxxi                 # sandwich core, fixed
  Ofs <- So - Sxx                            # moment offsets on the regressors
  kk <- length(xn) + 1L                      # + intercept (data are centered)
  j <- match(term, xn)
  Xc <- D[, xn, drop = FALSE]
  y <- D[, yn]
  vy_obs <- as.numeric(sum(y^2) / (n - 1))   # uncorrected response variance
  stat <- if (pivot == "homoscedastic") function(yv, vyv) {
    Sxy <- as.numeric(crossprod(Xc, yv)) / (n - 1)
    b <- as.numeric(Sxxi %*% Sxy)
    s2 <- vyv - 2 * sum(b * Sxy) + as.numeric(crossprod(b, So %*% b))
    bad <- !is.finite(s2) || s2 <= 0
    if (bad) s2 <- 1e-12
    list(b = b, se = sqrt(s2 * A[j, j] / (n - kk)), bad = bad)
  } else function(yv, vyv) {
    Sxy <- as.numeric(crossprod(Xc, yv)) / (n - 1)
    b <- as.numeric(Sxxi %*% Sxy)
    e <- as.numeric(yv - Xc %*% b)
    Sc <- Xc * e + matrix(as.numeric(Ofs %*% b), n, length(xn), byrow = TRUE)
    V <- Sxxi %*% crossprod(Sc) %*% Sxxi / (n - 1)^2
    bad <- !is.finite(V[j, j]) || V[j, j] <= 0
    list(b = b, se = sqrt(if (bad) 1e-12 else V[j, j]), bad = bad)
  }
  obs <- stat(y, vy_obs)
  tobs <- obs$b[j] / obs$se
  red <- setdiff(xn, term)
  ## E16 (2026-09-03, PI ruling): an EMPTY reduced model (one predictor, no
  ## covariates) has fitted values 0 - the residuals are the centered response
  fit <- if (length(red)) {
    rf <- .med_reg(Sraw, red, yn, n)
    as.numeric(D[, red, drop = FALSE] %*% rf$b)
  } else rep(0, n)
  r <- y - fit
  rows <- if (is.null(blocks)) list(seq_len(n)) else
    split(seq_len(n), blocks)
  if (is.null(blocks) && !is.null(strata))      # strata without blocks: permute within stratum
    rows <- split(seq_len(n), as.character(strata))
  ## exchange groups: blocks of equal size, and - when `strata` is given -
  ## of the same stratum (E15 ruling: the after-wave of the long builds;
  ## GATE_REPORT_delta section 1). `strata` is a per-row factor; a block
  ## must lie inside one stratum.
  grp <- if (block_mode == "exchange") {
    key <- lengths(rows)
    if (!is.null(strata)) {
      sv <- as.character(strata)
      bs <- vapply(rows, function(i) {
        u <- unique(sv[i])
        if (length(u) != 1L) stop("a permutation block spans more than one stratum", call. = FALSE)
        u }, character(1))
      key <- paste(key, bs)
    }
    split(seq_along(rows), key)
  } else NULL
  tn <- numeric(B); bad <- 0L
  for (bb in seq_len(B)) {
    ix <- integer(n)
    if (block_mode == "exchange" && !is.null(blocks)) {
      for (st in grp) ix[unlist(rows[st], use.names = FALSE)] <-
        unlist(rows[st[sample.int(length(st))]], use.names = FALSE)
    } else {
      for (rw in rows) ix[rw] <- rw[sample.int(length(rw))]
    }
    yp <- fit + r[ix]
    yp <- yp - mean(yp)
    sp <- stat(yp, sum(yp^2) / (n - 1))
    if (sp$bad) bad <- bad + 1L
    tn[bb] <- sp$b[j] / sp$se
  }
  b_all <- obs$b; names(b_all) <- xn
  list(t = tobs, p = (1 + sum(abs(tn) >= abs(tobs))) / (B + 1),
       b = b_all[[term]], se = obs$se, degenerate = bad)
}

## ---------------------------------------------------------------------------
## the core: latent moderated mediation on a prepared tone list
## ---------------------------------------------------------------------------

## tones: named list of n x p_s module-tone matrices (x systems, then m, w).
## dat: data.frame with outcome + covariates (already complete-case, and
## chip-transformed when chip is declared). Everything downstream of the
## map/alignment plumbing lives here so it can be validated standalone.
.med_core <- function(tones, dat, x, m, w, outcome, covariates,
                      measurement, B, blocks, model8, conditional_at,
                      alpha_ssc, boot = 0L, single_rel = NULL,
                      block_mode = "within", mod_weights = NULL,
                      pivot = "hc", strata = NULL) {
  n <- nrow(dat)
  systems <- c(x, m, if (!is.null(w)) w)
  meas <- lapply(systems, function(s)
    .med_measure(tones[[s]], measurement = measurement, label = s,
                 single_rel = single_rel, mod_w = mod_weights[[s]]))
  names(meas) <- systems
  F <- vapply(meas, `[[`, numeric(n), "f")
  colnames(F) <- systems
  voff <- vapply(meas, `[[`, numeric(1), "V")

  ## assemble the centered data matrix: latents, products, covariates, ys
  prods <- list()
  D <- F
  if (!is.null(w)) for (xx in x) {
    pn <- paste0(xx, ":", w)
    D <- cbind(D, F[, xx] * F[, w]); colnames(D)[ncol(D)] <- pn
    prods[[pn]] <- c(xx, w)
  }
  C <- if (length(covariates))
    scale(as.matrix(dat[, covariates, drop = FALSE]), scale = FALSE) else NULL
  yv <- dat[[outcome]] - mean(dat[[outcome]])
  D <- cbind(D, C, y = yv)
  if (!is.null(C)) colnames(D)[seq(ncol(F) + length(prods) + 1L,
                                   length.out = ncol(C))] <- covariates
  D <- scale(D, scale = FALSE)

  mm <- .med_moments(D, voffset = voff, products = prods, n = n,
                     alpha = alpha_ssc)

  ## M-model: m ~ x + w + x:w + covariates
  a_terms <- c(x, if (!is.null(w)) w, names(prods))
  mx <- c(a_terms, covariates)
  ## Y-model: y ~ m + x + w (+ products if model 8) + covariates
  yx <- c(m, x, if (!is.null(w)) w, if (model8) names(prods), covariates)

  paths <- list()
  for (tm in a_terms) paths[[paste0("a[", tm, "]")]] <-
    .med_fl(D, mm$S, mx, m, tm, n, B, blocks, Vy = voff[[m]],
            block_mode = block_mode, Sobs = mm$Sobs,
            pivot = pivot, strata = strata)
  paths[[paste0("b[", m, "]")]] <-
    .med_fl(D, mm$S, yx, "y", m, n, B, blocks, Vy = 0,
            block_mode = block_mode, Sobs = mm$Sobs,
            pivot = pivot, strata = strata)
  ## direct effects: the X (and W) coefficients of the SAME Y-model - the
  ## c' paths of the mediation triad, estimated controlling for M
  for (tm in c(x, if (!is.null(w)) w))
    paths[[paste0("c'[", tm, "]")]] <-
      .med_fl(D, mm$S, yx, "y", tm, n, B, blocks, Vy = 0,
              block_mode = block_mode, Sobs = mm$Sobs,
            pivot = pivot, strata = strata)
  if (model8) for (pn in names(prods)) paths[[paste0("c[", pn, "]")]] <-
    .med_fl(D, mm$S, yx, "y", pn, n, B, blocks, Vy = 0,
            block_mode = block_mode, Sobs = mm$Sobs,
            pivot = pivot, strata = strata)

  collin <- list(M_model = .med_collin(mm$S, mx),
                 Y_model = .med_collin(mm$S, yx))

  bkey <- paste0("b[", m, "]")
  idx <- list()
  if (!is.null(w)) {
    for (xx in x) {
      pn <- paste0(xx, ":", w); ak <- paste0("a[", pn, "]")
      idx[[paste0("IMM[", xx, "]")]] <- list(
        est = paths[[ak]]$b * paths[[bkey]]$b,
        p_js = max(paths[[ak]]$p, paths[[bkey]]$p),
        parts = c(ak, bkey))
    }
  } else {
    for (xx in x) {
      ak <- paste0("a[", xx, "]")
      idx[[paste0("IE[", xx, "]")]] <- list(
        est = paths[[ak]]$b * paths[[bkey]]$b,
        p_js = max(paths[[ak]]$p, paths[[bkey]]$p),
        parts = c(ak, bkey))
    }
  }
  pj <- vapply(idx, `[[`, numeric(1), "p_js")
  ph <- stats::p.adjust(pj, method = "holm")
  for (i in seq_along(idx)) idx[[i]]$p_holm <- ph[i]

  cond <- NULL
  if (!is.null(w)) {
    sdw <- sqrt(max(mm$S[w, w], 0))
    cond <- do.call(rbind, lapply(x, function(xx) {
      pn <- paste0(xx, ":", w)
      a1 <- paths[[paste0("a[", xx, "]")]]$b
      a3 <- paths[[paste0("a[", pn, "]")]]$b
      b1 <- paths[[bkey]]$b
      data.frame(x = xx, w_sd = conditional_at,
                 indirect = (a1 + a3 * conditional_at * sdw) * b1)
    }))
  }

  list(paths = paths, indices = idx, conditional = cond, measurement = meas,
       moments = mm, collinearity = collin, n = n, D = D,
       terms = list(mx = mx, yx = yx, a = a_terms, bkey = bkey))
}

## ---------------------------------------------------------------------------
## the user-facing engine
## ---------------------------------------------------------------------------

#' Latent moderated mediation over DMSA system tones
#'
#' Tests a moderated-mediation claim in which the predictors, moderator and
#' mediator are latent SYSTEM TONES, each measured by its module tones as
#' reflective indicators, and the outcome is an observed variable. The
#' structural model runs on measurement-error-corrected moments
#' (structural-after-measurement / Croon family), so the interaction paths
#' are not attenuated by tone unreliability; inference is permutation-native
#' (Freedman-Lane per path, joint-significance for products, Holm over the
#' declared index family).
#'
#' @param data data.frame of phenotype rows, or a \code{dmsa_import} object.
#' @param methylation Optional samples x probes matrix aligned to
#'   \code{data} rows (beta scale by default). When \code{NULL}, probe
#'   columns are taken from \code{data} itself.
#' @param map A \code{dmsa_sets()} cascade or anything it accepts
#'   (default \code{"alpha"}).
#' @param x Character vector of predictor system short names (1+).
#' @param m Mediator system short name.
#' @param w Moderator system short name, or \code{NULL} for plain mediation.
#' @param outcome Outcome column in \code{data}.
#' @param covariates Character vector of covariate columns.
#' @param measurement Measurement model for each system's module tones:
#'   \code{"tau"} (essentially tau-equivalent - equal loadings, free
#'   uniquenesses; the small-sample default per Strauss 2026),
#'   \code{"congeneric"} (free loadings; steps down to tau on
#'   non-convergence or a Heywood bound - bounds are never patched to
#'   zero), or \code{"single"} (reliability-corrected single indicator,
#'   Savalei 2019). Systems with fewer than 3 usable modules fall back to
#'   \code{"single"} automatically.
#' @param min_module_probes A module enters as an indicator only with at
#'   least this many usable aligned probes (default 3).
#' @param single_rel Optional FIXED reliability for any system that lands in
#'   single-indicator mode (Savalei 2019: fixed a-priori values outperform
#'   same-sample estimates in small samples). \code{NULL} estimates by
#'   Spearman-Brown instead.
#' @param module_weights How module tones combine in the SINGLE-INDICATOR
#'   composite: \code{"equal"} (default - standardized tones, equal a
#'   priori weight), \code{"sqrt_probes"} or \code{"probes"} (precision
#'   weights by each module's usable probe count). Factor rungs always
#'   weight modules by lambda/theta (Bartlett); probe-level weighting
#'   inside each module tone is upstream and always on.
#' @param measure_covariates Which covariates are residualized out of the
#'   module tones BEFORE the measurement stage (the composition
#'   pre-adjustment). \code{NULL} (default) uses \code{covariates};
#'   \code{character(0)} disables pre-adjustment (composition then acts
#'   only through the structural covariates); a subset (e.g. the cell
#'   fractions alone) brackets the under-/over-adjustment trade-off.
#' @param chip Optional column in \code{data} naming the array chip. When
#'   given, module tones are chip-quasi-demeaned (REML, R/ranef.R) BEFORE
#'   the measurement stage, and outcome/covariates before the structural
#'   stage.
#' @param block Optional column in \code{data} naming the exchangeability
#'   block (the family/couple identifier, cID). When given, Freedman-Lane
#'   relabels WHOLE blocks of equal size (row order preserved inside each
#'   block - the same relabelling \code{dmsa_model()} uses), preserving
#'   within-family dependence under the null, and the bootstrap resamples
#'   blocks rather than individuals. Incompatible with
#'   \code{perm_blocks = "chip"}.
#' @param perm_blocks \code{"none"} permutes residuals unrestricted after
#'   the chip transform; \code{"chip"} permutes within chip blocks
#'   (Winkler et al. 2015). The registered simulation gate selects the
#'   scheme; both are available so the gate can be run.
#' @param model \code{"7"} (moderation on the a-paths only) or \code{"8"}
#'   (also on the direct paths).
#' @param map_level \code{"full"} or \code{"confidence"} - which analysed
#'   direction-map level supplies the probe universe (as in
#'   \code{dmsa_frame()}); module membership joins from the cascade by gene.
#' @param B Permutations per path.
#' @param boot Bootstrap resamples for percentile CIs; the ENTIRE pipeline
#'   (measurement models included) is refit in every resample. \code{0}
#'   disables.
#' @param alpha_ssc Bogaert-Loh-Rosseel alpha for the small-sample offset
#'   correction; \code{0} (default) subtracts the full offset, with the
#'   lambda (positive-definiteness) shrink always active and reported.
#' @param conditional_at Moderator values (in latent SD units) at which
#'   conditional indirect effects are reported.
#' @param beta_input \code{TRUE} when \code{methylation} holds betas
#'   (converted to M-values internally).
#' @param seed Optional integer; the caller's RNG state is restored.
#' @param pivot Studentisation of the Freedman-Lane statistic:
#'   \code{"hc"} (default since 0.99.13, engine ruling E15) is the
#'   heteroscedasticity-consistent errors-in-variables sandwich - the
#'   estimating-equation sandwich of the corrected estimator with the meat
#'   recomputed per permutation (DiCiccio & Romano 2017; Buonaccorsi 2010);
#'   \code{"homoscedastic"} is the E14 sandwich of 0.99.11-0.99.12 and
#'   reproduces those results exactly. Under a heteroscedastic mediator
#'   disturbance the E14 interaction test ran at .09-.10; the HC pivot holds
#'   .05 at no cost under homoscedasticity (GATE_REPORT_hc).
#' @param strata Optional column in \code{data}: permutation stratum. With
#'   \code{block}, equal-size blocks exchange only within a stratum (each
#'   block must lie inside one stratum); without \code{block}, residuals
#'   permute within stratum. The registered use is the after-wave of the
#'   before/after builds.
#' @return An object of class \code{dmsa_mediate}: per-system measurement
#'   tables, path coefficient tables with permutation p per path, indices
#'   of moderated mediation with joint-significance + Holm p, conditional
#'   indirect effects, bootstrap CIs when requested, and the corrected
#'   moment machinery for inspection.
#' @examples
#' ## synthetic end-to-end run at analysis scale (see tests for the
#' ## simulation-gate harness)
#' \donttest{
#' set.seed(1)
#' n <- 120
#' eta <- matrix(rnorm(n * 3), n, 3)                 # hpa, oxytocin, kyn
#' mk <- function(f, p) sapply(seq_len(p), function(j)
#'   0.7 * f + rnorm(n, sd = sqrt(1 - 0.49)))
#' tones <- list(hpa = mk(eta[, 1], 6), oxytocin = mk(eta[, 2], 8),
#'               kyn = mk(0.4 * eta[, 1] * eta[, 2] + 0.3 * eta[, 1], 5))
#' dat <- data.frame(y = 0.5 * (0.4 * eta[, 1] * eta[, 2] + 0.3 * eta[, 1]) +
#'                     rnorm(n), age = rnorm(n))
#' fit <- dmsa:::.med_core(tones, dat, x = "hpa", m = "kyn", w = "oxytocin",
#'                         outcome = "y", covariates = "age",
#'                         measurement = "tau", B = 199, blocks = NULL,
#'                         model8 = FALSE, conditional_at = c(-1, 0, 1),
#'                         alpha_ssc = 0)
#' fit$indices
#' }
#' @export
dmsa_mediate <- function(data, methylation = NULL, map = "alpha",
                         x, m, w = NULL, outcome, covariates = character(),
                         measurement = c("tau", "congeneric", "single"),
                         min_module_probes = 3L, single_rel = NULL,
                         module_weights = c("equal", "sqrt_probes", "probes"),
                         measure_covariates = NULL,
                         chip = NULL, block = NULL,
                         perm_blocks = c("none", "chip"),
                         model = c("7", "8"),
                         map_level = c("full", "confidence"),
                         B = 1999L, boot = 0L, alpha_ssc = 0,
                         conditional_at = c(-1, 0, 1),
                         beta_input = TRUE, seed = NULL,
                         pivot = c("hc", "homoscedastic"), strata = NULL) {
  measurement <- match.arg(measurement)
  perm_blocks <- match.arg(perm_blocks)
  pivot <- match.arg(pivot)
  model <- match.arg(model)
  map_level <- match.arg(map_level)
  module_weights <- match.arg(module_weights)
  if (!is.null(seed)) {
    .old_seed <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(.old_seed))
      assign(".Random.seed", .old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv()))
        rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(seed)
  }

  ## -- unwrap a dmsa_import delivery, dmsa_frame-style ---------------------
  if (inherits(data, "dmsa_import")) {
    if (is.null(methylation)) methylation <- data$methylation
    data <- data$data
  }
  data <- as.data.frame(data)
  sets <- if (inherits(map, "dmsa_sets")) map else dmsa_sets(map)
  cas <- sets$cascade

  ## E8 fix (2026-08-29, PI-approved): REFUSE rather than silently substitute.
  ## A dmsa_sets object built from a CUSTOM cascade, or a cascade CSV path,
  ## used to fall through to the bundled Alpha coverage map - module
  ## membership from the user's biology, probe directions from Alpha's panel,
  ## and nothing on screen saying so. Direction must now be explicit.
  .map_ok <- identical(map, "alpha") ||
    (is.data.frame(map) &&
       all(c("probe", "column", "gene", "system_id", "best_direction")
           %in% names(map))) ||
    (inherits(map, "dmsa_sets") &&
       identical(attr(map, "source") %||% map$source %||% "",
                 .cas_builtin("cascade") %||% "__none__"))
  if (!.map_ok)
    stop("`map` carries module membership but no probe directions, and dmsa ",
         "will not silently substitute the bundled Alpha direction map for ",
         "a custom cascade.\nEither pass map = \"alpha\" to use the bundled ",
         "directions EXPLICITLY alongside your cascade, or supply a ",
         "direction-map data.frame (columns probe, column, gene, system_id, ",
         "best_direction, p_plus).", call. = FALSE)


  systems <- c(x, m, w)
  bad <- setdiff(systems, unique(cas$system_short))
  if (length(bad))
    stop("system short name(s) not in the map: ", paste(bad, collapse = ", "),
         "\nsee dmsa_systems()", call. = FALSE)
  if (m %in% c(x, w) || (!is.null(w) && any(x == w)))
    stop("x, m and w must be distinct systems", call. = FALSE)

  ## -- complete cases over the model variables -----------------------------
  if (!is.null(block) && perm_blocks == "chip")
    stop("`block` (family exchangeability) and perm_blocks = 'chip' cannot ",
         "be combined; family blocks are the exchangeable unit", call. = FALSE)
  need <- c(outcome, covariates, chip, block, strata)
  miss <- setdiff(need, names(data))
  if (length(miss))
    stop("column(s) not in data: ", paste(miss, collapse = ", "), call. = FALSE)
  keep <- stats::complete.cases(data[, need, drop = FALSE])

  ## -- the analysed direction map (dmsa_frame's authority) ------------------
  maps <- .frame_maps(if (inherits(map, "dmsa_sets")) "alpha" else map)
  mp <- maps[[map_level]]
  mp$system_id <- as.character(mp$system_id)

  ## -- methylation matrix, M-value scale -----------------------------------
  ME <- if (is.null(methylation)) {
    pc <- intersect(unique(c(mp$probe, mp$column)), names(data))
    if (!length(pc)) stop("no map probes found among data columns; pass ",
                          "`methylation`", call. = FALSE)
    as.matrix(data[, pc, drop = FALSE])
  } else as.matrix(methylation)
  if (nrow(ME) != nrow(data))
    stop("methylation has ", nrow(ME), " rows but data has ", nrow(data),
         call. = FALSE)
  mode(ME) <- "numeric"
  if (beta_input) ME <- dmsa_mvalues(ME)
  data <- data[keep, , drop = FALSE]; ME <- ME[keep, , drop = FALSE]
  n <- nrow(data)
  if (n < 30L) stop("only ", n, " complete rows; this engine is not built ",
                    "for n < 30", call. = FALSE)

  ## -- module tones per system ---------------------------------------------
  ## Probe universe and directions come from the ANALYSED map (dmsa_frame's
  ## authority, full or confidence level); module membership comes from the
  ## selection cascade, joined by gene.
  chipv <- if (!is.null(chip)) as.factor(data[[chip]]) else NULL
  sysid <- stats::setNames(as.character(sets$systems$system_id),
                           sets$systems$system_short)
  tones <- list(); coverage <- list()
  ## polarity table (gene -> w_g), subset per system before aligning: gene
  ## symbols can repeat across systems with different polarity
  pol_tab <- if (inherits(sets$polarity, "dmsa_polarity"))
    sets$polarity$polarity else sets$polarity
  for (s in systems) {
    sm <- mp[mp$system_id == sysid[[s]], , drop = FALSE]
    ps <- if (is.null(pol_tab)) NULL else {
      pp <- if ("system_id" %in% names(pol_tab))
        pol_tab[as.character(pol_tab$system_id) == sysid[[s]], , drop = FALSE]
      else pol_tab
      if (nrow(pp)) pp else NULL
    }
    gm <- unique(cas[cas$system_short == s,
                     c("module_id", "module", "gene"), drop = FALSE])
    sm$module_id <- gm$module_id[match(sm$gene, gm$gene)]
    sm$module <- gm$module[match(sm$gene, gm$gene)]
    ## a probe can appear under its bare CpG id or the manifest-suffixed
    ## analysis column name; match on whichever key hits
    sm$col <- ifelse(sm$probe %in% colnames(ME), sm$probe,
                     ifelse(sm$column %in% colnames(ME), sm$column,
                            NA_character_))
    sm <- sm[!is.na(sm$col) & !is.na(sm$module_id) & !duplicated(sm$col), ,
             drop = FALSE]
    tl <- list(); cov_rows <- list()
    for (md in unique(gm$module_id)) {
      cm <- sm[sm$module_id == md, , drop = FALSE]
      lab <- gm$module[match(md, gm$module_id)]
      if (nrow(cm) < min_module_probes) {
        cov_rows[[md]] <- data.frame(system = s, module = lab,
                                     probes = nrow(cm), used = FALSE)
        next
      }
      al <- if (is.null(ps)) NULL else
        tryCatch(suppressMessages(suppressWarnings(
          dmsa_align(data.frame(cpg = cm$probe, d = cm$best_direction,
                                p_plus = cm$p_plus),
                     genes = cm$gene, level = "system",
                     polarity = ps, missing_polarity = "drop"))),
          error = function(e) NULL)
      sc <- if (is.null(al) || sum(al$usable) < min_module_probes) NULL else
        tryCatch(suppressWarnings(
          dmsa_scores(ME[, cm$col[match(al$probe, cm$probe)], drop = FALSE],
                      al, flavours = "aligned")$aligned),
          error = function(e) NULL)
      if (!is.null(sc)) { tl[[lab]] <- sc; attr(tl, "np") <-
        c(attr(tl, "np"), nrow(cm)) }
      cov_rows[[md]] <- data.frame(system = s, module = lab,
                                   probes = nrow(cm), used = !is.null(sc))
    }
    if (length(tl) < 1L)
      stop("system '", s, "' has no usable module tone at ",
           min_module_probes, "+ probes under the '", map_level,
           "' map", call. = FALSE)
    np <- attr(tl, "np")
    mw <- switch(module_weights, equal = NULL, probes = np,
                 sqrt_probes = sqrt(np))
    TT <- do.call(cbind, tl)
    attr(TT, "mod_w") <- mw
    ## chip quasi-demeaning BEFORE the measurement stage
    if (!is.null(chipv)) {
      nc <- as.numeric(table(chipv)); names(nc) <- levels(chipv)
      gi <- as.character(chipv)
      X0 <- matrix(1, n, 1)
      for (j in seq_len(ncol(TT))) {
        gam <- tryCatch(.ri_reml(TT[, j], X0, chipv), error = function(e) 0)
        if (gam > 0) TT[, j] <- as.numeric(
          .ri_transform(TT[, j, drop = FALSE], chipv, nc, gi, gam))
      }
    }
    ## composition pre-adjustment: residualize module tones on the model
    ## covariates BEFORE the measurement stage (symmetric with the chip
    ## transform). In whole blood, cell-composition differences manufacture
    ## spurious between-gene covariance (Jaffe & Irizarry 2014, Genome
    ## Biology), so a latent tone fit on unadjusted module tones can be a
    ## cell-fraction factor in disguise.
    mcv <- if (is.null(measure_covariates)) covariates else measure_covariates
    if (length(mcv)) {
      CV <- cbind(1, as.matrix(data[, mcv, drop = FALSE]))
      mode(CV) <- "numeric"
      if (all(is.finite(CV))) {
        qrc <- qr(CV)
        mw_keep <- attr(TT, "mod_w")
        TT <- apply(TT, 2, function(v) qr.resid(qrc, v))
        attr(TT, "mod_w") <- mw_keep
      }
    }
    tones[[s]] <- TT
    coverage[[s]] <- do.call(rbind, cov_rows)
  }

  ## -- chip transform for outcome and covariates ---------------------------
  dat2 <- data
  if (!is.null(chipv)) {
    nc <- as.numeric(table(chipv)); names(nc) <- levels(chipv)
    gi <- as.character(chipv); X0 <- matrix(1, n, 1)
    for (v in c(outcome, covariates)) {
      vv <- as.numeric(dat2[[v]])
      gam <- tryCatch(.ri_reml(vv, X0, chipv), error = function(e) 0)
      if (gam > 0) dat2[[v]] <- as.numeric(
        .ri_transform(cbind(vv), chipv, nc, gi, gam))
    }
  }
  blockv <- if (!is.null(block)) as.factor(data[[block]]) else NULL
  blocks <- if (!is.null(blockv)) blockv else if (perm_blocks == "chip") {
    if (is.null(chipv)) stop("perm_blocks = 'chip' needs `chip`", call. = FALSE)
    chipv
  } else NULL
  bmode <- if (!is.null(blockv)) "exchange" else "within" 
  stratav <- if (!is.null(strata)) as.factor(data[[strata]]) else NULL

  mwl <- lapply(tones, attr, "mod_w")
  core <- .med_core(tones, dat2, x = x, m = m, w = w, outcome = outcome,
                    covariates = covariates, measurement = measurement,
                    B = B, blocks = blocks, model8 = model == "8",
                    conditional_at = conditional_at, alpha_ssc = alpha_ssc,
                    single_rel = single_rel, block_mode = bmode,
                    mod_weights = mwl, pivot = pivot, strata = stratav)

  ## -- full-pipeline nonparametric bootstrap (measurement refit per draw) --
  ci <- NULL
  if (boot > 0L) {
    est_of <- function(cc) {
      out <- c(vapply(cc$paths, function(p) p$b, numeric(1)),
               vapply(cc$indices, function(i) i$est, numeric(1)))
      out
    }
    obs <- est_of(core)
    BM <- matrix(NA_real_, boot, length(obs),
                 dimnames = list(NULL, names(obs)))
    brows <- if (!is.null(blockv)) split(seq_len(n), blockv) else NULL
    for (bb in seq_len(boot)) {
      ix <- if (is.null(brows)) sample.int(n, n, replace = TRUE) else
        unlist(brows[sample.int(length(brows), length(brows),
                                replace = TRUE)], use.names = FALSE)
      cc <- tryCatch(
        .med_core(lapply(tones, function(TT) TT[ix, , drop = FALSE]),
                  dat2[ix, , drop = FALSE], x = x, m = m, w = w,
                  outcome = outcome, covariates = covariates,
                  measurement = measurement, B = 0L, blocks = NULL,
                  model8 = model == "8", conditional_at = conditional_at,
                  alpha_ssc = alpha_ssc, single_rel = single_rel,
                  mod_weights = mwl),
        error = function(e) NULL)
      if (!is.null(cc)) {
        eb <- est_of(cc)
        BM[bb, names(eb)] <- eb
      }
    }
    ok <- stats::complete.cases(BM)
    ci <- list(
      lower = apply(BM[ok, , drop = FALSE], 2, stats::quantile, 0.025),
      upper = apply(BM[ok, , drop = FALSE], 2, stats::quantile, 0.975),
      n_ok = sum(ok), draws = boot)
  }

  structure(list(
    call = match.call(), n = n, B = B, systems = systems,
    x = x, m = m, w = w, outcome = outcome, model = model,
    measurement = lapply(core$measurement, function(ms)
      ms[c("model", "p", "V", "omega", "determinacy", "lambda", "theta",
           "notes")]),
    coverage = do.call(rbind, coverage),
    paths = core$paths, indices = core$indices,
    conditional = core$conditional, ci = ci,
    collinearity = core$collinearity,
    perm_blocks = if (!is.null(blockv)) "family-exchange" else perm_blocks,
    pivot = pivot, strata = strata,
    moments = list(shrink = core$moments$shrink,
                   lambda_triggered = core$moments$lambda_triggered)),
    class = "dmsa_mediate")
}

#' @export
print.dmsa_mediate <- function(x, ...) {
  cat("Latent DMSA moderated mediation (Hayes model ", x$model, ")\n", sep = "")
  cat(sprintf("  n = %d   B = %d permutations   blocks: %s   pivot: %s%s\n",
              x$n, x$B, x$perm_blocks, x$pivot %||% "homoscedastic",
              if (!is.null(x$strata)) paste0("   strata: ", x$strata) else ""))
  cat("  measurement (per system):\n")
  for (s in names(x$measurement)) {
    ms <- x$measurement[[s]]
    cat(sprintf("    %-12s %d modules  model %-10s omega %.2f%s%s\n",
                s, ms$p, ms$model, ms$omega,
                if (is.finite(ms$determinacy))
                  sprintf("  determinacy %.2f", ms$determinacy) else "",
                if (length(ms$notes)) "  [stepped]" else ""))
  }
  if (x$moments$lambda_triggered)
    cat(sprintf("  NOTE: offset shrunk to %.2f of full correction to keep the moment matrix PD\n",
                x$moments$shrink))
  if (!is.null(x$collinearity)) {
    cl <- x$collinearity
    for (mdl in names(cl)) {
      v <- cl[[mdl]]$vif
      if (is.null(v)) next
      top <- which.max(v)
      flag <- if (max(v) > 10 || cl[[mdl]]$kappa > 30) "  <-- COLLINEARITY WARNING" else ""
      cat(sprintf("  collinearity %-8s max VIF %.1f (%s)  condition no. %.1f%s\n",
                  mdl, max(v), names(v)[top], cl[[mdl]]$kappa, flag))
    }
  }
  cat("  structural paths (Freedman-Lane p):\n")
  for (nm in names(x$paths)) {
    p <- x$paths[[nm]]
    cat(sprintf("    %-16s b = %+.3f (SE %.3f)  p = %.4f\n",
                nm, p$b, p$se, p$p))
  }
  cat("  indices (joint significance, Holm over the family):\n")
  for (nm in names(x$indices)) {
    i <- x$indices[[nm]]
    cat(sprintf("    %-16s est = %+.4f  p_js = %.4f  p_holm = %.4f\n",
                nm, i$est, i$p_js, i$p_holm))
  }
  if (!is.null(x$ci))
    cat(sprintf("  bootstrap: %d/%d full-pipeline resamples usable (percentile CIs stored)\n",
                x$ci$n_ok, x$ci$draws))
  invisible(x)
}
