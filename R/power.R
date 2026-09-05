## ---------------------------------------------------------------------------
## dmsa_power_*(): power and sample size for the path methods (POWER v2,
## 2026-09-05). Effect sizes are stated on the LATENT (true-score) scale and
## reliabilities per role, exactly as the engine treats them: a .70
## reliability attenuates the observed signal by sqrt(.70) per error-prone
## side; the corrected-moments estimator undoes the attenuation in the
## ESTIMATE but cannot recover the lost information, so power is that of
## the attenuated correlation. Two engines: closed forms (Fisher z on the
## attenuated partial correlation; instant) and simulation THROUGH the
## package's own estimator (.mm_core_v / .med_fl): mode "t" reads the
## engine's HC errors-in-variables t (the E15 pivot) against a normal
## reference (seconds), mode "perm" runs the shipped Freedman-Lane test
## (minutes). The E15 gate showed the two agree at the level and within a
## few points of power. Engine untouched.
## ---------------------------------------------------------------------------

#' Power and sample size for a single (partial) slope between error-prone scores
#'
#' Two-sided test of a latent (true-score) correlation \code{r} between two
#' scores of reliabilities \code{rel_x}, \code{rel_y}, with \code{k} other
#' predictors in the equation (they cost degrees of freedom); Fisher-z
#' approximation on the attenuated correlation \code{r * sqrt(rel_x * rel_y)}.
#'
#' @param r latent correlation (or semi-partial correlation) to detect.
#' @param n sample size.
#' @param rel_x,rel_y reliabilities (1 = observed without error).
#' @param k number of other predictors in the equation.
#' @param alpha two-sided level.
#' @param power target power for \code{dmsa_n_slope}.
#' @return power (a probability), or the n reaching \code{power}.
#' @export
dmsa_power_slope <- function(r, n, rel_x = 1, rel_y = 1, k = 0, alpha = 0.05) {
  r_obs <- r * sqrt(rel_x * rel_y)
  z <- atanh(r_obs) * sqrt(pmax(n - k - 3, 1))
  stats::pnorm(z - stats::qnorm(1 - alpha / 2)) + stats::pnorm(-z - stats::qnorm(1 - alpha / 2))
}
#' @rdname dmsa_power_slope
#' @export
dmsa_n_slope <- function(r, power = 0.80, rel_x = 1, rel_y = 1, k = 0, alpha = 0.05) {
  r_obs <- r * sqrt(rel_x * rel_y)
  ceiling(((stats::qnorm(1 - alpha / 2) + stats::qnorm(power)) / atanh(r_obs))^2 + k + 3)
}

## one simulated cell of the moderated-mediation model through the package core
.power_modmed_cell <- function(n, a1, a3, b1, b3, a2, c1, c2, rxw, rel, w_on, quadratic, k_cov, families,
                               mode, B, alpha, moderator = TRUE) {
  ## family blocks: INDEPENDENT family components in X, M and Y (within-family
  ## correlation of every variable, ICC ~ .36 / .16 / .16, no confounding of
  ## the paths - the blocks change the exchangeable unit, not the estimand)
  if (families) { nf <- ceiling(n / 2); fu <- function() rep(stats::rnorm(nf), each = 2)[seq_len(n)]; ux <- fu(); um <- fu(); uy <- fu() } else { ux <- 0; um <- 0; uy <- 0 }
  X <- if (families) .6 * ux + .8 * stats::rnorm(n) else stats::rnorm(n); X <- X / stats::sd(X)
  W <- rxw * X + sqrt(1 - rxw^2) * stats::rnorm(n)
  M <- a1 * X + a2 * W + a3 * X * W; M <- M + stats::rnorm(n) * sqrt(max(1 - stats::var(M), .2)) + .4 * um
  Y <- b1 * M + c1 * X + c2 * W + b3 * M * W; Y <- Y + stats::rnorm(n) * sqrt(max(1 - stats::var(Y), .2)) + .4 * uy
  mk <- function(v, rho) {                     # the engine's composite role: f / sqrt(rel * var), V = (1 - rel) / rel
    if (rho >= 1) return(list(f = v - mean(v), V = 0))
    z <- v + stats::rnorm(n, 0, sqrt(stats::var(v) * (1 - rho) / rho)); f <- z - mean(z)
    list(f = f / sqrt(rho * stats::var(f)), V = (1 - rho) / rho) }
  meas <- list(X = mk(X, rel[["X"]]), M = mk(M, rel[["M"]]), W = mk(W, rel[["W"]]), Y = mk(Y, rel[["Y"]]))
  covM <- if (k_cov > 0) matrix(stats::rnorm(n * k_cov), n, k_cov, dimnames = list(NULL, paste0("c", seq_len(k_cov)))) else NULL
  fam <- if (families) factor(rep(seq_len(ceiling(n / 2)), each = 2)[seq_len(n)]) else NULL
  core <- .mm_core_v(meas, "X", "M", if (moderator) "W" else NULL, "Y", if (moderator) w_on else "a", FALSE, quadratic, covM, colnames(covM),
                     if (mode == "perm") B else 0L, fam, if (families) "exchange" else "within", c(-1, 0, 1), 0, pivot = "hc")
  pv <- function(k) { p <- core$paths[[k]]; if (is.null(p)) return(NA_real_); if (mode == "perm") p$p else 2 * stats::pnorm(-abs(p$t)) }
  p <- c(a1 = pv("a[X]"), a3 = pv("a[X:W]"), b1 = pv("b[M]"), b3 = pv("b[M:W]"))
  pu <- if (!moderator) NA_real_ else if (w_on == "both") min(stats::p.adjust(c(max(p["a3"], p["b1"]), max(p["b3"], p["a1"]), max(p["a3"], p["b3"])), "holm")) else
    if (w_on == "a") max(p["a3"], p["b1"]) else max(p["b3"], p["a1"])
  est <- core$estimates
  c(a1 = p[["a1"]] <= alpha, a3 = p[["a3"]] <= alpha, b1 = p[["b1"]] <= alpha, b3 = p[["b3"]] <= alpha,
    IE = max(p["a1"], p["b1"]) <= alpha, IMM = pu <= alpha,
    est_a1 = est[["a1"]], est_a3 = est[["a3"]], est_b1 = est[["b1"]], est_b3 = est[["b3"]], shrink = core$moments$shrink)
}

#' Power for moderated mediation (and plain mediation) through the package estimator
#'
#' Simulates the structural model \code{M = a1 X + a2 W + a3 XW + e},
#' \code{Y = b1 M + c1 X + c2 W + b3 MW + e} with unit-variance latents,
#' adds measurement error to every role by its reliability, and runs the
#' package's corrected-moments estimator (the core of
#' \code{\link{dmsa_modmed}}) on every replicate. Power per path, for the
#' indirect effect \code{IE} (joint significance of a1 and b1) and for the
#' index of moderated mediation \code{IMM} (\code{w_on = "a"}: a3 with b1;
#' \code{"b"}: a1 with b3; \code{"both"}: the Holm union over the three
#' products - the registered test). \code{mode = "t"} reads the engine's HC
#' errors-in-variables t against a normal reference (seconds per cell);
#' \code{mode = "perm"} runs the Freedman-Lane permutation test itself with
#' \code{B} permutations (minutes). Plain mediation: \code{a3 = b3 = 0},
#' \code{w_on = "a"} and read \code{IE}. Also returns the mean estimate of
#' every path (the corrected estimator is unbiased in the mean; the
#' bootstrap CI note of the gate reports applies) and the rate at which
#' the conditioning guard shrank the offsets.
#'
#' @param n sample size (rows).
#' @param a1,a3,b1,b3 latent-scale paths: X -> M, XW -> M, M -> Y, MW -> Y.
#' @param a2,c1,c2 nuisance paths W -> M, X -> Y, W -> Y.
#' @param rxw correlation of X and W.
#' @param rel named reliabilities of \code{X, M, W, Y} (1 = observed).
#' @param w_on where the moderator acts (\code{"a"}, \code{"b"}, \code{"both"});
#'   \code{moderator = FALSE} drops W from the model (plain mediation:
#'   \code{M ~ X}, \code{Y ~ M + X}, read \code{IE}).
#' @param quadratic run the estimator with the quadratic safeguard (its power price).
#' @param k_cov number of nuisance covariates in both equations.
#' @param families rows come in exchangeable pairs (family blocks; the
#'   permutation exchanges blocks, as in the registered runs).
#' @param reps replicates; \code{mode}, \code{B} the test engine; \code{alpha} the level; \code{seed}.
#' @return a \code{dmsa_power} object: \code{$power} (a1, a3, b1, b3, IE, IMM),
#'   \code{$estimates} (mean corrected estimates against the truth),
#'   \code{$design}.
#' @export
dmsa_power_modmed <- function(n, a1, a3 = 0, b1, b3 = 0, a2 = 0.2, c1 = 0.2, c2 = 0.1, rxw = 0.3,
                              rel = c(X = 0.7, M = 0.7, W = 1, Y = 1), w_on = c("both", "a", "b"),
                              quadratic = FALSE, k_cov = 0L, families = FALSE, reps = 500L,
                              mode = c("t", "perm"), B = 199L, alpha = 0.05, seed = 1L, moderator = TRUE) {
  w_on <- match.arg(w_on); mode <- match.arg(mode)
  if (!moderator) { a3 <- 0; b3 <- 0; a2 <- 0; c2 <- 0; w_on <- "a" }   # no W in the world or the model
  rel0 <- c(X = 0.7, M = 0.7, W = 1, Y = 1); rel0[names(rel)] <- rel; rel <- rel0
  if (!is.null(seed)) set.seed(seed)
  R <- t(vapply(seq_len(reps), function(r) tryCatch(
    .power_modmed_cell(n, a1, a3, b1, b3, a2, c1, c2, rxw, rel, w_on, quadratic, k_cov, families, mode, B, alpha, moderator),
    error = function(e) rep(NA_real_, 11)), numeric(11)))
  colnames(R) <- c("a1", "a3", "b1", "b3", "IE", "IMM", "est_a1", "est_a3", "est_b1", "est_b3", "shrink")
  ok <- !is.na(R[, "IE"]) & (!moderator | !is.na(R[, "IMM"]))
  pw <- colMeans(R[ok, 1:6, drop = FALSE], na.rm = TRUE)
  if (w_on == "b") pw[["a3"]] <- NA_real_; if (w_on == "a") pw[["b3"]] <- NA_real_
  if (!moderator) pw[c("a3", "b3", "IMM")] <- NA_real_
  est <- rbind(truth = c(a1 = a1, a3 = a3, b1 = b1, b3 = b3), mean_estimate = colMeans(R[ok, 7:10, drop = FALSE]))
  colnames(est) <- c("a1", "a3", "b1", "b3")
  structure(list(method = if (!moderator) "mediation" else "moderated mediation", power = pw, estimates = est,
                 guard_rate = mean(R[ok, "shrink"] < 1), reps_ok = sum(ok),
                 design = list(n = n, a1 = a1, a3 = a3, b1 = b1, b3 = b3, a2 = a2, c1 = c1, c2 = c2, rxw = rxw, rel = rel, w_on = w_on,
                               quadratic = quadratic, k_cov = k_cov, families = families, reps = reps, mode = mode, B = if (mode == "perm") B else NA, alpha = alpha, moderator = moderator)),
            class = "dmsa_power")
}

#' Sample size for a target power of the moderated-mediation index (or any estimand)
#'
#' Bisection on \code{n} over \code{\link{dmsa_power_modmed}} in the fast mode.
#' @inheritParams dmsa_power_modmed
#' @param target which estimand: \code{"IMM"}, \code{"IE"}, \code{"a1"}, \code{"a3"}, \code{"b1"}, \code{"b3"}.
#' @param power target power; \code{n_range} the search bracket; \code{reps} per evaluation.
#' @return the smallest n on a grid of 10 reaching \code{power} (with the power found).
#' @export
dmsa_n_modmed <- function(target = "IMM", power = 0.80, n_range = c(60, 2000), reps = 300L, ...) {
  f <- function(n) dmsa_power_modmed(n = n, reps = reps, mode = "t", ...)$power[[target]]
  lo <- n_range[1]; hi <- n_range[2]
  if (f(hi) < power) return(c(n = NA_real_, power = f(hi)))
  while (hi - lo > 10) { mid <- round((lo + hi) / 2 / 10) * 10; if (f(mid) >= power) hi <- mid else lo <- mid }
  c(n = hi, power = f(hi))
}

#' Power for the parent-to-child paired slope at one level of the hierarchy
#'
#' Closed form for \code{child_u ~ mother_u + father_u + k covariates} with
#' latent slopes \code{beta_mother}, \code{beta_father}, parents correlated
#' \code{r_mf}, every score at reliability \code{rel} (1 at the probe
#' level). With \code{units > 1} the level-local maxT family is bounded by
#' Bonferroni (\code{alpha / units}, the conservative side of maxT), giving
#' the power to NAME the unit within its level.
#' @param n triads; \code{units} the family size at the level; \code{k} covariates.
#' @param beta_mother,beta_father latent slopes; \code{r_mf} parents' correlation; \code{rel} reliability of the scores; \code{alpha}.
#' @return a data.frame: one row per parent with the raw and the level-named power.
#' @export
dmsa_power_transmit <- function(n, beta_mother, beta_father = 0, r_mf = 0.2, rel = 0.7, units = 1L, k = 4L, alpha = 0.05) {
  vy <- beta_mother^2 + beta_father^2 + 2 * beta_mother * beta_father * r_mf
  vy <- vy + (1 - min(.95, vy))
  sp <- c(mother = beta_mother, father = beta_father) * sqrt(1 - r_mf^2) / sqrt(vy)
  data.frame(parent = names(sp), semi_partial_r = as.numeric(sp),
             power_raw = dmsa_power_slope(as.numeric(sp), n, rel, rel, k + 1, alpha),
             power_named = dmsa_power_slope(as.numeric(sp), n, rel, rel, k + 1, alpha / max(units, 1)),
             units = units, stringsAsFactors = FALSE)
}

#' Power for dyadic concordance at one level, its outcome surface and its exposure form
#'
#' \code{dmsa_power_concordance}: within-couple correlation \code{r_couple}
#' of scores at reliability \code{rel} against the re-pairing null, per
#' unit, and the power to name a unit within a level of \code{units}
#' (Bonferroni bound of maxT). \code{dmsa_power_exposure}: the slope change
#' gamma3 read as the difference of two within-couple correlations
#' (\code{r0} unexposed, \code{r1} exposed, proportion \code{p1}), Fisher z.
#' \code{dmsa_power_surface}: simulation of the corrected dyadic surface
#' (curvature \code{a4} along the incongruence line, partner scores at
#' reliability \code{rel}, couple correlation \code{r_couple}) through the
#' MATH 1.6 offsets with the conditioning guard, the HC contrast statistic
#' against a normal reference; also the guard rate and the median ratio of
#' the corrected curvature to the truth.
#' @param n_couples couples; \code{r_couple} latent within-couple correlation; \code{rel} reliability; \code{units} family size; \code{alpha}.
#' @return \code{dmsa_power_concordance} / \code{dmsa_power_exposure}: named numeric (power, and the n for .80);
#'   \code{dmsa_power_surface}: named numeric (power, guard_rate, median_a4_ratio).
#' @export
dmsa_power_concordance <- function(n_couples, r_couple, rel = 0.7, units = 1L, alpha = 0.05) {
  c(power_raw = dmsa_power_slope(r_couple, n_couples, rel, rel, 0, alpha),
    power_named = dmsa_power_slope(r_couple, n_couples, rel, rel, 0, alpha / max(units, 1)),
    n_for_80 = dmsa_n_slope(r_couple, 0.80, rel, rel, 0, alpha))
}
#' @rdname dmsa_power_concordance
#' @param r0,r1 within-couple correlations without / with the exposure; \code{p1} proportion exposed.
#' @export
dmsa_power_exposure <- function(n_couples, r0, r1, p1 = 0.5, rel = 0.7, alpha = 0.05) {
  pw <- function(nc) { z0 <- atanh(r0 * rel); z1 <- atanh(r1 * rel); n0 <- nc * (1 - p1); n1 <- nc * p1
    se <- sqrt(1 / (n0 - 3) + 1 / (n1 - 3)); z <- (z1 - z0) / se
    stats::pnorm(z - stats::qnorm(1 - alpha / 2)) + stats::pnorm(-z - stats::qnorm(1 - alpha / 2)) }
  nn <- 20; while (pw(nn) < 0.80 && nn < 1e5) nn <- nn + 2
  c(power = pw(n_couples), n_for_80 = nn)
}
#' @rdname dmsa_power_concordance
#' @param n actor rows of the surface; \code{a4} latent curvature along the incongruence line; \code{b1,b2} actor and partner linear effects; \code{reps}, \code{seed};
#'   \code{mode = "perm"} (default) runs the Freedman-Lane contrast test with \code{B} permutations - the HC contrast statistic against a normal reference (\code{"t"}) is conservative for the corrected surface.
#' @export
dmsa_power_surface <- function(n, a4, r_couple = 0.1, rel = 0.7, b1 = 0.2, b2 = 0.1, reps = 500L, alpha = 0.05, seed = 2L,
                               mode = c("perm", "t"), B = 199L) {
  mode <- match.arg(mode)
  set.seed(seed); b3 <- a4 / 4; b5 <- a4 / 4; b4 <- -a4 / 2      # a4 = b3 - b4 + b5
  rej <- rep(NA, reps); ratio <- rep(NA_real_, reps); guard <- rep(NA, reps)
  xn <- c("act", "par", "act2", "actpar", "par2"); cv <- stats::setNames(c(0, 0, 1, -1, 1), xn)
  for (r in seq_len(reps)) {
    ea <- stats::rnorm(n); eb <- r_couple * ea + sqrt(1 - r_couple^2) * stats::rnorm(n)
    mu <- b1 * ea + b2 * eb + b3 * ea^2 + b4 * ea * eb + b5 * eb^2
    y <- mu + stats::rnorm(n) * sqrt(max(1 - stats::var(mu), .3))
    V <- (1 - rel) / rel; fa <- ea + stats::rnorm(n, 0, sqrt(V)); fb <- eb + stats::rnorm(n, 0, sqrt(V))
    D <- scale(cbind(act = fa, par = fb, act2 = fa^2, actpar = fa * fb, par2 = fb^2, y = y), scale = FALSE)
    mm <- tryCatch(.mm_moments_sq(D, voffset = c(act = V, par = V), products = list(act2 = c("act", "act"), actpar = c("act", "par"), par2 = c("par", "par")), n = n), error = function(e) NULL)
    if (is.null(mm)) next
    perms <- if (mode == "perm") .hier_perms(NULL, n, B) else list()
    fl <- tryCatch(.hier_fl(D, mm$S, mm$Sobs, xn, "y", NULL, n, perms, contrast = cv), error = function(e) NULL)
    if (is.null(fl)) next
    rej[r] <- if (mode == "perm") (1 + sum(abs(fl$tn) >= abs(fl$t))) / (B + 1) <= alpha else abs(fl$t) > stats::qnorm(1 - alpha / 2)
    ratio[r] <- fl$b / a4; guard[r] <- mm$shrink < 1
  }
  c(power = mean(rej, na.rm = TRUE), guard_rate = mean(guard, na.rm = TRUE), median_a4_ratio = stats::median(ratio, na.rm = TRUE), reps_ok = sum(!is.na(rej)))
}

#' @export
print.dmsa_power <- function(x, ...) {
  d <- x$design
  cat(sprintf("Power for %s through the package estimator (%s mode%s; %d reps, n = %d, alpha %.2f)\n", x$method, d$mode,
              if (d$mode == "perm") sprintf(", B = %d", d$B) else "", x$reps_ok, d$n, d$alpha))
  cat(sprintf("  paths (latent scale): a1 %.2f  a3 %.2f  b1 %.2f  b3 %.2f   w_on %s%s   rel X %.2f M %.2f W %.2f Y %.2f%s%s\n",
              d$a1, d$a3, d$b1, d$b3, d$w_on, if (d$quadratic) " + quadratic safeguard" else "", d$rel[["X"]], d$rel[["M"]], d$rel[["W"]], d$rel[["Y"]],
              if (d$k_cov > 0) sprintf("   %d covariates", d$k_cov) else "", if (d$families) "   family blocks" else ""))
  cat("  power: ", paste(sprintf("%s %.2f", names(x$power), x$power), collapse = "   "), "\n", sep = "")
  cat("  mean estimate / truth: ", paste(sprintf("%s %.3f/%.2f", colnames(x$estimates), x$estimates["mean_estimate", ], x$estimates["truth", ]), collapse = "   "),
      sprintf("   guard fired in %.0f%% of reps\n", 100 * x$guard_rate), sep = "")
  invisible(x)
}
