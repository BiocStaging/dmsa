# latent moderated mediation: measurement ladder, corrected moments,
# permutation inference, and the end-to-end engine

gen_latent <- function(n, a1 = 0.3, a3 = 0.35, b1 = 0.5, lam = 0.7,
                       p = c(hpa = 6, oxytocin = 8, kyn = 5)) {
  hpa <- rnorm(n); ot <- rnorm(n)
  kp <- a1 * hpa + a3 * hpa * ot +
    rnorm(n, sd = sqrt(max(0.05, 1 - a1^2 - a3^2)))
  y <- b1 * kp + rnorm(n)
  mk <- function(f, pp) sapply(seq_len(pp), function(j)
    lam * f + rnorm(n, sd = sqrt(1 - lam^2)))
  list(tones = list(hpa = mk(hpa, p[1]), oxytocin = mk(ot, p[2]),
                    kyn = mk(kp, p[3])),
       dat = data.frame(y = y, age = rnorm(n)))
}

test_that("one-factor ML recovers loadings, tau and congeneric", {
  set.seed(1)
  n <- 5000; f <- rnorm(n)
  lam0 <- c(.5, .6, .7, .8)
  Z <- sapply(lam0, function(l) l * f + rnorm(n, sd = sqrt(1 - l^2)))
  S <- cov(Z)
  fc <- dmsa:::.med_fa1(S, n, model = "congeneric")
  expect_true(fc$conv); expect_false(fc$heywood)
  expect_equal(fc$lambda, lam0, tolerance = 0.06)
  ft <- dmsa:::.med_fa1(S, n, model = "tau")
  expect_true(ft$conv)
  expect_equal(unique(round(ft$lambda, 10)), round(ft$lambda[1], 10))
  expect_equal(mean(ft$lambda), mean(lam0), tolerance = 0.08)
})

test_that("Bartlett mapping: sane V, omega, determinacy", {
  set.seed(2)
  n <- 400; f <- rnorm(n); lam <- rep(.7, 6)
  Z <- sapply(lam, function(l) l * f + rnorm(n, sd = sqrt(1 - l^2)))
  bl <- dmsa:::.med_bartlett(scale(Z, scale = FALSE), lam, 1 - lam^2)
  expect_gt(bl$V, 0); expect_lt(bl$V, 1)
  expect_gt(bl$omega, .8); expect_lt(bl$omega, 1)
  expect_gt(bl$determinacy, .9)
  expect_gt(cor(bl$f, f), .85)
})

test_that("measurement ladder: congeneric steps down, <3 indicators -> single", {
  set.seed(3)
  g <- gen_latent(120)
  ms <- dmsa:::.med_measure(g$tones$hpa, measurement = "tau", label = "hpa")
  expect_identical(ms$model, "tau")
  two <- g$tones$kyn[, 1:2]
  ms2 <- dmsa:::.med_measure(two, measurement = "tau", label = "kyn")
  expect_identical(ms2$model, "single")
  expect_true(length(ms2$notes) >= 1L)
  expect_gt(ms2$V, 0)
})

test_that("corrected moments subtract the right offsets and guard PD", {
  set.seed(4)
  n <- 300
  D <- cbind(a = rnorm(n), b = rnorm(n))
  D <- cbind(D, "a:b" = D[, "a"] * D[, "b"])
  D <- scale(D, scale = FALSE)
  mm <- dmsa:::.med_moments(D, voffset = c(a = .2, b = .3),
                            products = list("a:b" = c("a", "b")), n = n)
  expect_equal(mm$Sobs["a", "a"] - mm$S["a", "a"], .2, tolerance = 1e-10)
  expect_gt(mm$Sobs["a:b", "a:b"] - mm$S["a:b", "a:b"], 0)
  expect_false(mm$lambda_triggered)
  ## absurd offsets must trigger the lambda shrink, not a crash
  mm2 <- dmsa:::.med_moments(D, voffset = c(a = 5, b = 5),
                             products = list("a:b" = c("a", "b")), n = n)
  expect_true(mm2$lambda_triggered)
  ev <- eigen(mm2$S, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev), 0)
})

test_that("core recovers structural effects the naive composite attenuates", {
  set.seed(5)
  R <- 12; corr <- naiv <- matrix(NA, R, 2)
  for (r in seq_len(R)) {
    g <- gen_latent(250, lam = 0.6)
    cc <- dmsa:::.med_core(g$tones, g$dat, x = "hpa", m = "kyn",
                           w = "oxytocin", outcome = "y", covariates = "age",
                           measurement = "tau", B = 0L, blocks = NULL,
                           model8 = FALSE, conditional_at = 0, alpha_ssc = 0)
    cx <- rowMeans(g$tones$hpa); cw <- rowMeans(g$tones$oxytocin)
    cm <- rowMeans(g$tones$kyn)
    corr[r, ] <- c(cc$paths[["a[hpa:oxytocin]"]]$b, cc$paths[["b[kyn]"]]$b)
    naiv[r, ] <- c(coef(lm(cm ~ cx * cw + g$dat$age))[["cx:cw"]],
                   coef(lm(g$dat$y ~ cm + cx + cw + g$dat$age))[["cm"]])
  }
  ## b1: corrected mean within 15% of truth; naive composite off by more
  expect_lt(abs(mean(corr[, 2]) / 0.5 - 1), 0.15)
  expect_gt(abs(mean(naiv[, 2]) / 0.5 - 1),
            abs(mean(corr[, 2]) / 0.5 - 1))
})

test_that("joint significance stays conservative at the composite null", {
  set.seed(6)
  R <- 25; pj <- numeric(R)
  for (r in seq_len(R)) {
    g <- gen_latent(120, a3 = 0)          # corner: a3 = 0, b1 != 0
    cc <- dmsa:::.med_core(g$tones, g$dat, x = "hpa", m = "kyn",
                           w = "oxytocin", outcome = "y", covariates = "age",
                           measurement = "tau", B = 99L, blocks = NULL,
                           model8 = FALSE, conditional_at = 0, alpha_ssc = 0)
    pj[r] <- cc$indices[["IMM[hpa]"]]$p_js
  }
  expect_lte(mean(pj <= 0.05), 0.16)      # loose bound for 25 reps
})

test_that("dmsa_mediate runs end-to-end on the bundled map", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv",
                                 package = "dmsa")))
  set.seed(7)
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction) & mp$system_id %in% c(1, 2, 15, 24), ]
  pr <- unique(mp$probe)
  n <- 96
  B <- matrix(plogis(rnorm(n * length(pr), sd = 1.2)), n, length(pr),
              dimnames = list(NULL, pr))
  dat <- data.frame(EPDS = rnorm(n), age = rnorm(n),
                    chip = factor(rep(1:12, each = 8)))
  dat$EPDS[3] <- NA                       # complete-case path
  fit <- suppressWarnings(dmsa_mediate(
    dat, methylation = B, map = "alpha",
    x = c("hpa", "immune"), m = "kynurenine", w = "oxytocin",
    outcome = "EPDS", covariates = "age", chip = "chip",
    measurement = "tau", B = 99, boot = 0, seed = 1))
  expect_s3_class(fit, "dmsa_mediate")
  expect_equal(fit$n, n - 1)              # the NA row dropped
  expect_named(fit$indices, c("IMM[hpa]", "IMM[immune]"))
  ## noise betas: nothing may reach significance
  expect_true(all(vapply(fit$indices, function(i) i$p_holm, 1) > .05))
  expect_true(all(fit$coverage$probes[fit$coverage$used] >= 3))
  expect_output(print(fit), "Latent DMSA moderated mediation")
  ## block-restricted permutation and model 8 both run
  fit2 <- suppressWarnings(dmsa_mediate(
    dat, methylation = B, map = "alpha",
    x = "hpa", m = "kynurenine", w = "oxytocin",
    outcome = "EPDS", covariates = "age", chip = "chip",
    perm_blocks = "chip", model = "8", B = 49, seed = 2))
  expect_true(any(grepl("^c\\[", names(fit2$paths))))
  expect_true(is.finite(fit2$paths[["a[hpa:oxytocin]"]]$p))
})

test_that("plain mediation (w = NULL) and RNG restoration", {
  set.seed(8)
  g <- gen_latent(100)
  cc <- dmsa:::.med_core(g$tones, g$dat, x = "hpa", m = "kyn", w = NULL,
                         outcome = "y", covariates = "age",
                         measurement = "tau", B = 49L, blocks = NULL,
                         model8 = FALSE, conditional_at = 0, alpha_ssc = 0)
  expect_named(cc$indices, "IE[hpa]")
  expect_true(is.finite(cc$indices[["IE[hpa]"]]$p_js))
  ## seed restoration contract (as in every other engine)
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv",
                                 package = "dmsa")))
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction) & mp$system_id %in% c(1, 2, 15), ]
  pr <- unique(mp$probe); n <- 60
  B <- matrix(plogis(rnorm(n * length(pr))), n, length(pr),
              dimnames = list(NULL, pr))
  dat <- data.frame(y = rnorm(n), age = rnorm(n))
  set.seed(99); before <- .Random.seed
  invisible(suppressWarnings(dmsa_mediate(
    dat, methylation = B, map = "alpha", x = "hpa", m = "kynurenine",
    w = "oxytocin", outcome = "y", covariates = "age",
    B = 19, seed = 42)))
  expect_identical(.Random.seed, before)
})

test_that("single_rel fixes the fallback reliability a priori", {
  set.seed(9)
  g <- gen_latent(100)
  two <- g$tones$kyn[, 1:2]
  ms <- dmsa:::.med_measure(two, measurement = "tau", label = "kyn",
                            single_rel = 0.70)
  expect_identical(ms$model, "single")
  expect_equal(ms$omega, 0.70, tolerance = 1e-10)
  ## unit-latent-variance conditioning: V = (1 - rel) / rel exactly
  expect_equal(ms$V, 0.30 / 0.70, tolerance = 1e-8)
})

test_that("family-block exchange permutation and cluster bootstrap run", {
  set.seed(10)
  g <- gen_latent(120)
  fam <- factor(rep(1:60, each = 2))
  cc <- dmsa:::.med_core(g$tones, g$dat, x = "hpa", m = "kyn",
                         w = "oxytocin", outcome = "y", covariates = "age",
                         measurement = "tau", B = 99L, blocks = fam,
                         model8 = FALSE, conditional_at = 0, alpha_ssc = 0,
                         block_mode = "exchange")
  expect_true(is.finite(cc$paths[["a[hpa:oxytocin]"]]$p))
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv",
                                 package = "dmsa")))
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction) & mp$system_id %in% c(1, 2, 15), ]
  pr <- unique(mp$probe); n <- 80
  B <- matrix(plogis(rnorm(n * length(pr))), n, length(pr),
              dimnames = list(NULL, pr))
  dat <- data.frame(y = rnorm(n), age = rnorm(n),
                    cID = rep(1:40, each = 2),
                    chip = factor(rep(1:10, each = 8)))
  fit <- suppressWarnings(dmsa_mediate(
    dat, methylation = B, map = "alpha", x = "hpa", m = "kynurenine",
    w = "oxytocin", outcome = "y", covariates = "age", chip = "chip",
    block = "cID", B = 49, boot = 20, seed = 3))
  expect_identical(fit$perm_blocks, "family-exchange")
  expect_true(fit$ci$n_ok > 0)
  ## block + chip permutation must refuse
  expect_error(dmsa_mediate(
    dat, methylation = B, map = "alpha", x = "hpa", m = "kynurenine",
    w = "oxytocin", outcome = "y", covariates = "age", chip = "chip",
    block = "cID", perm_blocks = "chip", B = 19),
    "cannot be combined")
})

## E14 (2026-09-02, PI ruling): the Freedman-Lane pivot is the errors-in-
## variables sandwich, not the plain-OLS rss on a corrected coefficient.
test_that("E14: the FL pivot is the EIV sandwich and reduces to OLS for observed regressors", {
  set.seed(14)
  n <- 150
  X <- rnorm(n); W <- .3 * X + sqrt(.91) * rnorm(n); Y <- .4 * X + .2 * W + rnorm(n)
  D <- scale(cbind(X = X, W = W, y = Y), scale = FALSE)
  ## (i) no measurement error: the pivot must equal the classical OLS t exactly
  mm <- dmsa:::.med_moments(D, voffset = numeric(0), products = list(), n = n)
  fl <- dmsa:::.med_fl(D, mm$S, c("X", "W"), "y", "X", n, B = 0L, Sobs = mm$Sobs)
  ols <- summary(lm(Y ~ X + W))$coefficients
  expect_equal(unname(fl$b), unname(ols["X", "Estimate"]), tolerance = 1e-8)
  expect_equal(unname(fl$se), unname(ols["X", "Std. Error"]), tolerance = 1e-8)
  expect_equal(unname(fl$t), unname(ols["X", "t value"]), tolerance = 1e-8)
  ## (ii) with error on X (reliability .7): se equals the sandwich formula
  ##      sigma2_obs * [S*xx^-1 S_obs S*xx^-1]_jj / (n - k), sigma2_obs from the
  ##      UNCORRECTED response variance
  rel <- .7
  Xe <- X + rnorm(n, 0, sqrt(var(X) * (1 - rel) / rel)); Xe <- (Xe - mean(Xe)) / sqrt(rel * var(Xe))
  V <- (1 - rel) / rel
  D2 <- scale(cbind(X = Xe, W = W, y = Y), scale = FALSE)
  mm2 <- dmsa:::.med_moments(D2, voffset = c(X = V), products = list(), n = n)
  fl2 <- dmsa:::.med_fl(D2, mm2$S, c("X", "W"), "y", "X", n, B = 0L, Sobs = mm2$Sobs)
  xn <- c("X", "W"); Sxxi <- solve(mm2$S[xn, xn]); So <- mm2$Sobs[xn, xn]
  b <- as.numeric(Sxxi %*% mm2$S[xn, "y"])
  s2 <- mm2$Sobs["y", "y"] - 2 * sum(b * mm2$S[xn, "y"]) + as.numeric(t(b) %*% So %*% b)
  se_ref <- sqrt(s2 * (Sxxi %*% So %*% Sxxi)[1, 1] / (n - 3))
  expect_equal(unname(fl2$b), b[1], tolerance = 1e-10)
  expect_equal(unname(fl2$se), unname(se_ref), tolerance = 1e-10)
  expect_gt(s2, 0)
  ## the corrected coefficient is disattenuated relative to the raw one
  expect_gt(abs(fl2$b), abs(coef(lm(Y ~ Xe + W))[["Xe"]]) * 1.05)
})

test_that("E14: no degenerate permutations with error-prone regressors, level plausible", {
  set.seed(15)
  n <- 160; rel <- .7; B <- 199L
  pv <- numeric(20); deg <- 0L
  for (r in seq_len(20)) {
    X <- rnorm(n); W <- .3 * X + sqrt(.91) * rnorm(n); M <- .4 * X + .2 * W + rnorm(n) * .8
    noisy <- function(v) { e <- v + rnorm(n, 0, sqrt(var(v) * (1 - rel) / rel)); (e - mean(e)) / sqrt(rel * var(e)) }
    F <- cbind(X = noisy(X), W = noisy(W))
    D <- scale(cbind(F, `X:W` = F[, "X"] * F[, "W"], m = M - mean(M)), scale = FALSE)
    V <- (1 - rel) / rel
    mm <- dmsa:::.med_moments(D, voffset = c(X = V, W = V), products = list(`X:W` = c("X", "W")), n = n)
    fl <- dmsa:::.med_fl(D, mm$S, c("X", "W", "X:W"), "m", "X:W", n, B, Sobs = mm$Sobs)   # true a3 = 0
    pv[r] <- fl$p; deg <- deg + fl$degenerate
  }
  expect_identical(deg, 0L)
  expect_lte(mean(pv <= .05), 0.25)     # loose bound, 20 reps at a true null
  expect_true(all(pv > 0 & pv <= 1))
})
