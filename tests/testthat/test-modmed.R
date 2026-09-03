# dmsa_modmed(): the general O/L moderated-mediation front end. It must
# (i) reproduce dmsa_mediate() model 7 exactly when every role is a system,
# (ii) run every role mix (observed / composite / latent, latent outcome),
# (iii) reduce its square-capable moment correction to .med_moments() when
# no square is present, and (iv) bootstrap with the rung held fixed.

mk_data <- function(n = 120, seed = 11) {
  set.seed(seed)
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv", package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction) & mp$system_id %in% c(1, 2, 15, 24), ]
  pr <- unique(mp$probe)
  B <- matrix(plogis(rnorm(n * length(pr), sd = 1.2)), n, length(pr), dimnames = list(NULL, pr))
  dat <- data.frame(EPDS = rnorm(n), age = rnorm(n), support = rnorm(n),
                    chip = factor(rep(seq_len(n / 8), each = 8)), fam = factor(rep(seq_len(n / 2), each = 2)))
  list(B = B, dat = dat)
}

test_that("dmsa_modmed reproduces dmsa_mediate model 7 exactly", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_data()
  f1 <- suppressWarnings(dmsa_mediate(d$dat, methylation = d$B, x = "hpa", m = "kynurenine", w = "oxytocin",
    outcome = "EPDS", covariates = "age", chip = "chip", block = "fam", measurement = "congeneric", B = 99, seed = 1))
  f2 <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, x = "hpa", m = "kynurenine", w = "oxytocin", y = "EPDS",
    w_on = "a", covariates = "age", chip = "chip", block = "fam", measurement = "congeneric", B = 99, seed = 1))
  expect_s3_class(f2, "dmsa_modmed")
  expect_setequal(names(f1$paths), names(f2$paths))
  for (k in names(f1$paths)) {
    expect_equal(f2$paths[[k]]$b, f1$paths[[k]]$b, tolerance = 1e-12, info = k)
    expect_identical(f2$paths[[k]]$p, f1$paths[[k]]$p, info = k)
  }
  expect_equal(f2$indices$IMM$est, f1$indices[["IMM[hpa]"]]$est, tolerance = 1e-12)
  expect_identical(f2$indices$IMM$p_js, f1$indices[["IMM[hpa]"]]$p_js)
  expect_equal(f2$moments$shrink, f1$moments$shrink)
  expect_output(print(f2), "first stage")
})

test_that("role resolution: systems vs columns, composite reliability, refusals", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_data()
  expect_error(dmsa_modmed(d$dat, methylation = d$B, x = "hpa", m = "nope", y = "EPDS", B = 9), "neither a system")
  expect_error(dmsa_modmed(d$dat, methylation = d$B, x = "hpa", m = "kynurenine", y = "EPDS", roles = c(y = "latent"), B = 9), "cannot be 'latent'")
  expect_error(dmsa_modmed(d$dat, methylation = d$B, x = "hpa", m = "kynurenine", y = "EPDS", roles = c(x = "observed"), B = 9), "cannot be 'observed'")
  expect_error(dmsa_modmed(d$dat, methylation = d$B, x = "hpa", m = "hpa", y = "EPDS", B = 9), "distinct")
  f <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, x = "support", m = "hpa", w = "oxytocin", y = "EPDS",
    roles = c(m = "composite", w = "composite"), composite_rel = c(m = .70, w = .60), w_on = "b", covariates = "age", chip = "chip", B = 19, seed = 4))
  expect_match(f$roles[["support"]], "^O$")
  expect_match(f$roles[["hpa"]], "rel 0.70"); expect_match(f$roles[["oxytocin"]], "rel 0.60")
  expect_equal(f$measurement$hpa$V, (1 - .7) / .7)
  expect_true("b[hpa:oxytocin]" %in% names(f$paths)); expect_false("a[support:oxytocin]" %in% names(f$paths))
  expect_match(f$indices$IMM$kind, "second stage")
})

test_that("dual stage with a latent outcome and the quadratic safeguard", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_data()
  f <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, x = "hpa", m = "kynurenine", w = "support", y = "immune",
    w_on = "both", quadratic = TRUE, direct_moderation = TRUE, covariates = "age", chip = "chip", block = "fam", B = 19, seed = 5))
  expect_match(f$roles[["immune"]], "^L")
  expect_true(all(c("a[hpa:support]", "b[kynurenine:support]", "c'[hpa:support]", "a[hpa^2]", "a[support^2]") %in% names(f$paths)))
  expect_true(all(c("kynurenine^2", "hpa^2", "support^2") %in% f$terms$Y_model))
  expect_named(f$indices$IMM$products, c("I_a", "I_b", "I_ab"))
  expect_equal(f$indices$IMM$p_js, min(f$indices$IMM$p_holm))
  expect_equal(f$indices$IMM$est, f$conditional$indirect[3] - f$conditional$indirect[1])
  expect_output(print(f), "dual stage")
})

test_that("square-capable moments equal .med_moments without squares and add the MATH 1.6 offsets with them", {
  set.seed(6); n <- 300
  D <- scale(cbind(x = rnorm(n), w = rnorm(n), y = rnorm(n)), scale = FALSE)
  D <- cbind(D, `x:w` = D[, "x"] * D[, "w"]); D <- scale(D, scale = FALSE)
  vo <- c(x = .3, w = .2); pr <- list(`x:w` = c("x", "w"))
  m1 <- dmsa:::.med_moments(D, vo, pr, n); m2 <- dmsa:::.mm_moments_sq(D, vo, pr, n)
  expect_equal(m1$S, m2$S); expect_equal(m1$shrink, m2$shrink)
  D2 <- cbind(D, `x^2` = D[, "x"]^2); D2 <- scale(D2, scale = FALSE)
  m3 <- dmsa:::.mm_moments_sq(D2, vo, c(pr, list(`x^2` = c("x", "x"))), n)
  sxx <- m3$Sobs["x", "x"] - .3
  expect_equal(m3$O["x^2", "x^2"], 4 * sxx * .3 + 2 * .3^2)
  expect_equal(m3$O["x^2", "x:w"], 2 * m3$Sobs["x", "w"] * .3)
})

test_that("plain mediation and the rung-fixed bootstrap", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_data()
  f <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, x = "support", m = "hpa", w = NULL, y = "EPDS",
    covariates = "age", chip = "chip", block = "fam", B = 19, boot = 12, seed = 7))
  expect_null(f$w); expect_named(f$indices, "IE")
  expect_true(all(c("a[support]", "b[hpa]", "c'[support]") %in% names(f$paths)))
  expect_equal(f$ci$draws, 12); expect_true(f$ci$n_ok >= 10)
  expect_true(all(c("percentile", "basic") %in% names(f$ci)))
  expect_true(all(f$ci$percentile[, "lower"] <= f$ci$percentile[, "upper"]))
  expect_equal(unname(f$ci$estimate["IE"]), f$indices$IE$est)
  expect_equal(f$ci$rung, "fixed")
  expect_output(print(f), "plain mediation")
  ## RNG restoration
  set.seed(99); before <- .Random.seed
  invisible(suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, x = "support", m = "hpa", y = "EPDS", covariates = "age", B = 9, seed = 42)))
  expect_identical(.Random.seed, before)
  ## E16: an empty reduced model (one predictor, no covariates) now runs
  f16 <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, x = "support", m = "hpa", y = "EPDS", B = 19, seed = 16))
  expect_true(is.finite(f16$paths[["a[support]"]]$p))
})
