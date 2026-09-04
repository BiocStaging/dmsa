# dmsa_modmed() 0.99.14: vector roles (several X, several parallel M), the
# measurement rule v1.2 (sex + global mean), the direction-source default,
# the connected measurement of two latent mediators, the sensitivity
# battery, and the hand-over from dmsa_mediate() with several mediators.

mk_data_v <- function(n = 120, seed = 11) {
  set.seed(seed)
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv", package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction) & mp$system_id %in% c(1, 2, 15, 24), ]
  pr <- unique(mp$probe)
  B <- matrix(plogis(rnorm(n * length(pr), sd = 1.2)), n, length(pr), dimnames = list(NULL, pr))
  dat <- data.frame(EPDS = rnorm(n), age = rnorm(n), support = rnorm(n), stress = rnorm(n), sex = rep(1:2, n / 2),
                    chip = factor(rep(seq_len(n / 8), each = 8)), fam = factor(rep(seq_len(n / 2), each = 2)))
  list(B = B, dat = dat)
}

test_that("the vector core reproduces the scalar core exactly", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_data_v()
  ## same measurement, same permutations: the two cores must agree to the digit
  st <- dmsa:::.paths_setup(d$dat, d$B, map = "alpha", systems = c("hpa", "kynurenine", "oxytocin"), columns = "EPDS",
                            covariates = "age", chip = "chip", block = "fam", global_mean = FALSE)
  spec <- function(v) list(name = v, kind = if (v == "EPDS") "column" else "system", role = if (v == "EPDS") "observed" else "latent")
  vals <- c("hpa", "kynurenine", "oxytocin", "EPDS"); crel <- stats::setNames(rep(.7, 4), vals)
  meas <- lapply(vals, function(v) dmsa:::.mm_measure_role(v, spec(v), st$tones, st$dat2, "congeneric", NULL, crel, st$n)); names(meas) <- vals
  set.seed(5); c1 <- dmsa:::.mm_core(meas, "hpa", "kynurenine", "oxytocin", "EPDS", "both", TRUE, TRUE, st$covM, "age", 49, st$blocks, "exchange", c(-1, 0, 1), 0)
  set.seed(5); c2 <- dmsa:::.mm_core_v(meas, "hpa", "kynurenine", "oxytocin", "EPDS", "both", TRUE, TRUE, st$covM, "age", 49, st$blocks, "exchange", c(-1, 0, 1), 0)
  expect_identical(names(c1$paths), names(c2$paths))
  for (k in names(c1$paths)) { expect_equal(c2$paths[[k]]$b, c1$paths[[k]]$b, tolerance = 1e-12, info = k); expect_identical(c2$paths[[k]]$p, c1$paths[[k]]$p, info = k) }
  expect_equal(c2$indices$IMM$est, c1$indices$IMM$est, tolerance = 1e-12)
  expect_identical(c2$indices$IMM$p_holm, c1$indices$IMM$p_holm)
  expect_equal(c2$estimates[c("a1", "a3", "b1", "b3", "IE", "dtheta")], c1$estimates[c("a1", "a3", "b1", "b3", "IE", "dtheta")], tolerance = 1e-12)
  expect_equal(c2$conditional, c1$conditional, tolerance = 1e-12)
})

test_that("several X and several parallel M: one M-model per mediator, products per pair, Holm over the family", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_data_v()
  f <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = c("support", "stress"), m = c("hpa", "kynurenine"), w = "oxytocin", y = "EPDS",
                                    w_on = "a", covariates = "age", chip = "chip", block = "fam", B = 19, seed = 3, global_mean = FALSE, theta_mm = "none"))
  expect_s3_class(f, "dmsa_modmed")
  expect_true(all(c("a[support->hpa]", "a[stress->hpa]", "a[support->kynurenine]", "a[support:oxytocin->hpa]", "b[hpa]", "b[kynurenine]", "c'[support]", "c'[stress]") %in% names(f$paths)))
  expect_setequal(names(f$indices), c(paste0("IE[", c("support->hpa", "support->kynurenine", "stress->hpa", "stress->kynurenine"), "]"),
                                      paste0("IMM[", c("support->hpa", "support->kynurenine", "stress->hpa", "stress->kynurenine"), "]")))
  ie <- f$indices[grep("^IE", names(f$indices))]
  expect_equal(vapply(ie, function(i) i$p_holm_family, 1), stats::p.adjust(vapply(ie, function(i) i$p_js, 1), "holm"))
  expect_equal(f$indices[["IE[support->hpa]"]]$est, f$paths[["a[support->hpa]"]]$b * f$paths[["b[hpa]"]]$b)
  expect_equal(nrow(f$conditional), 4 * 3); expect_true(all(c("x", "m") %in% names(f$conditional)))
  expect_output(print(f), "support -> hpa")
  ## dmsa_mediate() hands several mediators over
  g <- suppressWarnings(dmsa_mediate(d$dat, methylation = d$B, map = "alpha", x = "support", m = c("hpa", "kynurenine"), w = "oxytocin", outcome = "EPDS",
                                     covariates = "age", chip = "chip", block = "fam", measurement = "congeneric", B = 9, seed = 3, global_mean = FALSE))
  expect_s3_class(g, "dmsa_modmed"); expect_equal(g$w_on, "a"); expect_equal(g$m, c("hpa", "kynurenine"))
  ## plain parallel mediation (no W), with the bootstrap
  h <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = "support", m = c("hpa", "kynurenine"), y = "EPDS",
                                    covariates = "age", chip = "chip", block = "fam", B = 9, boot = 5, seed = 3, global_mean = FALSE, theta_mm = "none"))
  expect_true(all(grepl("^IE", names(h$indices)))); expect_equal(h$ci$draws, 5)
  expect_true(all(c("IE[support->hpa]", "IE[support->kynurenine]") %in% rownames(h$ci$percentile)))
})

test_that("measurement rule v1.2 is the default and is switchable", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_data_v()
  f1 <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = "support", m = "hpa", y = "EPDS", covariates = "age", sex = "sex", B = 9, seed = 1))
  expect_match(f1$rule, "residualised on age \\+ sex \\+ global mean")
  expect_match(f1$rule, "equations carry age \\+ sex \\+ global mean")
  expect_match(f1$global_mean, "methylation columns")
  expect_true("global_mean" %in% f1$terms$Y_model); expect_true("sex" %in% f1$terms$M_model)
  f0 <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = "support", m = "hpa", y = "EPDS", covariates = "age", B = 9, seed = 1, global_mean = FALSE))
  expect_match(f0$rule, "residualised on age; equations carry age$")
  expect_false("global_mean" %in% f0$terms$Y_model)
  expect_false(isTRUE(all.equal(f1$paths[["b[hpa]"]]$b, f0$paths[["b[hpa]"]]$b)))
  ## the technical terms of the measurement rule can differ from the structural covariates
  f2 <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = "support", m = "hpa", y = "EPDS", covariates = "age", measure_covariates = "stress", sex = "sex", B = 9, seed = 1))
  expect_match(f2$rule, "residualised on stress \\+ sex \\+ global mean; equations carry age \\+ sex \\+ global mean")
  ## dmsa_mediate() carries the same rule
  g1 <- suppressWarnings(dmsa_mediate(d$dat, methylation = d$B, map = "alpha", x = "hpa", m = "kynurenine", w = "oxytocin", outcome = "EPDS", covariates = "age", sex = "sex", measurement = "congeneric", B = 9, seed = 1))
  expect_match(g1$rule, "sex \\+ global mean"); expect_output(print(g1), "rule:")
  expect_error(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = "support", m = "hpa", y = "EPDS", sex = "nope", B = 9), "sex column")
})

test_that("the direction source: auto resolves through cpgdirection or the snapshot, and is recorded", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_data_v()
  f <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, x = "support", m = "hpa", y = "EPDS", B = 9, seed = 1, global_mean = FALSE))
  expect_true(is.character(f$map_source))
  if (!is.null(dmsa:::.cpgd("cpg_gene_pairs"))) expect_match(f$map_source, "cpgdirection pairs") else expect_match(f$map_source, "snapshot")
  fa <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = "support", m = "hpa", y = "EPDS", B = 9, seed = 1, global_mean = FALSE))
  expect_match(fa$map_source, "coverage_v4_full")
  expect_error(dmsa_modmed(d$dat, methylation = d$B, map = "nope", x = "support", m = "hpa", y = "EPDS", B = 9), "`map` must be")
})

test_that("connected measurement of two latent mediators (MATH 1.7) is fitted and recorded", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  skip_if_not_installed("lavaan")
  d <- mk_data_v()
  f <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = "support", m = c("hpa", "kynurenine"), y = "EPDS", covariates = "age",
                                    B = 9, seed = 3, global_mean = FALSE, theta_mm = "auto"))
  expect_true(is.list(f$theta_mm)); expect_named(f$theta_mm, "hpa,kynurenine")
  th <- f$theta_mm[[1]]
  expect_true(is.character(th$note)); expect_true(is.logical(th$apply))
  if (is.finite(th$theta)) expect_equal(th$apply, abs(th$ratio) > .10)
  expect_output(print(f), "connected measurement")
})

test_that("the sensitivity battery runs in the call and returns one table", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_data_v()
  f <- suppressWarnings(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = "support", m = "hpa", w = "oxytocin", y = "EPDS", covariates = "age", chip = "chip", block = "fam",
                                    B = 9, seed = 1, sensitivities = c("homoscedastic", "unrestricted", "composite_rel", "quadratic", "global_mean_off"), B_sens = 9))
  sn <- f$sensitivities
  expect_true(is.data.frame(sn)); expect_equal(attr(sn, "B"), 9)
  expect_setequal(unique(sn$model), c("homoscedastic", "unrestricted", "composite_rel", "quadratic", "global_mean_off"))
  expect_true(all(c("IE", "IMM") %in% sn$term[sn$model == "quadratic"]))
  expect_true(all(sn$p[!is.na(sn$p)] > 0 & sn$p[!is.na(sn$p)] <= 1))
  expect_output(print(f), "sensitivities \\(B = 9\\)")
  expect_error(dmsa_modmed(d$dat, methylation = d$B, map = "alpha", x = "support", m = "hpa", y = "EPDS", B = 9, sensitivities = "nope"), "arg")
})
