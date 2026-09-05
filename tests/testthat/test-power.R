# dmsa_power_*(): planning through the package estimator (0.99.15).
# (i) closed forms: monotone in n and in the effect, attenuated by
# reliability, n_for inverts power; (ii) the simulation cell runs through
# .mm_core_v in both modes and holds the level at the null; (iii) the
# hierarchy forms and the surface return coherent numbers.

test_that("closed forms: attenuation, monotonicity, inversion", {
  expect_lt(dmsa_power_slope(.3, 100, .7, .7), dmsa_power_slope(.3, 100, 1, 1))
  expect_lt(dmsa_power_slope(.2, 100), dmsa_power_slope(.3, 100))
  expect_lt(dmsa_power_slope(.3, 100), dmsa_power_slope(.3, 200))
  n <- dmsa_n_slope(.3, .8, .7, .7, k = 2)
  expect_gte(dmsa_power_slope(.3, n, .7, .7, k = 2), .80)
  expect_lt(dmsa_power_slope(.3, n - 5, .7, .7, k = 2), .80)
  expect_equal(dmsa_power_slope(0, 100), .05, tolerance = 1e-6)
})

test_that("the moderated-mediation cell runs through the package core in both modes", {
  p <- dmsa_power_modmed(120, a1 = .4, a3 = .4, b1 = .4, w_on = "a", rel = c(X = .7, M = .9, W = 1), reps = 40, seed = 1)
  expect_s3_class(p, "dmsa_power")
  expect_named(p$power, c("a1", "a3", "b1", "b3", "IE", "IMM"))
  expect_true(is.na(p$power[["b3"]])); expect_true(all(p$power[c("a1", "a3", "b1", "IE", "IMM")] >= 0 & p$power[c("a1", "a3", "b1", "IE", "IMM")] <= 1))
  expect_equal(p$reps_ok, 40)
  expect_equal(unname(p$estimates["truth", ]), c(.4, .4, .4, 0))
  expect_output(print(p), "moderated mediation")
  pp <- dmsa_power_modmed(120, a1 = .4, a3 = .4, b1 = .4, w_on = "a", rel = c(X = .7, M = .9, W = 1), reps = 10, mode = "perm", B = 19, families = TRUE, seed = 1)
  expect_equal(pp$design$mode, "perm"); expect_equal(pp$design$B, 19); expect_true(pp$design$families)
  expect_true(all(pp$power[c("a1", "IE", "IMM")] >= 0))
  ## plain mediation reads IE; second stage names b3
  pm <- dmsa_power_modmed(120, a1 = .4, b1 = .4, reps = 20, moderator = FALSE)
  expect_equal(pm$method, "mediation"); expect_true(is.na(pm$power[["a3"]])); expect_true(is.na(pm$power[["IMM"]])); expect_false(is.na(pm$power[["IE"]]))
  pb <- dmsa_power_modmed(120, a1 = .4, b1 = .4, b3 = .4, w_on = "b", rel = c(X = .7, M = .7, W = .7), reps = 20, quadratic = TRUE, k_cov = 2L)
  expect_true(is.na(pb$power[["a3"]])); expect_false(is.na(pb$power[["b3"]])); expect_true(pb$design$quadratic)
})

test_that("the level holds at the null and a large effect is detected", {
  p0 <- dmsa_power_modmed(150, 0, 0, 0, 0, w_on = "both", rel = c(X = .7, M = .7, W = .7), reps = 300, seed = 7)
  expect_lt(max(p0$power[c("a1", "a3", "b1", "b3")]), .10)
  expect_lt(p0$power[["IMM"]], .06)
  p1 <- dmsa_power_modmed(300, .5, .5, .5, .5, w_on = "both", rel = c(X = .7, M = .7, W = .7), reps = 60, seed = 7)
  expect_gt(p1$power[["IMM"]], .8)
  nn <- dmsa_n_modmed("IMM", .5, n_range = c(60, 400), reps = 60, a1 = .5, a3 = .5, b1 = .5, b3 = .5, w_on = "both", rel = c(X = .7, M = .7, W = .7), seed = 2)
  expect_true(is.finite(nn[["n"]])); expect_gte(nn[["power"]], .5)
})

test_that("hierarchy forms: named power is below raw power and the surface returns its guard", {
  tr <- dmsa_power_transmit(69, .45, .45, units = 741)
  expect_equal(nrow(tr), 2); expect_true(all(tr$power_named < tr$power_raw))
  tr1 <- dmsa_power_transmit(69, .45, .45, units = 1); expect_equal(tr1$power_named, tr1$power_raw)
  co <- dmsa_power_concordance(194, .3, units = 9)
  expect_lt(co[["power_named"]], co[["power_raw"]]); expect_gte(dmsa_power_slope(.3, co[["n_for_80"]], .7, .7), .8)
  ex <- dmsa_power_exposure(192, 0, .3, p1 = .53); expect_true(ex[["power"]] > 0 && ex[["power"]] < 1); expect_gt(ex[["n_for_80"]], 192)
  sf <- dmsa_power_surface(120, .5, reps = 40)
  expect_named(sf, c("power", "guard_rate", "median_a4_ratio", "reps_ok")); expect_true(sf[["power"]] >= 0 && sf[["power"]] <= 1)
})
