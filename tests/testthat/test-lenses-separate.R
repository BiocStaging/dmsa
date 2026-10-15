## Spec 63: the three lenses tested SEPARATELY.
##
## Every other test exercises the lenses together, so a regression that made
## one lens return another's statistic - or swapped two labels - would pass
## as long as the any-lens union still fired, while every report attributed
## its findings to the wrong lens. These tests plant a signal whose SHAPE
## picks out a lens by that lens's own definition (R/triangulate.R header):
##
##   diffuse      directionless quadratic - fires on signal spread over the
##                probes whichever way each one points
##   coherence /  the aligned linear contrast, with 1/se^2 (coherence) or
##   composite    equal (composite) probe weights - both die when the aligned
##                signs cancel, both gain sqrt(k) on dense consistent signal
##
## Shapes were calibrated by simulation (dmsa_patch/undermine_lenses*.R,
## 2026-09-02): "paired cancel" = every probe tracks the outcome but the
## aligned signs cancel pairwise (diffuse only); "dense tiny consistent" =
## 60 probes each carrying an aligned effect too small for the quadratic
## lens to see (directional only). Flat weighting is used because the
## reliability weights can break an exact sign cancellation (reported to the
## PI as E12).

.lens_shape <- function(kind, seed) {
  set.seed(seed)
  if (kind == "diffuse") {
    n <- 160; k <- 10
    d <- rep(c(1, -1), length.out = k)          # mixed direction calls
    sg <- rep(c(1, -1), each = k / 2)           # aligned signs cancel pairwise
    y <- stats::rnorm(n); cov1 <- stats::rnorm(n)
    M <- sapply(seq_len(k), function(j) d[j] * sg[j] * 0.3 * y + stats::rnorm(n))
  } else if (kind == "directional") {
    n <- 200; k <- 60
    d <- rep(c(1, -1), length.out = k)
    y <- stats::rnorm(n); cov1 <- stats::rnorm(n)
    M <- sapply(seq_len(k), function(j) d[j] * 0.03 * y + stats::rnorm(n))
  } else {
    n <- 160; k <- 10
    d <- rep(c(1, -1), length.out = k)
    y <- stats::rnorm(n); cov1 <- stats::rnorm(n)
    M <- sapply(seq_len(k), function(j) stats::rnorm(n))
  }
  colnames(M) <- sprintf("cg%07d", seq_len(k))
  list(M = M, data = data.frame(y = y, cov1 = cov1), d = d, k = k)
}

.lens_tri <- function(s, B = 199, seed = 1) {
  al <- dmsa_align(data.frame(cpg = colnames(s$M), d = s$d,
                              p_plus = ifelse(s$d > 0, .95, .05)),
                   genes = rep("G1", s$k), level = "gene")
  dmsa_triangulate(s$M, s$data, rhs = c("y", "cov1"), term = "y",
                   units = rep("G1", s$k), alignment = al, B = B, seed = seed,
                   weighting = "flat")
}

test_that("spec 63: a sign-cancelling signal is seen by the diffuse lens alone", {
  r <- .lens_tri(.lens_shape("diffuse", seed = 2))
  expect_lt(r$p_diffuse, 0.05)
  expect_gt(r$p_coherence, 0.05)
  expect_gt(r$p_composite, 0.05)
  expect_equal(which.min(c(r$p_coherence, r$p_composite, r$p_diffuse)), 3L)
})

test_that("spec 63: a dense, tiny, consistent signal is seen by the directional lenses, not the diffuse one", {
  r <- .lens_tri(.lens_shape("directional", seed = 4))
  expect_lt(r$p_coherence, 0.05)
  expect_lt(r$p_composite, 0.05)
  expect_gt(r$p_diffuse, 0.05)
})

test_that("spec 63: pure noise fires no lens", {
  r <- .lens_tri(.lens_shape("null", seed = 7))
  expect_gt(r$p_coherence, 0.05)
  expect_gt(r$p_composite, 0.05)
  expect_gt(r$p_diffuse, 0.05)
})

## ---- report level: the carrying lens is the lens that fired --------------
## The Results sentence says "carried by the <lens> lens". These pin that the
## label, the units table and the prose all point at the lens the planted
## shape targets - the sentence that goes into papers.

.lens_frame <- function(kind, seed, od) {
  set.seed(seed)
  if (kind == "diffuse") { n <- 160; k <- 10; eff <- 0.3
    sg <- rep(c(1, -1), each = k / 2) } else { n <- 200; k <- 60; eff <- 0.035
    sg <- rep(1, k) }
  d <- rep(c(1, -1), length.out = k)
  dat <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                    cID = rep(seq_len(n / 2), each = 2))
  pr <- sprintf("cg%07d", seq_len(k + 6))
  genes <- c(rep("NR3C1", k), rep("FKBP5", 3), rep("CRH", 3))
  dd <- c(d, rep(c(1, -1, 1), 2))
  ys <- scale(dat$y)[, 1]
  for (j in seq_len(k))
    dat[[pr[j]]] <- stats::plogis(dd[j] * sg[j] * eff * ys + stats::rnorm(n))
  for (j in (k + 1):(k + 6)) dat[[pr[j]]] <- stats::plogis(stats::rnorm(n))
  pairs <- data.frame(cpg_id = pr, target_gene = genes, best_direction = dd,
                      probability_plus1 = ifelse(dd > 0, .95, .05),
                      usable = TRUE, abstain_reason = NA_character_,
                      best_evidence = "smr_high", direction_tier = "S1",
                      stringsAsFactors = FALSE)
  old <- options(dmsa.pair_table = pairs)
  f <- suppressMessages(dmsa_frame(
    dat, methylation = pr, direction_source = "cpgdirection",
    outcome = "y", covariates = "cov1", blocks = "cID",
    chip = FALSE, weighting = "flat", B = 199, progress = FALSE,
    beep = FALSE, outdir = od))
  options(old)
  f
}

test_that("spec 63: a diffuse-only gene is NAMED and reported as carried by the diffuse lens", {
  od <- tempfile("dmsa_lens63_diffuse")
  f <- .lens_frame("diffuse", seed = 11, od)
  r <- suppressWarnings(suppressMessages(dmsa_report(f)))
  u <- r$results[r$results$level == "gene" & r$results$unit == "NR3C1", ]
  expect_true(isTRUE(u$selected))
  expect_identical(u$best_lens, "diffuse")
  expect_lt(u$p_diffuse_adj, 0.05)
  expect_gt(u$p_coherence_adj, 0.05)
  expect_gt(u$p_composite_adj, 0.05)
  sm <- readLines(file.path(od, "summary.md"))
  ## naming sentence re-ruled 2026-10-15: "NR3C1 (y, q = ..., unit p ...,
  ## carried by the <lens> lens; ...)"
  expect_true(any(grepl("NR3C1 \\(y, q = [0-9.]+, unit p [0-9.]+, carried by the diffuse lens", sm)))
  ## the units table on disk agrees with the object
  ut <- utils::read.csv(file.path(od, "tables", "units.csv"))
  expect_identical(ut$best_lens[ut$level == "gene" & ut$unit == "NR3C1"], "diffuse")
  ## spec 62 badge pin: a confirmed gene carries its dagger(s) in the gene
  ## table and the headline says so - as real UTF-8 characters whatever the
  ## session locale (summary.md used to hold "<U+2020>" under C/POSIX)
  raw <- readLines(file.path(od, "summary.md"), warn = FALSE, encoding = "UTF-8")
  expect_true(all(validUTF8(raw)))
  expect_false(any(grepl("<U+20", raw, fixed = TRUE)))
  row <- grep("^\\| \\*\\*NR3C1\\*\\*", raw, value = TRUE)
  expect_length(row, 1L)
  if (isTRUE(u$exact_confirmed)) {
    expect_true(grepl("\u2020", row, useBytes = FALSE))
    expect_true(any(grepl("EXACT-CONFIRMED", raw)))
  }
  if (isTRUE(u$omnibus_confirmed)) {
    expect_true(grepl("\u2021", row))
    expect_true(any(grepl("OMNIBUS-CONFIRMED", raw)))
  }
  expect_true(isTRUE(u$exact_confirmed) || isTRUE(u$omnibus_confirmed))
})

test_that("spec 63: a dense consistent gene is NAMED by a directional lens, not the diffuse one", {
  od <- tempfile("dmsa_lens63_dir")
  f <- .lens_frame("directional", seed = 3, od)
  r <- suppressWarnings(suppressMessages(dmsa_report(f)))
  u <- r$results[r$results$level == "gene" & r$results$unit == "NR3C1", ]
  expect_true(isTRUE(u$selected))
  expect_true(u$best_lens %in% c("coherence", "composite"))
  expect_lt(min(u$p_coherence_adj, u$p_composite_adj), 0.05)
  expect_gt(u$p_diffuse_adj, 0.05)
  sm <- readLines(file.path(od, "summary.md"))
  expect_true(any(grepl(sprintf("NR3C1 \\(y, q = [0-9.]+, unit p [0-9.]+, carried by the %s lens",
                                u$best_lens), sm)))
  expect_false(any(grepl("NR3C1 \\(y, q = [0-9.]+, unit p [0-9.]+, carried by the diffuse lens", sm)))
})
