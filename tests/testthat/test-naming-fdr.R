## PI ruling 2026-10-15 (ERRCTL_DECISION_v2.md): a unit is NAMED when its
## dependence-calibrated Benjamini-Hochberg q-value on the unit p is at or
## below `fdr` within its level-local family. The pre-registered any-lens
## rule is kept and reported beside it; EXACT-CONFIRMED means the joint
## union test (p_unit_adj) also clears alpha.

.fdr_fixture <- function(n = 60, k = 3, genes = "FKBP5") {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene %in% genes], k)
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  list(data = d, cols = cols, n = n)
}

.null_family <- function(B, U, rho, seed) {
  set.seed(seed)
  C <- matrix(rho, U, U); diag(C) <- 1
  abs(matrix(stats::rnorm(B * U), B, U) %*% chol(C))
}
.p_against <- function(obs, nul)            # same convention as p_unit
  vapply(seq_along(obs), function(u)
    (1 + sum(nul[, u] >= obs[u])) / (nrow(nul) + 1), numeric(1))

test_that("q_unit is BH, q_unit_dep is monotone in it, NA stays NA", {
  nul <- .null_family(999, 6, .4, 1)
  p <- c(.001, .02, NA, .2, .5, .9)
  qd <- .q_unit_dep(p, nul)
  ok <- is.finite(p)
  expect_equal(qd$q[ok], stats::p.adjust(p[ok], "BH"))
  expect_true(is.na(qd$q[3]) && is.na(qd$q_dep[3]))
  expect_true(all(qd$q_dep[ok] > 0 & qd$q_dep[ok] <= 1))
  o <- order(qd$q[ok])
  expect_true(all(diff(qd$q_dep[ok][o]) >= 0))
  expect_length(qd$qmin, 999)
  ## no testable unit: everything NA, nothing named
  q0 <- .q_unit_dep(c(NA_real_, NA_real_), nul[, 1:2])
  expect_true(all(is.na(q0$q_dep)))
})

test_that("calibrated naming fires in a fraction fdr of null families, under dependence", {
  for (rho in c(0, .8)) {
    nul <- .null_family(499, 12, rho, 7)
    obs <- .null_family(300, 12, rho, 8)
    fired <- vapply(seq_len(nrow(obs)), function(i) {
      qd <- .q_unit_dep(.p_against(obs[i, ], nul), nul)
      any(qd$q_dep <= .10)
    }, logical(1))
    ## 300 null families: binomial SE ~ .017 around .10 (bounds ~3 SE)
    expect_gt(mean(fired), .045)
    expect_lt(mean(fired), .155)
  }
})

test_that("dmsa_triangulate returns q_unit and q_unit_dep for every testable unit", {
  set.seed(1)
  n <- 100
  d <- rep(c(1, -1), each = 3)
  f <- rnorm(n)
  M <- sapply(rep(d, 2), function(dj) dj * .5 * f + sqrt(1 - .5^2) * rnorm(n))
  colnames(M) <- paste0("cg", 1:12)
  dat <- data.frame(y = .4 * f + rnorm(n), cv = rnorm(n))
  units <- rep(c("G1", "G2"), each = 6)
  al <- dmsa_align(data.frame(cpg = colnames(M), d = rep(d, 2),
                              p_plus = ifelse(rep(d, 2) > 0, .9, .1)),
                   genes = units, level = "gene")
  r <- dmsa_triangulate(M, dat, rhs = c("y", "cv"), term = "y",
                        units = units, alignment = al, B = 99, seed = 1)
  expect_true(all(c("q_unit", "q_unit_dep") %in% names(r)))
  ok <- r$n_probes > 0
  expect_true(all(is.finite(r$q_unit_dep[ok])))
  expect_true(all(r$q_unit[ok] >= r$p_unit[ok] - 1e-12))
  expect_equal(r$q_unit[ok], stats::p.adjust(r$p_unit[ok], "BH"))
  expect_length(attr(r, "unit_null_qmin"), 99)
})

test_that("the report names by the calibrated q and keeps the pre-registered rule beside it", {
  set.seed(1); fx <- .fdr_fixture()
  od <- tempfile("dmsa_fdr_rule")
  f <- dmsa_frame(fx$data, methylation = fx$cols, direction_source = "bundled",
                  outcome = "y", covariates = "cov1", blocks = "cID",
                  chip = FALSE, B = 49, fdr = .10, plots = FALSE,
                  tables = TRUE, summary = TRUE, progress = FALSE,
                  beep = FALSE, table_type = "html", outdir = od)
  expect_equal(f$fdr, .10)
  r <- suppressWarnings(suppressMessages(dmsa_report(f)))
  u <- r$results
  expect_true(all(c("q_unit", "q_unit_dep", "prereg_any_lens",
                    "null_alarm_bh") %in% names(u)))
  expect_identical(u$selected,
                   is.finite(u$q_unit_dep) & u$q_unit_dep <= .10 &
                     u$n_probes > 0)
  expect_identical(u$prereg_any_lens, u$any_lens_hit)
  expect_true(all(!u$exact_confirmed | (u$selected & u$p_unit_adj < .05)))
  sm <- readLines(file.path(od, "summary.md"), warn = FALSE)
  expect_true(any(grepl("Pre-registered any-lens rule", sm, fixed = TRUE)))
  expect_true(any(grepl("`q_unit_dep`", sm, fixed = TRUE)))
  ut <- utils::read.csv(file.path(od, "tables", "units.csv"))
  expect_true(all(c("q_unit_dep", "prereg_any_lens") %in% names(ut)))
})

test_that("fdr is validated", {
  set.seed(1); fx <- .fdr_fixture()
  expect_error(dmsa_frame(fx$data, methylation = fx$cols,
                          direction_source = "bundled", outcome = "y",
                          covariates = "cov1", blocks = "cID", chip = FALSE,
                          fdr = 1.5, progress = FALSE, beep = FALSE,
                          outdir = tempfile("x")), "`fdr`")
})
