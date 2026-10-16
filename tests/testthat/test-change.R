# dmsa_change (the 2-wave identity) and the mod2 three-way through the frame.

.sim_change <- function(n = 160, seed = 5) {
  set.seed(seed)
  n <- 2L * (n %/% 2L)
  G <- 4L; K <- c(6L, 4L, 3L, 5L)
  gene <- rep(paste0("g", 1:G), K); P <- length(gene)
  d <- sample(c(-1, 1), P, TRUE)
  E <- rnorm(n); y <- rnorm(n); cid <- rep(seq_len(n / 2), each = 2)
  M0 <- matrix(rnorm(n * P), n); M1 <- M0 * .6 + matrix(rnorm(n * P), n) * .8
  ## plant: for g1, the CHANGE tracks y x E through the alignment
  idx <- which(gene == "g1")
  for (j in idx) M1[, j] <- M1[, j] + .35 * y * E * d[j]
  al <- dmsa_align(data.frame(cpg = sprintf("cg%04d", 1:P), d = d,
                              p_plus = ifelse(d > 0, .9, .1)),
                   genes = gene, level = "gene")
  list(M0 = M0, M1 = M1, gene = gene, al = al,
       dat = data.frame(y = y, E = E, m2 = rnorm(n), cv = rnorm(n)),
       cid = cid)
}

test_that("dmsa_change finds the planted dS x E unit and adjusts monotonely", {
  s <- .sim_change()
  r <- suppressMessages(
    dmsa_change(s$M0, s$M1, s$dat, outcome = "y", exposure = "E",
                units = s$gene, alignment = s$al, covariates = "cv",
                block = s$cid, B = 199, seed = 1))
  expect_s3_class(r, "dmsa_change")
  expect_true(all(is.finite(r$p)))
  g1 <- r[r$unit == "g1", ]
  expect_true(g1$p_adj < .05)
  expect_true(all(r$p_adj[r$unit != "g1"] > g1$p_adj))
  expect_true(all(diff(r$p_adj[order(r$p)]) >= -1e-12))  # step-down monotone
  expect_output(print(r), "common scale")
})

test_that("dmsa_change runs the outcome role and the mod2 three-way", {
  s <- .sim_change()
  r2 <- suppressMessages(
    dmsa_change(s$M0, s$M1, s$dat, exposure = "E", role = "outcome",
                units = s$gene, alignment = s$al, covariates = "cv",
                block = s$cid, B = 99, seed = 1))
  expect_true(all(is.finite(r2$p)))
  r3 <- suppressMessages(
    dmsa_change(s$M0, s$M1, s$dat, outcome = "y", exposure = "E", mod2 = "m2",
                units = s$gene, alignment = s$al, covariates = "cv",
                block = s$cid, B = 99, seed = 1))
  expect_match(attr(r3, "term"), "mod2")
  expect_true(all(is.finite(r3$p)))
  expect_error(suppressMessages(
    dmsa_change(s$M0, s$M1, s$dat, outcome = "y", exposure = "nope",
                units = s$gene, alignment = s$al, B = 99)), "exposure")
})

test_that("the frame's mod2 three-way detects a planted interaction", {
  set.seed(3)
  n <- 160
  map <- data.frame(gene = rep(c("G1", "G2", "G3"), each = 4))
  map$probe <- sprintf("cg%06d", 1:12)
  map$column <- paste0(map$probe, "_", map$gene)
  map$system_id <- 1L; map$system <- "SimSys"
  map$best_direction <- rep(c(-1, 1), 6)
  map$p_plus <- ifelse(map$best_direction > 0, .9, .1)
  map$best_tier <- "A"; map$smr_tier <- ""
  d <- data.frame(out1 = rnorm(n), m1 = rnorm(n), m2 = rnorm(n),
                  cov1 = rnorm(n), cID = rep(1:(n / 2), each = 2))
  y3 <- scale(d$m1)[, 1] * scale(d$m2)[, 1]
  for (i in seq_len(nrow(map))) {
    base <- rnorm(n)
    if (map$gene[i] == "G1")
      base <- base + .5 * d$out1 * y3 * sign(map$best_direction[i])
    d[[map$column[i]]] <- stats::plogis(base)
  }
  fr <- dmsa_frame(d, map = map, outcome = "out1", covariates = "cov1",
                   blocks = "cID", moderation = TRUE, mod = "m1",
                   mod2 = "m2", B = 199, seed = 1,
                   outdir = file.path(tempdir(), "mod2out"))
  r <- dmsa_report(fr)
  M <- r$moderation
  g1 <- M[M$level == "gene" & M$unit == "G1", ]
  expect_true(nrow(g1) == 1 && g1$p_composite_adj < .05)
  expect_true(all(M[M$level == "gene" & M$unit != "G1", "p_composite"] > .05))
})


test_that("change scores use one common wave scale", {
  set.seed(41)
  M0 <- cbind(rnorm(80), rnorm(80))
  M1 <- M0 + 1.25
  st <- dmsa:::.dmsa_change_stack(M0, M1, winsor = NULL)
  sc <- dmsa:::.dmsa_change_project(st, c(1, 1))
  expect_gt(mean(sc$S1 - sc$S0), 0.5)
  ## A separate within-wave z transform would erase this common shift.
  old0 <- as.numeric(scale(rowMeans(scale(M0))))
  old1 <- as.numeric(scale(rowMeans(scale(M1))))
  expect_lt(abs(mean(old1 - old0)), 1e-10)
})

test_that("fixed within-pair roles are refused outside paired mode", {
  set.seed(42)
  np <- 24; n <- 2 * np
  gene <- rep(c("g1", "g2"), each = 3); P <- length(gene)
  dir <- rep(c(1, -1, 1), 2)
  M0 <- matrix(rnorm(n * P), n)
  M1 <- .7 * M0 + matrix(rnorm(n * P), n) * .7
  al <- dmsa_align(data.frame(cpg = paste0("cg", seq_len(P)), d = dir,
                              p_plus = ifelse(dir > 0, .9, .1)),
                   genes = gene, level = "gene")
  dat <- data.frame(mother = rep(c(0, 1), np), cid = rep(seq_len(np), each = 2))
  expect_error(
    suppressMessages(dmsa_change(M0, M1, dat, exposure = "mother",
      role = "outcome", units = gene, alignment = al, block = dat$cid,
      B = 49, seed = 1)),
    "fixed two-level within-block role")
})

test_that("paired change mode detects a planted focal-role program", {
  set.seed(43)
  np <- 60; n <- 2 * np
  gene <- rep(c("g1", "g2", "g3"), each = 4); P <- length(gene)
  dir <- rep(c(1, -1, 1, -1), 3)
  role <- rep(c(0, 1), np)
  cid <- rep(seq_len(np), each = 2)
  M0 <- matrix(rnorm(n * P), n)
  M1 <- .65 * M0 + matrix(rnorm(n * P), n) * .75
  j1 <- which(gene == "g1")
  M1[role == 1, j1] <- M1[role == 1, j1, drop = FALSE] +
    matrix(rep(.85 * dir[j1], each = sum(role == 1)), nrow = sum(role == 1))
  al <- dmsa_align(data.frame(cpg = paste0("cg", seq_len(P)), d = dir,
                              p_plus = ifelse(dir > 0, .95, .05)),
                   genes = gene, level = "gene")
  dat <- data.frame(mother = role, cid = cid, dcell = rnorm(n, 0, .05))
  r <- suppressMessages(dmsa_change(
    M0, M1, dat, exposure = "mother", role = "outcome",
    units = gene, alignment = al, covariates = "dcell",
    block = dat$cid, paired = TRUE, pair_focal = 1,
    weighting = "flat", B = 399, seed = 7))
  expect_s3_class(r, "dmsa_change")
  expect_true(isTRUE(attr(r, "paired")))
  expect_match(attr(r, "inference"), "wild bootstrap")
  expect_false(isTRUE(attr(r, "baseline_adjusted")))
  expect_lt(r$p_adj[r$unit == "g1"], .05)
  expect_true(all(diff(r$p_adj[order(r$p)]) >= -1e-12))
})

test_that("outcome-side reliability weighting is refused until calibrated", {
  s <- .sim_change()
  expect_error(
    suppressMessages(dmsa_change(s$M0, s$M1, s$dat, exposure = "E",
      role = "outcome", units = s$gene, alignment = s$al,
      weighting = "reliability", B = 19)),
    "requires weighting = 'flat'")
})

test_that("duplicate CpG x unit rows are refused", {
  M0 <- matrix(rnorm(40), 10, 4); M1 <- matrix(rnorm(40), 10, 4)
  units <- rep("g1", 4)
  al <- dmsa_align(data.frame(cpg = c("cg1", "cg1", "cg2", "cg3"),
                              d = c(1, 1, -1, 1),
                              p_plus = c(.9, .9, .1, .9)),
                   genes = units, level = "gene")
  dat <- data.frame(y = rnorm(10), E = rnorm(10))
  expect_error(dmsa_change(M0, M1, dat, outcome = "y", exposure = "E",
                           units = units, alignment = al, B = 0),
               "same CpG appears more than once")
})
