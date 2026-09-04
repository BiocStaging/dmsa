# dmsa_transmit() / dmsa_concordance(): the full-hierarchy path functions
# (0.99.14). (i) the shared-stream statistic .hier_fl() reproduces the
# engine's .med_fl() exactly (b, t and p under the same seed); (ii) the two
# functions run on synthetic frames at every level and return coherent
# tables (maxT p >= raw p, p_tree >= maxT p); (iii) the refusals (batch has
# no default; roles); (iv) dmsa_contract() carries the rule.

mk_hier <- function(n_fam = 40, seed = 3) {
  set.seed(seed)
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv", package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction) & mp$system_id %in% c(1, 2), ]
  pr <- unique(mp$probe)
  ## parents: 2 per family (sex 2 = woman, 1 = man), on 8-sample chips, on 3 plates
  np <- 2 * n_fam
  P <- matrix(plogis(rnorm(np * length(pr), sd = 1.2)), np, length(pr), dimnames = list(NULL, pr))
  par <- data.frame(cID = rep(seq_len(n_fam), each = 2), sex = rep(c(2, 1), n_fam), age = rnorm(np), Epi = runif(np),
                    chip = factor(rep(seq_len(np / 8), each = 8)), plate = factor(rep(1:3, length.out = n_fam)[rep(seq_len(n_fam), each = 2)]),
                    y = rnorm(np), pill = rep(rbinom(n_fam, 1, .5), each = 2))
  ## children: one per family, share a little of the mother's methylation
  C <- matrix(plogis(rnorm(n_fam * length(pr), sd = 1.2)), n_fam, length(pr), dimnames = list(NULL, pr))
  C <- plogis(qlogis(C) + 0.4 * qlogis(P[seq(1, np, 2), ]))
  chi <- data.frame(cID = seq_len(n_fam), sex = rbinom(n_fam, 1, .5) + 1, bw = rnorm(n_fam), Epi = runif(n_fam),
                    chip = factor(rep(seq_len(ceiling(n_fam / 8)), each = 8)[seq_len(n_fam)]))
  list(P = P, par = par, C = C, chi = chi)
}

test_that(".hier_fl reproduces .med_fl exactly on a shared stream", {
  set.seed(7); n <- 80; fam <- factor(rep(seq_len(40), each = 2))
  D <- scale(cbind(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n), c1 = rnorm(n)), scale = FALSE)
  mm <- dmsa:::.med_moments(D, voffset = c(x1 = 0.3), products = list(), n = n)
  set.seed(11); r1 <- dmsa:::.med_fl(D, mm$S, c("x1", "x2", "c1"), "y", "x1", n, 99, fam, Vy = 0, block_mode = "exchange", Sobs = mm$Sobs, pivot = "hc")
  perms <- dmsa:::.hier_perms(fam, n, 99, seed = 11)
  r2 <- dmsa:::.hier_fl(D, mm$S, mm$Sobs, c("x1", "x2", "c1"), "y", "x1", n, perms)
  expect_equal(r2$b, r1$b, tolerance = 1e-12)
  expect_equal(r2$t, r1$t, tolerance = 1e-10)
  expect_identical((1 + sum(abs(r2$tn) >= abs(r2$t))) / 100, r1$p)
  ## the contrast form with a unit vector is the single-term form
  cv <- c(x1 = 1, x2 = 0, c1 = 0)
  r3 <- dmsa:::.hier_fl(D, mm$S, mm$Sobs, c("x1", "x2", "c1"), "y", NULL, n, perms, contrast = cv)
  expect_equal(r3$t, r2$t, tolerance = 1e-12); expect_equal(r3$tn, r2$tn, tolerance = 1e-12)
  ## maxT >= raw p, and both in (0, 1]
  TN <- cbind(r2$tn, r3$tn); tt <- c(r2$t, r3$t)
  expect_true(all(dmsa:::.hier_maxT(tt, TN) >= dmsa:::.hier_raw_p(tt, TN)))
})

test_that("dmsa_transmit runs the four levels and returns coherent tables", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_hier()
  child <- dmsa_contract(d$chi, methylation = d$C, covariates = "bw", measure_covariates = "Epi", sex = "sex", chip = "chip")
  parents <- dmsa_contract(d$par, methylation = d$P, covariates = "age", measure_covariates = "Epi", sex = "sex", chip = "chip")
  tr <- suppressWarnings(dmsa_transmit(child, parents, link = "cID", parent_role = "sex", roles = c(mother = 2, father = 1),
                                       systems = "oxytocin", map = "alpha", B = 19, seed = 1))
  expect_s3_class(tr, "dmsa_transmit")
  expect_setequal(unique(tr$units$level), c("probe", "gene", "module", "system"))
  expect_equal(tr$design$n, 40)
  u <- tr$units
  expect_true(all(u$p_mother >= 1 / 20 & u$p_mother <= 1))
  expect_true(all(u$p_mother_maxT >= u$p_mother - 1e-12)); expect_true(all(u$p_father_maxT >= u$p_father - 1e-12))
  expect_true(all(u$p_mother_tree >= u$p_mother_maxT - 1e-12)); expect_true(all(u$p_father_tree >= u$p_father_maxT - 1e-12))
  expect_true(all(u$shrink[u$level == "probe"] == 1))
  expect_match(tr$design$rule_child, "sex \\+ global mean")
  expect_equal(colnames(tr$units)[1:3], c("system", "level", "unit"))
  expect_true(all(c("global_mean_child", "global_mean_mother", "global_mean_father") %in% tr$design$covariates))
  expect_output(print(tr), "Parent -> child transmission")
  ## the planted maternal resemblance is visible at the probe level: mother slopes exceed father slopes on average
  pr <- u[u$level == "probe", ]
  expect_gt(mean(pr$b_mother), mean(pr$b_father))
  ## a level subset and the rule switches
  tr2 <- suppressWarnings(dmsa_transmit(dmsa_contract(d$chi, methylation = d$C, covariates = "bw", global_mean = FALSE),
                                        dmsa_contract(d$par, methylation = d$P, covariates = "age", global_mean = FALSE),
                                        systems = "oxytocin", levels = c("module", "system"), map = "alpha", B = 9))
  expect_setequal(unique(tr2$units$level), c("module", "system"))
  expect_false(any(grepl("global_mean", tr2$design$covariates)))
  expect_error(dmsa_transmit(child, parents, roles = c(mum = 2, dad = 1), systems = "oxytocin", map = "alpha", B = 9), "mother and father")
})

test_that("dmsa_concordance: batch has no default, step 0 at every level, outcome surface and exposure form", {
  skip_if_not(nzchar(system.file("extdata", "coverage_v4_full.csv", package = "dmsa")))
  d <- mk_hier()
  parents <- dmsa_contract(d$par, methylation = d$P, covariates = "age", measure_covariates = "Epi", sex = "sex", chip = "chip")
  expect_error(dmsa_concordance(parents, couple = "cID", partner_by = "sex", roles = c(woman = 2, man = 1), systems = "oxytocin", map = "alpha", B = 9, NP = 20), "no default")
  s0 <- suppressWarnings(dmsa_concordance(parents, couple = "cID", partner_by = "sex", roles = c(woman = 2, man = 1), systems = "oxytocin",
                                          batch = "plate", map = "alpha", B = 9, NP = 39, seed = 1))
  expect_s3_class(s0, "dmsa_concordance")
  expect_equal(s0$design$couples, 40); expect_equal(s0$design$batches, 3)
  expect_setequal(unique(s0$step0$level), c("probe", "gene", "module", "system"))
  expect_true(all(s0$step0$p_maxT >= s0$step0$p - 1e-12)); expect_true(all(s0$step0$p_tree >= s0$step0$p_maxT - 1e-12))
  expect_true(all(abs(s0$step0$r) <= 1))
  expect_output(print(s0), "Dyadic concordance")
  ## no shared batch declared: one group
  s1 <- suppressWarnings(dmsa_concordance(parents, couple = "cID", partner_by = "sex", roles = c(woman = 2, man = 1), systems = "oxytocin",
                                          batch = FALSE, levels = "module", map = "alpha", B = 9, NP = 19))
  expect_equal(s1$design$batches, 1); expect_match(s1$design$batch, "none declared")
  ## outcome surface and exposure form
  s2 <- suppressWarnings(dmsa_concordance(parents, couple = "cID", partner_by = "sex", roles = c(woman = 2, man = 1), systems = "oxytocin",
                                          batch = "plate", outcome = "y", exposure = "pill", exposed = "woman", levels = c("module", "system"),
                                          map = "alpha", B = 19, NP = 19, seed = 2))
  expect_true(is.data.frame(s2$surface)); expect_true(all(c("a4", "a4_uncorrected", "p_maxT", "p_uncorrected_maxT", "p_tree") %in% names(s2$surface)))
  expect_true(all(s2$surface$p_maxT >= s2$surface$p - 1e-12))
  expect_true(is.data.frame(s2$exposure)); expect_true(all(c("gamma3", "p_maxT", "p_tree") %in% names(s2$exposure)))
  expect_equal(attr(s2$exposure, "design")$exposed_role, "woman")
  expect_output(print(s2), "surface -> y"); expect_output(print(s2), "exposure pill of the woman")
  expect_error(dmsa_concordance(parents, couple = "cID", partner_by = "sex", roles = c(woman = 2, man = 1), systems = "oxytocin",
                                batch = "plate", exposure = "pill", exposed = "nobody", map = "alpha", B = 9, NP = 9), "exposed")
})

test_that("dmsa_contract carries the rule and refuses the wrong input", {
  ct <- dmsa_contract(data.frame(a = 1), covariates = "a", sex = "s")
  expect_s3_class(ct, "dmsa_contract"); expect_true(ct$global_mean); expect_equal(ct$tissue, "blood")
  expect_error(dmsa:::.as_contract(1:3, "child"), "dmsa_contract")
})

test_that("the path functions resolve the direction map through cpgdirection when it is installed", {
  skip_if_not(!is.null(dmsa:::.cpgd("cpg_gene_pairs")))
  d <- mk_hier(n_fam = 20)
  st <- dmsa:::.paths_setup(d$par, d$P, map = "auto", systems = "oxytocin", covariates = "age", chip = "chip")
  expect_match(st$map_source, "cpgdirection pairs")
  expect_match(st$map_source, "cascade pairs only")
  expect_true(all(c("probe", "gene", "best_direction", "p_plus") %in% names(st$mp)))
})

test_that("the regenerated cpgdirection snapshot is the fallback of the path functions", {
  sn <- dmsa:::.paths_snapshot("blood")
  skip_if(is.null(sn))
  expect_true(all(c("probe", "column", "gene", "system_id", "system", "module_id", "best_direction", "p_plus", "tissue") %in% names(sn)))
  expect_true(all(sn$tissue == "blood")); expect_gt(nrow(sn), 9000)
  expect_true(all(sn$best_direction %in% c(-1, 1)))
  sn2 <- dmsa:::.paths_snapshot("nasal_epithelium"); expect_true(all(sn2$tissue == "nasal_epithelium"))
  cf <- dmsa:::.paths_snapshot("blood", "confidence"); expect_true(all(cf$direction_tier %in% c("M", "S1")))
})
