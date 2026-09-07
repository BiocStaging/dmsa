## A small synthetic HPA frame built from the bundled cascade, so that the
## path functions see real system / module / gene names.
## Directions repeat per module (four of five agree), so the majority builder
## - the default since 0.99.17 - keeps enough probes to score every module.
.rp_fixture <- function(n = 80L, n_mod = 2L, per_mod = 5L, suffix = FALSE,
                        seed = 3L) {
  set.seed(seed)
  cas <- dmsa_sets("alpha")$cascade
  h <- cas[cas$system_short == "hpa", ]
  mods <- utils::head(unique(h$module_id), n_mod)
  pick <- do.call(rbind, lapply(mods, function(m)
    utils::head(unique(h[h$module_id == m, c("cpg", "gene", "module_id")]),
                per_mod)))
  mm <- unique(h[, c("module_id", "module")])
  map <- data.frame(system_id = 2, system = "HPA axis & glucocorticoid signalling",
                    probe = pick$cpg, column = pick$cpg, gene = pick$gene,
                    module_id = pick$module_id,
                    module = mm$module[match(pick$module_id, mm$module_id)],
                    best_direction = rep(c(-1, -1, -1, -1, 1), length.out = nrow(pick)),
                    p_plus = rep(c(.1, .1, .1, .1, .9), length.out = nrow(pick)),
                    stringsAsFactors = FALSE)
  dat <- data.frame(x = stats::rnorm(n), cov1 = stats::rnorm(n),
                    cID = rep(seq_len(n / 2), each = 2))
  eta <- stats::rnorm(n)
  for (i in seq_len(nrow(map)))
    dat[[map$column[i]]] <- stats::plogis(
      stats::rnorm(n) + .9 * eta * map$best_direction[i])
  dat$y <- .5 * eta + .3 * dat$x + stats::rnorm(n)
  if (suffix)
    names(dat)[names(dat) %in% map$column] <- paste0(map$column, "_TC21_T1_G")
  list(data = dat, map = map)
}

test_that("dmsa_report() writes tables, figures and summary.md for a path method", {
  f <- .rp_fixture()
  r <- dmsa_modmed(f$data, map = f$map, x = "x", m = "hpa", y = "y",
                   covariates = "cov1", block = "cID", B = 19L, seed = 1,
                   single_rel = .7, global_mean = FALSE)
  expect_s3_class(r, "dmsa_modmed")

  out <- tempfile()
  rep <- dmsa_report(r, outdir = out)
  expect_s3_class(rep, "dmsa_path_report")
  expect_true(file.exists(file.path(out, "summary.md")))
  expect_true(file.exists(file.path(out, "tables", "paths.csv")))
  expect_true(length(list.files(file.path(out, "figures"),
                                pattern = "\\.png$")) >= 1L)

  ## the summary opens with the call that produced it, and carries the
  ## numbers the object already holds - nothing is recomputed
  md <- readLines(file.path(out, "summary.md"))
  expect_true(any(grepl("dmsa_modmed", md)))
  tb <- utils::read.csv(file.path(out, "tables", "paths.csv"))
  expect_equal(nrow(tb), length(r$paths))
  expect_equal(tb$b[match(names(r$paths)[1], tb$path)],
               r$paths[[1]]$b, tolerance = 1e-8)

  ## an existing report is never silently replaced, and a refusal leaves
  ## what was there untouched
  n_before <- length(list.files(out, recursive = TRUE))
  expect_error(dmsa_report(r, outdir = out), "already contains a DMSA report")
  expect_equal(length(list.files(out, recursive = TRUE)), n_before)
  expect_s3_class(dmsa_report(r, outdir = out, overwrite = TRUE),
                  "dmsa_path_report")
})

test_that("the registered defaults are the shipped defaults", {
  ## PI ruling 2026-09-07: a single-indicator tone carries the FIXED .70, and
  ## the module tone is the majority builder - the registered recipe. Before
  ## 0.99.17 the reliability was estimated and floored at .05 on real data,
  ## which made every path through that tone unreadable.
  expect_equal(formals(dmsa_modmed)$single_rel, 0.70)
  expect_equal(formals(dmsa_mediate)$single_rel, 0.70)
  expect_equal(eval(formals(dmsa_modmed)$module_score)[1], "majority")
  expect_equal(eval(formals(dmsa_transmit)$module_score)[1], "majority")
  expect_equal(eval(formals(dmsa_concordance)$module_score)[1], "majority")

  f <- .rp_fixture()
  r <- dmsa_modmed(f$data, map = f$map, x = "x", m = "hpa", y = "y",
                   block = "cID", B = 19L, seed = 1, global_mean = FALSE)
  expect_equal(r$measurement$hpa$omega, 0.70, tolerance = 1e-8)

  ## NA restores the estimated reliability
  r2 <- dmsa_modmed(f$data, map = f$map, x = "x", m = "hpa", y = "y",
                    block = "cID", B = 19L, seed = 1, global_mean = FALSE,
                    single_rel = NA)
  expect_false(isTRUE(all.equal(r2$measurement$hpa$omega, 0.70)))
})

test_that("the paths table carries the standardized path beside the raw one", {
  f <- .rp_fixture()
  r <- dmsa_modmed(f$data, map = f$map, x = "x", m = "hpa", y = "y",
                   block = "cID", B = 19L, seed = 1, global_mean = FALSE)
  d <- dmsa:::.mr_paths_df(r)
  expect_true("b_std" %in% names(d))
  ## a latent role has corrected variance 1, so b[m] standardizes by the
  ## outcome's SD alone
  i <- match("b[hpa]", d$path)
  expect_equal(d$b_std[i], d$b[i] / r$scales[["y"]], tolerance = 1e-10)
  ## and an observed-to-observed path by both
  j <- match("c'[x]", d$path)
  expect_equal(d$b_std[j], d$b[j] * r$scales[["x"]] / r$scales[["y"]],
               tolerance = 1e-10)
})

test_that("dmsa_report() dispatches, and refuses what it cannot report", {
  expect_error(dmsa_report(42), "dmsa_frame or the result of")
  expect_true("dmsa_report.dmsa_frame" %in%
                as.character(utils::methods("dmsa_report")))
})

test_that("methylation columns are matched by canonical CpG id", {
  ## a frame whose columns carry a study suffix (cg01234567_TC21_T1_GENE)
  ## must reach the same pairs as one whose columns are bare CpG ids: before
  ## 0.99.17 the suffixed frame matched no pair at all and every path method
  ## stopped at "no usable module tone".
  a <- .rp_fixture(suffix = FALSE)
  b <- .rp_fixture(suffix = TRUE)     # same seed, same numbers, other names
  st_bare <- dmsa:::.paths_setup(a$data, NULL, map = a$map, systems = "hpa",
                                 columns = c("x", "y"), block = "cID",
                                 global_mean = FALSE)
  st_suff <- dmsa:::.paths_setup(b$data, NULL, map = b$map, systems = "hpa",
                                 columns = c("x", "y"), block = "cID",
                                 global_mean = FALSE)
  expect_equal(colnames(st_suff$ME), colnames(st_bare$ME))
  expect_equal(unname(st_suff$tones[["hpa"]]),
               unname(st_bare$tones[["hpa"]]), tolerance = 1e-10)

  ## the two designs of one CpG are one measurement: averaged, not counted
  ## twice (the frame's collapse_site policy)
  d <- a$data
  first <- a$map$column[1]
  d[[paste0(first, "_B")]] <- d[[first]]
  names(d)[names(d) == first] <- paste0(first, "_A")
  st_dup <- dmsa:::.paths_setup(d, NULL, map = a$map, systems = "hpa",
                                columns = c("x", "y"), block = "cID",
                                global_mean = FALSE)
  expect_equal(ncol(st_dup$ME), ncol(st_bare$ME))
  expect_equal(unname(st_dup$tones[["hpa"]]),
               unname(st_bare$tones[["hpa"]]), tolerance = 1e-10)
})
