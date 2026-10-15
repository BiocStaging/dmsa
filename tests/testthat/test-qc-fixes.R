## Regression pins for the 2026-08-28 whole-package QC pass.

.qc_fixture <- function(n = 60, k = 3, genes = "FKBP5") {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene %in% genes], k)
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  list(data = d, cols = cols, n = n)
}

.qc_run <- function(fx, data = fx$data, methylation = fx$cols, ...)
  dmsa_frame(data, methylation = methylation, direction_source = "bundled",
             outcome = "y",
             covariates = "cov1", blocks = "cID", chip = FALSE,
             B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
             progress = FALSE, beep = FALSE, outdir = tempfile("dmsa_qc"), ...)

## ---- sample alignment ------------------------------------------------------

test_that("a NUMERIC id column is auto-detected from matrix rownames", {
  set.seed(1); fx <- .qc_fixture()
  d <- fx$data; d$pid <- 1000 + seq_len(fx$n)
  M <- as.matrix(d[fx$cols]); rownames(M) <- as.character(d$pid)
  sh <- sample(fx$n); Msh <- M[sh, , drop = FALSE]
  expect_message(f <- .qc_run(fx, data = d, methylation = Msh),
                 "auto-detected")
  expect_true(any(grepl("different order", f$corrections$issue)))
})

test_that("rownames with no matching data column are an error, not positional", {
  set.seed(1); fx <- .qc_fixture()
  M <- as.matrix(fx$data[fx$cols])
  rownames(M) <- paste0("S", seq_len(fx$n))   # matches nothing in `data`
  e <- expect_error(.qc_run(fx, methylation = M))
  expect_match(conditionMessage(e), "no column of `data` matches")
  expect_match(conditionMessage(e), "No analysis was run")
})

## ---- logical outcome -------------------------------------------------------

test_that("a logical outcome runs as 0/1 instead of dissolving into NA", {
  set.seed(1); fx <- .qc_fixture()
  d <- fx$data; d$hit <- d$y > 0
  f <- dmsa_frame(d, methylation = fx$cols, direction_source = "bundled",
                  outcome = "hit",
                  covariates = "cov1", blocks = "cID", chip = FALSE,
                  B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
                  progress = FALSE, beep = FALSE, outdir = tempfile("qc"))
  expect_s3_class(f, "dmsa_frame")
  expect_identical(f$outcome_type, "logistic")   # two-level autofix still fires
  expect_true(any(grepl("TRUE/FALSE", unlist(f$corrections))))
})

## ---- mixed matrix column naming -------------------------------------------

test_that("a matrix mixing column names and bare probe ids loses nothing", {
  set.seed(1); fx <- .qc_fixture()
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  bare <- mp$probe[match(fx$cols, mp$column)]
  M <- as.matrix(fx$data[fx$cols])
  colnames(M) <- c(fx$cols[1], bare[2], fx$cols[3])   # mixed naming
  f <- .qc_run(fx, methylation = M)
  expect_s3_class(f, "dmsa_frame")
  expect_equal(nrow(f$map), 3L)                        # nothing became all-NA
})

## ---- outcome_type positional refusal survives the default sentinel --------

test_that("a positional 3-vector equal to the default is still refused", {
  set.seed(1); fx <- .qc_fixture()
  d <- fx$data; d$y2 <- stats::rnorm(fx$n); d$y3 <- stats::rnorm(fx$n)
  e <- expect_error(.qc_run(fx, data = d,
                            outcome_type = c("gaussian", "logistic",
                                             "multinomial")))
  expect_match(conditionMessage(e), "does not name them")
})

## ---- cascade: renamed alpha keeps polarity ---------------------------------

test_that("renaming the bundled alpha cascade does not drop its polarity", {
  s1 <- dmsa_sets("alpha")
  s2 <- dmsa_sets("alpha", name = "my renamed run")
  expect_false(is.null(s1$polarity))
  expect_false(is.null(s2$polarity))
  expect_equal(nrow(s2$polarity), nrow(s1$polarity))
})

## ---- id columns never parsed as numeric ------------------------------------

test_that("module ids like 2.10 survive a reference round-trip", {
  ref <- dmsa_reference(
    systems = data.frame(gene = paste0("G", 1:6), system_id = "1",
                         system = "S",
                         module_id = c("2.1", "2.1", "2.1", "2.10", "2.10",
                                       "2.10"),
                         module = c(rep("m one", 3), rep("m ten", 3))),
    polarity = data.frame(gene = paste0("G", 1:6), system_id = "1",
                          w_g = 1))
  d <- tempfile("refdir"); dmsa_reference_write(ref, d)
  back <- dmsa_reference_read(d)
  expect_setequal(unique(back$systems$module_id), c("2.1", "2.10"))
})

## ---- alpha_reference(modules=) is a per-gene override ----------------------

test_that("a partial modules table re-partitions only the genes it names", {
  ref0 <- alpha_reference()
  g <- ref0$systems$gene[!is.na(ref0$systems$module_id)][1]
  sid <- ref0$systems$system_id[ref0$systems$gene == g][1]
  ref1 <- alpha_reference(modules = data.frame(
    gene = g, system_id = sid, module_id = "99.1", module = "custom"))
  s1 <- ref1$systems
  expect_identical(s1$module_id[s1$gene == g & s1$system_id == sid], "99.1")
  ## every OTHER gene keeps its codebook module
  keep0 <- ref0$systems[!(ref0$systems$gene == g &
                          ref0$systems$system_id == sid), ]
  keep1 <- s1[!(s1$gene == g & s1$system_id == sid), ]
  expect_equal(sum(is.na(keep1$module_id)), sum(is.na(keep0$module_id)))
})

## ---- reference_csv: role matching is exact ---------------------------------

test_that("brake-of-driver never becomes an anchor", {
  f <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    system = "S1", gene = paste0("G", 1:4),
    w_g = c(1, -1, 1, 1),
    role = c("driver", "brake-of-driver", "support", "support")),
    f, row.names = FALSE)
  ref <- dmsa_reference_csv(f, quiet = TRUE, min_genes = 3L)
  expect_identical(ref$anchors$gene, "G1")
})

## ---- report: overwrite guard covers a summary-only report ------------------

test_that("a summary-only report is protected from overwrite too", {
  set.seed(1); fx <- .qc_fixture()
  od <- tempfile("qc_ow")
  f <- dmsa_frame(fx$data, methylation = fx$cols, direction_source = "bundled",
                  outcome = "y",
                  covariates = "cov1", blocks = "cID", chip = FALSE,
                  B = 19, plots = FALSE, tables = FALSE, summary = TRUE,
                  progress = FALSE, beep = FALSE, outdir = od)
  r1 <- dmsa_report(f)
  expect_true(file.exists(file.path(od, "summary.md")))
  expect_error(dmsa_report(f), "already contains")
  expect_s3_class(dmsa_report(f, overwrite = TRUE), "dmsa_report")
})

## ---- report: spec 27 per-lens flags ----------------------------------------

test_that("the units table carries per-lens survivor flags and p_unit", {
  set.seed(1); fx <- .qc_fixture()
  f <- dmsa_frame(fx$data, methylation = fx$cols, direction_source = "bundled",
                  outcome = "y",
                  covariates = "cov1", blocks = "cID", chip = FALSE,
                  B = 19, plots = FALSE, tables = TRUE, summary = FALSE,
                  progress = FALSE, beep = FALSE, table_type = "html",
                  outdir = tempfile("qc27"))
  r <- dmsa_report(f)
  u <- r$results
  expect_true(all(c("selected_coherence", "selected_composite",
                    "selected_diffuse", "n_lenses_hit", "any_lens_hit",
                    "p_unit_adj") %in% names(u)))
  ## RE-RULED 2026-10-15 (supersedes E2): naming is the dependence-calibrated
  ## BH q on p_unit at fdr; the pre-registered any-lens rule is kept beside
  ## it as prereg_any_lens (= any_lens_hit), with its realized FWER disclosed.
  expect_identical(u$selected, is.finite(u$q_unit_dep) &
                     u$q_unit_dep <= f$fdr & u$n_probes > 0)
  expect_identical(u$prereg_any_lens, u$any_lens_hit)
  expect_true(all(is.finite(u$p_unit_adj[u$n_probes > 0])))
  expect_true(all(u$n_lenses_hit ==
                  u$selected_coherence + u$selected_composite +
                  u$selected_diffuse))
})
