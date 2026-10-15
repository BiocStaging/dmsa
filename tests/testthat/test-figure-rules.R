## Spec 29-38 (figure rules) + the gene-results summary section, added after a
## real Alpha battery run put ~70 genes from 2 systems in ONE uncolored panel
## 6024 px tall, and summary.md never stated a single gene number (PI,
## 2026-08-29). .rp_fig_overview() is tested directly on a fabricated results
## table: one panel per system, at most 40 units per page, per-system accent
## colors, near-threshold shading, and the files it claims to write.

.fr_fake <- function() list(plot_type = "png", palette = "viridis",
                            alpha = 0.05, correction = "maxT", labels = NULL)

.res_fake <- function(n1 = 52, n2 = 12) {
  n <- n1 + n2
  data.frame(
    outcome = "y", level = "gene", n_probes = 1L,
    system_id = rep(c("s1", "s2"), c(n1, n2)),
    system = rep(c("HPA axis & glucocorticoid signalling",
                   "Oxytocin, vasopressin & CICR"), c(n1, n2)),
    unit = sprintf("GENE%02d", seq_len(n)),
    p_unit_adj = c(c(.05, .06, .30, .40), runif(n1 - 4, .3, 1),
                   runif(n2, .3, 1)),
    p_omnibus = runif(n, .01, 1),
    ## naming/ordering/shading are keyed to the BEST LENS adjusted p (the
    ## any-lens rule): the first two genes are named, the next two shaded
    p_coherence_adj = c(c(.02, .03, .10, .19), runif(n1 - 4, .25, 1),
                        runif(n2, .25, 1)),
    p_composite_adj = runif(n, .25, 1),
    p_diffuse_adj = runif(n, .25, 1),
    selected = c(rep(TRUE, 2), rep(FALSE, n - 2)),
    stringsAsFactors = FALSE)
}

test_that("overview: one panel per system, 40-unit pages, named files", {
  skip_if_not_installed("png")   ## the page-height assertions read the PNG back
  set.seed(7)
  fr <- .fr_fake(); res <- .res_fake(52, 12)
  base <- file.path(tempfile("figrules"), "overview_y")
  dir.create(dirname(base), recursive = TRUE)
  files <- .rp_fig_overview(fr, res, "y", base)

  ## 52 genes -> 2 pages for system 1; 12 -> 1 page for system 2
  expect_length(files, 3L)
  expect_true(all(file.exists(files)))
  expect_equal(sum(grepl("HPA_axis", files)), 2L)
  expect_equal(sum(grepl("_p1\\.png$", files)), 1L)
  expect_equal(sum(grepl("_p2\\.png$", files)), 1L)
  expect_equal(sum(grepl("Oxytocin", files)), 1L)
  ## no page may exceed 40 rows: at 0.24 in/row + 1.9 in margins and 300 dpi
  ## the tallest legal page is (40*.24+1.9)*300 = 3450 px
  for (f in files) {
    d <- dim(png::readPNG(f))
    expect_lte(d[1], 3450L)
  }
})

test_that("overview: system slug and page suffix appear only when needed", {
  set.seed(7)
  fr <- .fr_fake()
  res <- .res_fake(15, 12)
  res <- res[res$system_id == "s1", , drop = FALSE]   # ONE system, ONE page
  base <- file.path(tempfile("figrules"), "overview_y")
  dir.create(dirname(base), recursive = TRUE)
  files <- .rp_fig_overview(fr, res, "y", base)
  expect_length(files, 1L)
  expect_identical(basename(files), "overview_y.png")
})

test_that("overview: pages are split by naming-statistic rank within a system", {
  skip_if_not_installed("png")   ## the page-height assertions read the PNG back
  set.seed(7)
  fr <- .fr_fake(); res <- .res_fake(52, 12)
  base <- file.path(tempfile("figrules"), "overview_y")
  dir.create(dirname(base), recursive = TRUE)
  files <- .rp_fig_overview(fr, res, "y", base)
  ## page 1 of system 1 holds its 40 smallest naming statistics (best-lens
  ## adjusted p): it must be the taller file of the two pages (40 rows vs 12)
  p1 <- files[grepl("_p1\\.png$", files)]; p2 <- files[grepl("_p2\\.png$", files)]
  expect_gt(dim(png::readPNG(p1))[1], dim(png::readPNG(p2))[1])
})

test_that("overview: an empty result set writes nothing and returns NULL", {
  fr <- .fr_fake(); res <- .res_fake(4, 2)
  res$n_probes <- 0L
  base <- file.path(tempfile("figrules"), "overview_y")
  dir.create(dirname(base), recursive = TRUE)
  expect_null(.rp_fig_overview(fr, res, "y", base))
  expect_length(list.files(dirname(base)), 0L)
})

## ---- gene-significant / probe-silent locus panels (PI, 2026-08-29) --------
## Gene-level significance can only come from the probes, pooled; when no
## probe clears alpha on its own the locus panel must still be drawn, with
## the probes greyed and the pooling stated - not left for the reader to
## reconcile "significant gene" with "no significant probe".

test_that(".rp_probe_fits returns the nominal per-probe p of its own fit", {
  set.seed(11)
  n <- 60
  d <- data.frame(y = rnorm(n), cov1 = rnorm(n))
  M <- cbind(cgA = d$y * .5 + rnorm(n, 0, .5), cgB = rnorm(n))
  fr <- list(data = d, M = M, outcome = "y", covariates = "cov1", chip = "")
  f <- .rp_probe_fits(fr, "y", c("cgA", "cgB"))
  ## against lm() directly
  m <- lm(M[, "cgA"] ~ y + cov1, data = d)
  expect_equal(unname(f$p["cgA"]), summary(m)$coefficients["y", 4],
               tolerance = 1e-10)
  expect_lt(f$p[["cgA"]], .05)
  expect_gt(f$p[["cgB"]], .05)
})

test_that("a locus panel greys probes whose nominal p is above alpha", {
  set.seed(11)
  pr <- data.frame(probe = paste0("cg", 1:4), b = c(.02, -.01, .015, -.02),
                   se = rep(.02, 4), d = c(1, -1, 1, -1),
                   p = c(.3, .6, .2, .4))
  out <- suppressMessages(
    dmsa_plot_locus(pr, gene = "FAKE1", signal_p = .05,
                    file = tempfile(fileext = ".png")))
  expect_false(any(out$sig))          # all grey
  pr$p[2] <- .01
  out2 <- suppressMessages(
    dmsa_plot_locus(pr, gene = "FAKE1", signal_p = .05,
                    file = tempfile(fileext = ".png")))
  expect_identical(out2$sig, c(FALSE, TRUE, FALSE, FALSE))
})

test_that("the report's locus panel is drawn even when every probe is grey", {
  set.seed(11)
  n <- 60
  d <- data.frame(y = rnorm(n), cov1 = rnorm(n))
  M <- vapply(1:3, function(i) rnorm(n), numeric(n))   # pure noise probes
  colnames(M) <- paste0("cg000fake", 1:3)
  fr <- list(data = d, M = M, outcome = "y", covariates = "cov1", chip = "",
             alpha = .05, labels = NULL, plot_type = "png",
             sets_source = NULL, gene_models = NULL,
             sets = list(s1 = list(map = data.frame(
               gene = "FAKE1", probe = colnames(M), column = colnames(M),
               best_direction = c(1, -1, 1), stringsAsFactors = FALSE))))
  fp <- file.path(tempfile("locusgrey"), "locus_FAKE1_y")
  dir.create(dirname(fp), recursive = TRUE)
  r <- suppressMessages(.rp_fig_locus(fr, "y", "s1", "FAKE1", TRUE, fp))
  expect_identical(as.character(r), fp)
  expect_true(file.exists(paste0(fp, ".png")))
  ## spec 52: with no engine fit handed over, the panel says it fell back
  rows <- attr(r, "probe_rows")
  expect_equal(nrow(rows), 3L)
  expect_true(all(rows$estimator == "ols_fallback"))
})

## ---- a NAMED gene always gets its locus panel (PI, 2026-08-29) -------------
## "if a gene is named, its plot must be printed": the report draws a locus
## panel for every any-lens-named gene, falling back to the bare panel when
## the full one fails, and summary.md points to each panel by file name.

test_that("every named gene's locus panel is written and named in summary.md", {
  set.seed(21)
  pairs <- data.frame(
    cpg_id = sprintf("cg%07d", 1:3),
    target_gene = c("NR3C1", "FKBP5", "CRH"),
    best_direction = c(-1, 1, 1),
    probability_plus1 = c(.05, .9, .9),
    usable = TRUE, abstain_reason = NA_character_,
    best_evidence = "smr_high", direction_tier = "S1",
    stringsAsFactors = FALSE)
  n <- 80
  d <- data.frame(y = rnorm(n), cov1 = rnorm(n), cID = rep(1:40, each = 2))
  for (cl in pairs$cpg_id) d[[cl]] <- stats::plogis(rnorm(n))
  ## plant a strong aligned signal on NR3C1's probe so it is NAMED
  d$cg0000001 <- stats::plogis(-1.4 * scale(d$y)[, 1] + rnorm(n, 0, .4))
  od <- tempfile("dmsa_locus_named")
  old <- options(dmsa.pair_table = pairs); on.exit(options(old))
  f <- suppressMessages(dmsa_frame(
    d, methylation = pairs$cpg_id, direction_source = "cpgdirection",
    outcome = "y", covariates = "cov1", blocks = "cID",
    chip = FALSE, B = 99, progress = FALSE, beep = FALSE, outdir = od))
  r <- suppressMessages(dmsa_report(f))
  hits <- r$results[r$results$selected & r$results$level == "gene", ]
  skip_if(nrow(hits) == 0, "signal did not survive at B = 99 on this seed")
  for (i in seq_len(nrow(hits)))
    expect_true(file.exists(file.path(od, "figures",
      sprintf("locus_%s_%s.png", hits$unit[i], hits$outcome[i]))),
      info = hits$unit[i])
  sm <- readLines(file.path(od, "summary.md"))
  expect_true(any(grepl("Each named gene has its locus panel", sm)))
  expect_true(any(grepl(sprintf("locus_%s_", hits$unit[1]), sm)))
  ## the headline uses the naming rule's wording (re-ruled 2026-10-15:
  ## calibrated FDR), with the joint union test and the pre-registered
  ## any-lens rule reported beside it
  expect_true(any(grepl("named at FDR", sm)))
  expect_true(any(grepl("joint union adj", sm)))
  expect_true(any(grepl("Pre-registered any-lens rule", sm)))
  ## the badge tier (PI ruling 2026-08-29): exact_confirmed rides the units
  ## table, implies selected, and the run's measured FWER is stated
  u <- r$results
  expect_true(all(c("p_union_exact", "exact_confirmed", "p_omnibus_adj",
                    "omnibus_confirmed", "fwer_realized") %in% names(u)))
  expect_true(all(!u$exact_confirmed | u$selected))
  expect_true(all(!u$omnibus_confirmed | u$selected))
  expect_true(all(u$fwer_realized[is.finite(u$fwer_realized)] >= 0 &
                  u$fwer_realized[is.finite(u$fwer_realized)] <= 1))
  expect_true(any(grepl("realized family-wise error of the pre-registered any-lens rule", sm)))
})

## ---- moderation: runs WITH the main analysis; binary-aware figure ---------
## (PI, 2026-08-29: a moderation run produced no gene/module/system results
## at all, drew linear-probability scatters for a 0/1 outcome, and a
## Johnson-Neyman continuum over a two-level sex moderator.)

test_that("moderation = TRUE no longer replaces the main analysis", {
  set.seed(31)
  pairs <- data.frame(cpg_id = sprintf("cg%07d", 1:3),
    target_gene = c("NR3C1", "FKBP5", "CRH"), best_direction = c(-1, 1, 1),
    probability_plus1 = c(.05, .9, .9), usable = TRUE,
    abstain_reason = NA_character_, best_evidence = "smr_high",
    direction_tier = "S1", stringsAsFactors = FALSE)
  n <- 120
  d <- data.frame(pills = rep(c(0, 1), each = n / 2),
                  sex = rep(c(1, 2), n / 2), cov1 = rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in pairs$cpg_id) d[[cl]] <- stats::plogis(rnorm(n))
  d$cg0000001 <- stats::plogis(-.8 * scale(d$pills)[, 1] * (d$sex - 1.5) +
                               rnorm(n, 0, .6))
  od <- tempfile("dmsa_modmain")
  old <- options(dmsa.pair_table = pairs); on.exit(options(old))
  f <- suppressMessages(dmsa_frame(d, methylation = pairs$cpg_id,
    direction_source = "cpgdirection", predictors = "pills",
    covariates = "cov1", blocks = "cID", chip = FALSE,
    moderation = TRUE, mod = "sex",
    mod_levels = c("1" = "Husband", "2" = "Wife"),
    B = 49, progress = FALSE, beep = FALSE, outdir = od))
  r <- suppressWarnings(dmsa_report(f))
  ## BOTH batteries in one report
  expect_true(sum(r$results$level == "gene", na.rm = TRUE) > 0)
  expect_true(nrow(r$moderation) > 0)
  ## the main overview figure exists alongside any moderation figures
  expect_true(any(grepl("^overview_", list.files(file.path(od, "figures")))))
})

test_that("a binary outcome under frame_role='outcome' is not refused", {
  ## the gaussian gate applies only when the outcome is the RESPONSE of the
  ## moderated model (frame_role = 'predictor'); with predictors= the tone
  ## score is the response and a two-level oc is a plain group contrast
  set.seed(31)
  map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
                    probe = "cg01", column = "cg01",
                    best_direction = -1, p_plus = .1)
  d <- data.frame(pills = rep(c(0, 1), 30), sex = rep(c(1, 2), each = 30),
                  cov1 = rnorm(60), cID = rep(1:30, each = 2),
                  cg01 = stats::plogis(rnorm(60)))
  f <- suppressMessages(dmsa_frame(d, map = map, predictors = "pills",
    covariates = "cov1", blocks = "cID", moderation = TRUE,
    mod = "sex", B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
    progress = FALSE, beep = FALSE, outdir = tempfile()))
  expect_s3_class(suppressWarnings(dmsa_report(f)), "dmsa_report")
})
