# The two-call interface: validation, test drive, and a tiny simulated run.

.sim_frame_data <- function(n = 120, seed = 42) {
  set.seed(seed)
  map <- expand.grid(gene = paste0("G", 1:6), k = 1:4,
                     stringsAsFactors = FALSE)
  map$system_id <- ifelse(map$gene %in% paste0("G", 1:3), 1L, 2L)
  map$system <- ifelse(map$system_id == 1, "Sim system one", "Sim system two")
  map$probe <- sprintf("cg%07d", seq_len(nrow(map)))
  map$column <- paste0(map$probe, "_", map$gene)
  map$best_direction <- rep(c(-1, 1), length.out = nrow(map))
  ## one untestable gene: G6 has no direction calls at all
  map$best_direction[map$gene == "G6"] <- NA
  map$p_plus <- ifelse(map$best_direction > 0, .9, .1)
  map$best_tier <- "A"; map$smr_tier <- ""
  d <- data.frame(out1 = rnorm(n), grp = rep(1:2, length.out = n),
                  cov1 = rnorm(n), cov2 = rnorm(n),
                  m1 = rnorm(n), cID = rep(seq_len(n / 2), each = 2))
  ## beta-scale methylation columns; G1 tracks out1 through its alignment
  for (i in seq_len(nrow(map))) {
    sig <- if (map$gene[i] == "G1")
      0.35 * d$out1 * sign(map$best_direction[i]) else 0
    d[[map$column[i]]] <- stats::plogis(rnorm(n, sd = 1) + sig)
  }
  list(data = d, map = map)
}

test_that("frame validation catches the layout errors", {
  s <- .sim_frame_data()
  expect_error(dmsa_frame(s$data, map = s$map, outcome = "out1",
                          covariates = c("cov1", "cov2"),
                          moderation = TRUE),
               "requires `mod`")
  expect_error(dmsa_frame(s$data, map = s$map, outcome = "out1",
                          covariates = c("cov1", "cov2"),
                          mod2 = "m1"),
               "without `mod`")
  expect_error(dmsa_frame(s$data, map = s$map, outcome = "nope",
                          covariates = "cov1"), "not in `data`")
  expect_error(dmsa_frame(s$data, map = s$map, outcome = "out1",
                          covariates = "cov1", systems = 99),
               "matched nothing")
  expect_error(dmsa_frame(s$data, map = s$map, outcome = "out1",
                          covariates = "cov1", module = TRUE),
               "module")
})

test_that("the test drive corrects loudly and the pilot runs", {
  s <- .sim_frame_data()
  s$data$cov_dup <- s$data$cov1                       # aliased column
  fr <- expect_warning(
    dmsa_frame(s$data, map = s$map, outcome = "out1",
               covariates = c("cov1", "cov2", "cov_dup"),
               blocks = "cID", B = 99, seed = 1),
    regexp = NA)                                       # no warning expected
  expect_s3_class(fr, "dmsa_frame")
  expect_true(any(grepl("aliased", fr$corrections$issue)))
  expect_true(any(grepl("beta", fr$corrections$issue)))
  expect_true(is.list(fr$pilot) && fr$pilot$ok)
  expect_output(print(fr), "DMSA frame")
  expect_output(dmsa_test_drive(fr), "Test drive")
})

test_that("binary outcome under gaussian autofixes to logistic", {
  s <- .sim_frame_data()
  fr <- dmsa_frame(s$data, map = s$map, outcome = "grp",
                   covariates = c("cov1", "cov2"), blocks = "cID",
                   B = 99)
  expect_identical(fr$outcome_type, "logistic")
  ## Assert the CLAIM, not the prose. This used to grep for the literal word
  ## "binary" in the issue text, so rewording the note - which was necessary,
  ## because "switched to logistic" wrongly implied a logistic model is fitted -
  ## failed a test whose behaviour had not changed.
  i <- which(fr$corrections$field == "grp")
  expect_gte(length(i), 1L)
  expect_match(paste(fr$corrections$issue[i], collapse = " "),
               "two-level|binary", ignore.case = TRUE)
  expect_match(paste(fr$corrections$action[i], collapse = " "),
               "contrast", ignore.case = TRUE)
})

test_that("a tiny report runs end to end and finds the planted gene", {
  s <- .sim_frame_data()
  od <- file.path(tempdir(), paste0("dmsa_out_", as.integer(runif(1, 1, 1e6))))
  fr <- dmsa_frame(s$data, map = s$map, outcome = "out1",
                   covariates = c("cov1", "cov2"), blocks = "cID",
                   B = 199, seed = 1, outdir = od, plot_type = "png")
  r <- dmsa_report(fr)
  expect_s3_class(r, "dmsa_report")
  expect_true(file.exists(file.path(od, "summary.md")))
  expect_true(file.exists(file.path(od, "tables", "units.csv")))
  u <- utils::read.csv(file.path(od, "tables", "units.csv"))
  g1 <- u[u$level == "gene" & u$unit == "G1", ]
  expect_true(nrow(g1) == 1 && g1$p_omnibus < .05)
  ## the direction-less gene is tracked as untestable, never analysed silently
  expect_true("G6" %in% fr$untestable$gene)
  expect_false("G6" %in% u$unit)
  smy <- readLines(file.path(od, "summary.md"))
  expect_true(any(grepl("untestable", smy)))
  expect_output(print(r), "DMSA report")
})

test_that("moderation runs composite-only through the new layout", {
  s <- .sim_frame_data()
  od <- file.path(tempdir(), paste0("dmsa_mod_", as.integer(runif(1, 1, 1e6))))
  fr <- dmsa_frame(s$data, map = s$map, outcome = "out1",
                   covariates = c("cov1", "cov2"), blocks = "cID",
                   moderation = TRUE, mod = "m1", B = 99, seed = 1,
                   outdir = od)
  r <- dmsa_report(fr)
  expect_false(is.null(r$moderation))
  expect_true(all(c("p_composite", "p_composite_adj") %in%
                    names(r$moderation)))
  expect_true(all(r$moderation$p_composite_adj >= r$moderation$p_composite -
                    1e-12, na.rm = TRUE))
})

test_that("dmsa_model returns its permutation null when asked", {
  set.seed(9)
  d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60))
  f <- dmsa_model(y ~ x + z, d, "x", B = 99, seed = 1, nulls = TRUE)
  expect_length(f$null_t, 99)
  f0 <- dmsa_model(y ~ x + z, d, "x", B = 99, seed = 1)
  expect_null(f0$null_t)
  expect_equal(f$p_perm, f0$p_perm)
})

# ---------------------------------------------------------------------------
# Progress bar and completion beep. Both are conveniences, so the tests that
# matter are the ones proving they are never load-bearing and never noisy where
# noise is wrong.
# ---------------------------------------------------------------------------

test_that("the work counter tracks the declared design, not a guess", {
  fk <- list(outcome = c("y1", "y2"),
             sets = list(`1` = NULL, `2` = NULL, `3` = NULL),
             modules = data.frame(system_id = c(1, 1, 2), gene = letters[1:3]))
  ## 2 outcomes x (1 system call + 3 gene calls + 2 module systems) = 12
  expect_equal(.rp_work(fk, c("system", "gene", "module")), 12L)
  expect_equal(.rp_work(fk, "gene"), 6L)
  expect_equal(.rp_work(fk, "system"), 2L)
  ## a design with nothing countable still returns a usable maximum
  expect_gte(.rp_work(list(outcome = "y", sets = list(), modules = NULL),
                      character(0)), 1L)
})

test_that("progress = FALSE is a silent no-op, not a hidden bar", {
  b <- .rp_bar(10, on = FALSE)
  expect_silent({ b$tick(); b$tick(3); b$done() })
  ## and an impossible total cannot start one either
  expect_silent({ z <- .rp_bar(0, on = TRUE); z$tick(); z$done() })
  expect_silent({ z <- .rp_bar(NA_real_, on = TRUE); z$tick(); z$done() })
})

test_that("a real bar draws, advances and closes without erroring", {
  out <- utils::capture.output({
    b <- .rp_bar(4, on = TRUE); b$tick(); b$tick(2); b$done()
  })
  expect_true(any(grepl("=|%", paste(out, collapse = ""))))
  ## ticking past the end clamps rather than erroring (a live bar prints, so
  ## capture it - the assertion is that it does not throw)
  expect_error(utils::capture.output({
    b2 <- .rp_bar(2, on = TRUE); b2$tick(99); b2$done() }), NA)
})

test_that("beep never errors, whatever it is handed", {
  ## beepr is a Suggests. Absent it must do nothing; present it must not throw
  ## on any of the accepted forms. Either way this cannot break a report.
  expect_silent(.rp_beep(FALSE))
  expect_silent(.rp_beep(NULL))
  skip_if(requireNamespace("beepr", quietly = TRUE), "beepr would make noise")
  expect_silent(.rp_beep(TRUE))
  expect_silent(.rp_beep(8))
})

test_that("the frame carries progress and beep, and report can override", {
  expect_true("progress" %in% names(formals(dmsa_frame)))
  expect_true("beep" %in% names(formals(dmsa_frame)))
  ## dmsa_report() is a generic since 0.99.17 (it also reports the four
  ## structural methods), so the two live on the frame METHOD
  expect_true(all(c("progress", "beep") %in%
                    names(formals(dmsa_report.dmsa_frame))))
  ## both default to interactive(), so a check run is silent and bar-free
  expect_identical(eval(formals(dmsa_frame)$beep), interactive())
  expect_identical(eval(formals(dmsa_frame)$progress), interactive())
})
