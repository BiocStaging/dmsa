## Spec 62 (continued): the display RULES the PI ruled on, pinned.
##
## test-figure-rules.R pins the overview mechanics and the locus panel. These
## pin what the moderation figures actually SAY and how they are laid out:
##   * a two-level moderator draws ONE panel, each level's own slope and p in
##     the legend, no Johnson-Neyman continuum ("a significant interaction
##     means, by definition, that the slopes differ")
##   * a two-level outcome is displayed as fitted probabilities (logistic
##     display; the tested statistic stays the linear product)
##   * a continuous moderator keeps the quantile ladder and the JN band
##   * the overview carries only its title
##   * the confirmation badges reach the summary beside the named gene
##
## Figures are written as PDF with compression off, so their text can be
## read back: R's pdf device writes each string as a (..) Tj or a kerned
## [(..) n (..)] TJ operator. No pixels, no OCR.

.fd_pdf_text <- function(file) {
  l <- readLines(file, warn = FALSE)
  l <- l[grepl("T[jJ]\\s*$", l)]
  vapply(l, function(x) {
    m <- regmatches(x, gregexpr("\\((?:[^()\\\\]|\\\\.)*\\)", x, perl = TRUE))[[1]]
    gsub("\\\\([()\\\\])", "\\1", paste(substr(m, 2, nchar(m) - 1), collapse = ""))
  }, character(1), USE.NAMES = FALSE)
}
.fd_pdf_width_in <- function(file) {
  l <- grep("MediaBox", readLines(file, warn = FALSE), value = TRUE)[1]
  as.numeric(strsplit(sub(".*MediaBox \\[([^]]*)\\].*", "\\1", l), " +")[[1]])[3] / 72
}
.fd_pairs <- function() data.frame(
  cpg_id = sprintf("cg%07d", 1:3), target_gene = c("NR3C1", "FKBP5", "CRH"),
  best_direction = c(-1, 1, 1), probability_plus1 = c(.05, .9, .9),
  usable = TRUE, abstain_reason = NA_character_, best_evidence = "smr_high",
  direction_tier = "S1", stringsAsFactors = FALSE)
.fd_frame <- function(d, od, ...) {
  pairs <- .fd_pairs()
  old <- options(dmsa.pair_table = pairs); on.exit(options(old))
  oldpdf <- grDevices::pdf.options(); grDevices::pdf.options(compress = FALSE)
  on.exit(do.call(grDevices::pdf.options, oldpdf), add = TRUE)
  f <- suppressMessages(dmsa_frame(d, methylation = pairs$cpg_id,
    direction_source = "cpgdirection", covariates = "cov1",
    blocks = "cID", chip = FALSE, B = 49, progress = FALSE,
    beep = FALSE, outdir = od, plot_type = "pdf", ...))
  suppressWarnings(suppressMessages(dmsa_report(f)))
  f
}

test_that("spec 62: two-level moderator x two-level outcome = ONE panel, per-level b/p, probabilities, no JN", {
  set.seed(31); n <- 120
  d <- data.frame(pills = rep(c(0, 1), each = n / 2), sex = rep(c(1, 2), n / 2),
                  cov1 = stats::rnorm(n), cID = rep(seq_len(n / 2), each = 2))
  for (cl in .fd_pairs()$cpg_id) d[[cl]] <- stats::plogis(stats::rnorm(n))
  d$cg0000001 <- stats::plogis(-.8 * scale(d$pills)[, 1] * (d$sex - 1.5) +
                               stats::rnorm(n, 0, .6))
  od <- tempfile("dmsa_fd_bin")
  .fd_frame(d, od, predictors = "pills", moderation = TRUE, mod = "sex",
            mod_levels = c("1" = "Husband", "2" = "Wife"))
  fig <- file.path(od, "figures", "moderation_NR3C1_pills.pdf")
  expect_true(file.exists(fig))
  ## one panel: the narrow device, not the two-panel one
  expect_lt(.fd_pdf_width_in(fig), 8)
  tx <- .fd_pdf_text(fig)
  ## each level's own slope and p, on the log-odds scale, in the legend
  expect_true(any(grepl("^Husband \\[1\\]:  b [-+][0-9.]+, p .*\\(log-odds\\)$", tx)))
  expect_true(any(grepl("^Wife \\[2\\]:  b [-+][0-9.]+, p .*\\(log-odds\\)$", tx)))
  ## the display is a probability of the outcome
  expect_true(any(grepl("^P\\(pills = 1\\)$", tx)))
  expect_true(any(grepl("Fitted probability of pills = 1 per level of sex", tx)))
  ## the caption states the rule
  expect_true(any(grepl("The tested interaction IS the difference between the two slopes", tx)))
  ## and nothing from the continuum panel
  expect_false(any(grepl("Johnson|Shaded: slope differs|pct\\)", tx)))
  expect_true(any(grepl("^NR3C1 x sex -> pills", tx)))
})

test_that("spec 62: two-level moderator x continuous outcome = one panel, linear slopes, no log-odds", {
  set.seed(33); n <- 120
  d <- data.frame(y = stats::rnorm(n), sex = rep(c(1, 2), n / 2),
                  cov1 = stats::rnorm(n), cID = rep(seq_len(n / 2), each = 2))
  for (cl in .fd_pairs()$cpg_id) d[[cl]] <- stats::plogis(stats::rnorm(n))
  d$cg0000001 <- stats::plogis(-.8 * scale(d$y)[, 1] * (d$sex - 1.5) +
                               stats::rnorm(n, 0, .6))
  od <- tempfile("dmsa_fd_lin")
  .fd_frame(d, od, predictors = "y", moderation = TRUE, mod = "sex",
            mod_levels = c("1" = "Husband", "2" = "Wife"))
  fig <- file.path(od, "figures", "moderation_NR3C1_y.pdf")
  expect_true(file.exists(fig))
  expect_lt(.fd_pdf_width_in(fig), 8)
  tx <- .fd_pdf_text(fig)
  expect_true(any(grepl("^Husband \\[1\\]:  b [-+][0-9.]+, p ", tx)))
  expect_true(any(grepl("^Wife \\[2\\]:  b [-+][0-9.]+, p ", tx)))
  expect_false(any(grepl("log-odds|^P\\(", tx)))
  expect_false(any(grepl("Johnson|Shaded: slope differs|pct\\)", tx)))
})

test_that("spec 62: a continuous moderator keeps the quantile ladder and the Johnson-Neyman panel", {
  set.seed(32); n <- 120
  d <- data.frame(y = stats::rnorm(n), age = stats::rnorm(n),
                  cov1 = stats::rnorm(n), cID = rep(seq_len(n / 2), each = 2))
  for (cl in .fd_pairs()$cpg_id) d[[cl]] <- stats::plogis(stats::rnorm(n))
  d$cg0000001 <- stats::plogis(-.8 * scale(d$y)[, 1] * d$age + stats::rnorm(n, 0, .6))
  od <- tempfile("dmsa_fd_cont")
  .fd_frame(d, od, predictors = "y", moderation = TRUE, mod = "age")
  fig <- file.path(od, "figures", "moderation_NR3C1_y.pdf")
  expect_true(file.exists(fig))
  ## two panels: the wide device
  expect_gt(.fd_pdf_width_in(fig), 9)
  tx <- .fd_pdf_text(fig)
  expect_true(any(grepl("^age = .* \\(10th pct\\)$", tx)))
  expect_true(any(grepl("^age = .* \\(50th pct\\)$", tx)))
  expect_true(any(grepl("^age = .* \\(90th pct\\)$", tx)))
  expect_true(any(grepl("^Simple slopes at 3 levels of age", tx)))
  expect_true(any(grepl("^slope of the tone score$", tx)))
  expect_true(any(grepl("^Shaded: slope differs from zero", tx)))
  expect_true(any(grepl("^age \\(observed range\\)$", tx)))
  expect_false(any(grepl("log-odds|The tested interaction IS the difference", tx)))
})

test_that("spec 62: the overview carries only its title", {
  set.seed(31); n <- 120
  d <- data.frame(pills = rep(c(0, 1), each = n / 2), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in .fd_pairs()$cpg_id) d[[cl]] <- stats::plogis(stats::rnorm(n))
  od <- tempfile("dmsa_fd_ov")
  .fd_frame(d, od, predictors = "pills")
  ov <- list.files(file.path(od, "figures"), pattern = "^overview_.*\\.pdf$",
                   full.names = TRUE)
  expect_length(ov, 1L)
  tx <- .fd_pdf_text(ov)
  ## title, axis, lens legend and unit labels - and no explanatory subtitle
  expect_true(any(grepl("^pills - HPA axis", tx)))
  ## axis and legend re-ruled 2026-10-15: the calibrated q (the naming
  ## statistic) is drawn beside the three per-lens adjusted p's
  expect_true(any(grepl("^calibrated q / family-adjusted p$", tx)))
  expect_true(all(c("calibrated q", "coherence", "composite", "diffuse",
                    "NR3C1", "FKBP5", "CRH") %in% tx))
  expect_false(any(grepl("Bold|named by|shaded|dagger|exact-confirmed", tx, ignore.case = TRUE)))
})
