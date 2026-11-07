## E17 (PI ruling, 2026-11-07): when a direction table carries
## cpgdirection's best_accuracy, P(d = +1) is best_accuracy for a +1 call and
## 1 - best_accuracy for a -1 call, ahead of every other column. The old
## ladder (p_plus -> probability_plus1 -> best_confidence -> smr_tier) serves
## only tables without that column.

.ba_base <- function() data.frame(
  cpg_id = paste0("cg", 1:6),
  best_direction = c(1, -1, 1, -1, 1, -1),
  stringsAsFactors = FALSE)

test_that("best_accuracy sets P(d = +1) directly, by the sign of the call", {
  d <- .ba_base()
  d$best_accuracy <- c(.97, .97, .67, .67, .5, .5)
  al <- dmsa_align(d, genes = rep("OXTR", 6))
  expect_equal(al$p_plus, c(.97, .03, .67, .33, .5, .5))
  expect_true(all(al$usable))
})

test_that("best_accuracy outranks every other confidence column", {
  ## a measured call reads best_confidence = 1 whatever its validated
  ## accuracy; an SMR row borrows a tier midpoint; neither may win once the
  ## validated figure is on the row
  d <- .ba_base()
  d$p_plus <- rep(.99, 6)
  d$probability_plus1 <- rep(.01, 6)
  d$best_confidence <- rep(1, 6)
  d$smr_tier <- rep("S1", 6)
  d$best_accuracy <- c(.97, .97, .80, .80, .62, .62)
  al <- dmsa_align(d, genes = rep("OXTR", 6))
  expect_equal(al$p_plus, c(.97, .03, .80, .20, .62, .38))
})

test_that("without best_accuracy the old ladder is unchanged", {
  d <- .ba_base()
  d$best_confidence <- c(.6, .6, NA, NA, 1, 1)
  d$smr_tier <- c(NA, NA, "S2", "S3", NA, NA)
  al <- dmsa_align(d, genes = rep("OXTR", 6))
  expect_equal(round(al$p_plus, 3), c(.8, .2, .85, .295, 1, 0))
})

test_that("a direction with no best_accuracy value is unusable, not filled from the ladder", {
  d <- .ba_base()
  d$best_confidence <- rep(1, 6)
  d$smr_tier <- rep("S1", 6)
  d$best_accuracy <- c(.9, .9, NA, NA, .7, .7)
  expect_warning(al <- dmsa_align(d, genes = rep("OXTR", 6)),
                 "no confidence value")
  expect_equal(sum(al$usable), 4L)
  expect_true(all(is.na(al$p_plus[3:4])))
  expect_equal(al$reason[3:4], c("no_confidence", "no_confidence"))
})

test_that("an abstained row stays abstained and does not warn", {
  d <- .ba_base()
  d$best_direction[2] <- NA
  d$best_accuracy <- c(.9, NA, .9, .9, .9, .9)
  expect_no_warning(al <- dmsa_align(d, genes = rep("OXTR", 6)))
  expect_true(is.na(al$p_plus[2]))
  expect_false(al$usable[2])
  expect_equal(sum(al$usable), 5L)
})

test_that("best_accuracy outside [0, 1] is clamped", {
  d <- .ba_base()[1:2, ]
  d$best_accuracy <- c(1.2, -0.1)
  al <- dmsa_align(d, genes = rep("OXTR", 2))
  expect_equal(al$p_plus, c(1, 1))
})

## the pair runtime hands the raw pair columns to dmsa_align() once, so the
## same rule must reach the engine's map
.ba_pairs <- function() data.frame(
  cpg_id         = c("cg0000001", "cg0000001", "cg0000002", "cg0000003",
                     "cg0000004"),
  target_gene    = c("OXTR", "AVP", "NR3C1", "CRH", "FKBP5"),
  best_direction = c(-1, 1, -1, NA, 1),
  probability_plus1 = c(0.05, 0.9, 0.08, NA, 0.8),
  best_confidence = c(1, NA, NA, NA, .6),
  smr_tier       = c(NA, "S3", "S1", NA, NA),
  best_accuracy  = c(.97, .71, .95, NA, .67),
  usable         = c(TRUE, TRUE, TRUE, FALSE, TRUE),
  abstain_reason = c(NA, NA, NA, "tissue_conflict", NA),
  best_evidence  = c("measured", "smr_weak", "smr_high", "tissue_conflict",
                     "catalogue_single"),
  direction_tier = c("M", "S3", "S1", "U", "B"),
  stringsAsFactors = FALSE)

test_that("the pair runtime takes p_plus from best_accuracy", {
  cols <- c("cg0000001", "cg0000002", "cg0000003", "cg0000004")
  rt <- .frame_pairs_runtime(cols, alpha_reference(), pairs = .ba_pairs())
  mp <- rt$full
  key <- paste(mp$probe, mp$gene)
  got <- setNames(mp$p_plus, key)
  expect_equal(unname(got[c("cg0000001 OXTR", "cg0000001 AVP",
                            "cg0000002 NR3C1", "cg0000004 FKBP5")]),
               c(.03, .71, .05, .67))
  expect_true("best_accuracy" %in% names(mp))
  expect_false("CRH" %in% mp$gene)
})

test_that("a usable pair without a best_accuracy value is reported, not weighted", {
  cols <- c("cg0000001", "cg0000002", "cg0000004")
  pt <- .ba_pairs()
  pt$best_accuracy[5] <- NA
  rt <- .frame_pairs_runtime(cols, alpha_reference(), pairs = pt)
  expect_false("FKBP5" %in% rt$full$gene)
  led <- rt$ledger
  r <- led[led$cpg_id == "cg0000004" & led$target_gene == "FKBP5", ]
  expect_false(r$used)
  expect_match(r$reason_not_used, "no recoverable confidence")
})

test_that("a pair table without best_accuracy keeps the old weights", {
  cols <- c("cg0000001", "cg0000002", "cg0000004")
  pt <- .ba_pairs()
  pt$best_accuracy <- NULL
  rt <- .frame_pairs_runtime(cols, alpha_reference(), pairs = pt)
  got <- setNames(rt$full$p_plus, paste(rt$full$probe, rt$full$gene))
  expect_equal(unname(got[c("cg0000001 OXTR", "cg0000001 AVP",
                            "cg0000002 NR3C1", "cg0000004 FKBP5")]),
               c(.05, .9, .08, .8))
})

test_that("a user map carrying best_accuracy sets p_plus in .frame_maps()", {
  mp <- data.frame(
    probe = paste0("cg", 1:4), column = paste0("cg", 1:4),
    gene = c("OXTR", "AVP", "NR3C1", "FKBP5"),
    system_id = "1", system = "Sys A",
    best_direction = c(1, -1, 1, -1),
    p_plus = c(.99, .01, .99, .01),
    best_accuracy = c(.97, .97, .62, NA),
    stringsAsFactors = FALSE)
  m <- .frame_maps(mp)
  expect_equal(nrow(m$full), 3L)
  expect_equal(m$full$p_plus, c(.97, .03, .62))

  ## without the column, p_plus is read as before
  mp$best_accuracy <- NULL
  m2 <- .frame_maps(mp)
  expect_equal(nrow(m2$full), 4L)
  expect_equal(m2$full$p_plus, mp$p_plus)
})

test_that("the bundled maps carry p_plus only, so E17 leaves them unchanged", {
  m <- .frame_read_map("alpha")
  expect_false("best_accuracy" %in% names(m))
})
