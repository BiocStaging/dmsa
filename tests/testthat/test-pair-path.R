## Spec 1/3/13-17/59: the CpG x gene PAIR path through dmsa_frame().
## Runs on every machine: the pair table is injected through
## options(dmsa.pair_table = ...), the documented seam, so no test here
## depends on cpgdirection being installed.

.pp_ref <- function() alpha_reference()

## a pair fixture over real reference genes:
##   cg0000001  -> OXTR (tier M, measured)   + AVP (tier S3)   co-effect
##   cg0000002  -> NR3C1 (tier S1)
##   cg0000003  -> CRH  (abstained: tissue conflict)
##   cg0000004  -> FKBP5 (tier B)
.pp_pairs <- function() data.frame(
  cpg_id        = c("cg0000001", "cg0000001", "cg0000002", "cg0000003",
                    "cg0000004"),
  target_gene   = c("OXTR", "AVP", "NR3C1", "CRH", "FKBP5"),
  best_direction = c(-1, 1, -1, NA, 1),
  probability_plus1 = c(0.05, 0.9, 0.08, NA, 0.8),
  usable        = c(TRUE, TRUE, TRUE, FALSE, TRUE),
  abstain_reason = c(NA, NA, NA, "tissue_conflict", NA),
  best_evidence = c("measured", "smr_weak", "smr_high", "tissue_conflict",
                    "catalogue_single"),
  direction_tier = c("M", "S3", "S1", "U", "B"),
  stringsAsFactors = FALSE)

.pp_data <- function(n = 60, cols) {
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  d
}

.pp_run <- function(d, cols = NULL, ...) {
  old <- options(dmsa.pair_table = .pp_pairs()); on.exit(options(old))
  dmsa_frame(d, methylation = cols, direction_source = "cpgdirection",
             outcome = "y", covariates = "cov1", blocks = "cID",
             chip = FALSE, B = 19, plots = FALSE, tables = FALSE,
             summary = FALSE, progress = FALSE, beep = FALSE,
             outdir = tempfile("dmsa_pp"), ...)
}

test_that("the pair path builds a frame: one map row per usable pair", {
  set.seed(1)
  cols <- c("cg0000001", "cg0000002", "cg0000003", "cg0000004")
  f <- .pp_run(.pp_data(cols = cols), cols)
  expect_identical(f$direction_source, "cpgdirection")
  ## 4 usable pairs (CRH abstained), incl. the co-effect
  expect_equal(nrow(f$map), 4L)
  expect_setequal(f$map$gene, c("OXTR", "AVP", "NR3C1", "FKBP5"))
})

test_that("co-effects are preserved: one CpG maps to TWO genes as two pairs", {
  set.seed(1)
  cols <- c("cg0000001", "cg0000002")
  f <- .pp_run(.pp_data(cols = cols), cols)
  co <- f$map[f$map$probe == "cg0000001", ]
  expect_equal(nrow(co), 2L)
  expect_setequal(co$gene, c("OXTR", "AVP"))
  ## opposite directions across a CpG's genes are BOTH retained
  expect_setequal(co$best_direction, c(-1, 1))
  ## and both rows point at the same physical measurement
  expect_equal(length(unique(co$column)), 1L)
})

test_that("the ledger records every pair with a reason for the unused", {
  set.seed(1)
  cols <- c("cg0000001", "cg0000002", "cg0000003", "cg0000004")
  f <- .pp_run(.pp_data(cols = cols), cols)
  led <- f$pair_ledger
  expect_equal(nrow(led), 5L)                       # every discovered pair
  crh <- led[led$target_gene == "CRH", ]
  expect_false(crh$used)
  expect_identical(crh$reason_not_used, "tissue_conflict")
  expect_true(all(led$used[led$target_gene != "CRH"]))
})

test_that("spec 59: full analyses every usable pair, confidence only M/S1", {
  set.seed(1)
  cols <- c("cg0000001", "cg0000002", "cg0000004")
  f_full <- .pp_run(.pp_data(cols = cols), cols)               # default full
  f_conf <- .pp_run(.pp_data(cols = cols), cols, cpg_map = "confidence")
  expect_setequal(f_full$map$gene, c("OXTR", "AVP", "NR3C1", "FKBP5"))
  ## confidence: OXTR (M) and NR3C1 (S1) only - a strict subset
  expect_setequal(f_conf$map$gene, c("OXTR", "NR3C1"))
  expect_true(all(paste(f_conf$map$probe, f_conf$map$gene) %in%
                  paste(f_full$map$probe, f_full$map$gene)))
})

test_that("replicate designs of one CpG collapse into one averaged column", {
  set.seed(1)
  cols <- c("cg0000001_TC21_T1_OXTR", "cg0000001_BC21_T1_OXTR", "cg0000002")
  d <- .pp_data(cols = cols)
  f <- .pp_run(d, cols)
  ## the two designs are ONE measurement named by the canonical CpG
  expect_true("cg0000001" %in% f$map$column)
  expect_false(any(cols[1:2] %in% f$map$column))
  m <- f$measurements
  expect_equal(sort(m$replicates[[which(m$cpg_id == "cg0000001")]]), sort(cols[1:2]))
  ## and the frame's M-values equal logit(mean beta) of the two designs
  avg <- rowMeans(cbind(d[[cols[1]]], d[[cols[2]]]))
  want <- log2(pmin(pmax(avg, 1e-4), 1 - 1e-4) /
               (1 - pmin(pmax(avg, 1e-4), 1 - 1e-4)))
  got <- f$M[, which(colnames(f$M) == "cg0000001")[1]]
  expect_equal(unname(got), unname(want), tolerance = 1e-8)
})

test_that("untestable genes are the reference genes with no usable pair", {
  set.seed(1)
  cols <- c("cg0000001", "cg0000002")
  f <- .pp_run(.pp_data(cols = cols), cols)
  ## CRH (abstained upstream, not even discovered for these cols) and every
  ## other reference gene of the touched systems is untestable
  expect_true(all(c("system_id", "system", "gene") %in% names(f$untestable)))
  expect_true(nrow(f$untestable) > 0)
  expect_false(any(f$untestable$gene %in% f$map$gene))
})

test_that("dmsa_coverage reports reference vs testable, and errors on bundled", {
  set.seed(1)
  cols <- c("cg0000001", "cg0000002", "cg0000004")
  f <- .pp_run(.pp_data(cols = cols), cols)
  cv <- dmsa_coverage(f)
  expect_s3_class(cv, "dmsa_coverage")
  expect_true(all(cv$n_genes_testable + cv$n_genes_untestable ==
                  cv$n_genes_reference))
  expect_true(all(cv$family_size_gene == cv$n_genes_testable))
  gv <- dmsa_coverage(f, "gene")
  expect_true(all(gv$testable))
  pv <- dmsa_coverage(f, "pairs")
  expect_equal(nrow(pv), nrow(f$pair_ledger))
  ## a bundled frame has no ledger and says so
  mp0 <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                     package = "dmsa"))
  mp0 <- mp0[is.finite(mp0$best_direction), ]
  bcols <- utils::head(mp0$column[mp0$gene == "FKBP5"], 3)
  db <- .pp_data(cols = bcols)
  fb <- dmsa_frame(db, methylation = bcols, direction_source = "bundled",
                   outcome = "y", covariates = "cov1", blocks = "cID",
                   chip = FALSE, B = 19, plots = FALSE, tables = FALSE,
                   summary = FALSE, progress = FALSE, beep = FALSE,
                   outdir = tempfile("dmsa_ppb"))
  expect_error(dmsa_coverage(fb), "no pair ledger")
})

test_that("the report writes the pair ledger CSV (spec 17)", {
  set.seed(1)
  cols <- c("cg0000001", "cg0000002")
  od <- tempfile("dmsa_pp17")
  old <- options(dmsa.pair_table = .pp_pairs()); on.exit(options(old))
  f <- dmsa_frame(.pp_data(cols = cols), methylation = cols,
                  direction_source = "cpgdirection", outcome = "y",
                  covariates = "cov1", blocks = "cID", chip = FALSE,
                  B = 19, plots = FALSE, tables = TRUE, summary = FALSE,
                  table_type = "html", progress = FALSE, beep = FALSE,
                  outdir = od)
  r <- dmsa_report(f)
  lf <- file.path(od, "tables", "cpg_gene_pair_ledger.csv")
  expect_true(file.exists(lf))
  led <- utils::read.csv(lf)
  expect_equal(nrow(led), nrow(f$pair_ledger))
  expect_true(all(c("pair_id", "used", "reason_not_used") %in% names(led)))
})

test_that("a user map table forces the bundled source", {
  set.seed(1)
  mp0 <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                     package = "dmsa"))
  mp0 <- mp0[is.finite(mp0$best_direction), ]
  sm <- mp0[mp0$gene == "FKBP5", ][1:3, ]
  d <- .pp_data(cols = sm$column)
  f <- dmsa_frame(d, methylation = sm$column, map = sm, outcome = "y",
                  covariates = "cov1", blocks = "cID", chip = FALSE,
                  B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
                  progress = FALSE, beep = FALSE, outdir = tempfile("dmsa_ppu"))
  expect_identical(f$direction_source, "bundled")
  expect_null(f$pair_ledger)
  ## asking for BOTH is contradictory and refused
  expect_error(
    dmsa_frame(d, methylation = sm$column, map = sm,
               direction_source = "cpgdirection", outcome = "y",
               covariates = "cov1", blocks = "cID", chip = FALSE,
               B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
               progress = FALSE, beep = FALSE, outdir = tempfile("x")),
    "pick one source")
})

test_that("non-site columns are never measurements on the pair path", {
  set.seed(1)
  cols <- c("cg0000001", "cg0000002")
  d <- .pp_data(cols = cols)
  d$rs1234 <- stats::plogis(stats::rnorm(nrow(d)))   # a variant probe
  old <- options(dmsa.pair_table = .pp_pairs()); on.exit(options(old))
  f <- dmsa_frame(d, direction_source = "cpgdirection", outcome = "y",
                  covariates = "cov1", blocks = "cID", chip = FALSE,
                  B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
                  progress = FALSE, beep = FALSE, outdir = tempfile("dmsa_pps"))
  expect_false("rs1234" %in% f$measurements$measurement_id)
  expect_false("y" %in% f$measurements$measurement_id)
})

test_that("every submitted CpG has a stated fate (integration-brief rule 4)", {
  set.seed(1)
  ## cg0000009 maps ONLY to a gene outside the reference -> its fate must be
  ## named, not silently absent
  pr <- rbind(.pp_pairs(),
              data.frame(cpg_id = "cg0000009", target_gene = "NOTAREFGENE",
                         best_direction = 1, probability_plus1 = 0.9,
                         usable = TRUE, abstain_reason = NA,
                         best_evidence = "measured", direction_tier = "M",
                         stringsAsFactors = FALSE))
  cols <- c("cg0000001", "cg0000002", "cg0000009")
  d <- .pp_data(cols = cols)
  old <- options(dmsa.pair_table = pr); on.exit(options(old))
  f <- dmsa_frame(d, methylation = cols, direction_source = "cpgdirection",
                  outcome = "y", covariates = "cov1", blocks = "cID",
                  chip = FALSE, B = 19, plots = FALSE, tables = FALSE,
                  summary = FALSE, progress = FALSE, beep = FALSE,
                  outdir = tempfile("dmsa_ppf"))
  expect_true("cg0000009" %in% attr(f$pair_ledger, "cpgs_outside_reference"))
  expect_false("cg0000009" %in% f$map$probe)
  i <- grep("outside reference", f$corrections$field)
  expect_gte(length(i), 1L)
  ## and the coverage print reconciles submitted = mapped + outside + ...
  cv <- dmsa_coverage(f)
  expect_equal(attr(cv, "n_cpgs_mapped") +
               length(attr(cv, "cpgs_outside_reference")),
               attr(cv, "n_cpgs_submitted"))
})

test_that("the adapter's live call is 1-on-1 with cpgdirection >= 2.5", {
  skip_if(is.null(dmsa:::.cpgd("cpg_gene_pairs")),
          "cpgdirection not installed")
  ## the brief's canonical arguments must all exist in the installed function
  fn <- dmsa:::.cpgd("cpg_gene_pairs")
  fm <- names(formals(fn))
  expect_true(all(c("cpgs", "genes", "gene_mode", "annotation_mode",
                    "tissue", "probe_qc", "verbose") %in% fm))
  ## and the live result honours the column contract the adapter validates
  m <- dmsa:::.frame_measurements("cg25140571_TC21_T1_OXTR")
  p <- dmsa:::.frame_cpg_gene_pairs(m, selected_genes = "OXTR")
  expect_true(all(c("cpg_id", "target_gene", "best_direction", "usable",
                    "best_evidence", "direction_tier",
                    "probability_plus1") %in% names(p)))
  ox <- p[p$target_gene == "OXTR", ]
  expect_equal(nrow(ox), 1L)
  expect_equal(ox$best_direction, -1)      # the documented verbatim case
  expect_true(ox$usable)
})

## ---- probe-QC note scoped to the frame (PI ruling, 2026-08-29) ------------
## The exclusion count must answer "of the systems I am analysing", not "of
## everything in the file". Injected tables carry the excluded CpGs' own
## mapping through attr qc_excluded_pairs; without it the note states the
## file-wide scope in words instead of implying a frame scope.

test_that("the QC note is scoped to the analysed systems when mapping is known", {
  set.seed(5)
  cols <- c("cg0000001", "cg0000002")
  d <- .pp_data(cols = cols)
  pt <- .pp_pairs()
  attr(pt, "qc_excluded_cpgs") <- c("cg0000098", "cg0000099")
  ## one excluded CpG maps to a reference gene (OXTR); the other outside
  attr(pt, "qc_excluded_pairs") <- data.frame(
    cpg_id = c("cg0000098", "cg0000099"),
    target_gene = c("OXTR", "NOT_A_GENE"), stringsAsFactors = FALSE)
  old <- options(dmsa.pair_table = pt); on.exit(options(old))
  f <- suppressMessages(dmsa_frame(
    d, methylation = cols, direction_source = "cpgdirection",
    outcome = "y", covariates = "cov1", blocks = "cID",
    chip = FALSE, B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
    progress = FALSE, beep = FALSE, outdir = tempfile("dmsa_qcs")))
  qn <- f$corrections[f$corrections$field == "probe QC", ]
  expect_equal(nrow(qn), 1L)
  expect_match(qn$issue, "^1 CpG\\(s\\) mapping into the analysed system")
  expect_match(qn$action, "2 excluded across the whole submitted file")
  expect_match(qn$action, "cg0000098")     # only the in-scope example named
})

test_that("without pair mapping the QC note states the file-wide scope", {
  set.seed(5)
  cols <- c("cg0000001", "cg0000002")
  d <- .pp_data(cols = cols)
  pt <- .pp_pairs()
  attr(pt, "qc_excluded_cpgs") <- c("cg0000098", "cg0000099")
  old <- options(dmsa.pair_table = pt); on.exit(options(old))
  f <- suppressMessages(dmsa_frame(
    d, methylation = cols, direction_source = "cpgdirection",
    outcome = "y", covariates = "cov1", blocks = "cID",
    chip = FALSE, B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
    progress = FALSE, beep = FALSE, outdir = tempfile("dmsa_qcf")))
  qn <- f$corrections[f$corrections$field == "probe QC", ]
  expect_equal(nrow(qn), 1L)
  expect_match(qn$issue, "across the whole submitted file")
  expect_match(qn$issue, "NOT scoped")
})

test_that("LIVE: a masked in-reference probe is counted in scope only when its system is analysed", {
  skip_if(is.null(tryCatch(getExportedValue("cpgdirection", "cpg_gene_pairs"),
                           error = function(e) NULL)),
          "cpgdirection not installed")
  set.seed(9); n <- 60
  ## cg00242824 is fully masked (Zhou general mask) and maps to TNFRSF1B
  ## (Immune, inflammation & HLA) by the EPIC v2 manifest; cg26261055 ->
  ## CRHBP (HPA) is clean. (cg00143991, used before, reached TNFRSF1A only
  ## through a lookup row on another chromosome, gone in cpgdirectionData
  ## 0.99.5.)
  d <- data.frame(y = rnorm(n), cov1 = rnorm(n), cID = rep(1:30, each = 2))
  for (cl in c("cg26261055", "cg25140571", "cg00242824"))
    d[[cl]] <- stats::plogis(rnorm(n))
  fB <- suppressMessages(dmsa_frame(
    d, methylation = c("cg26261055", "cg25140571", "cg00242824"),
    direction_source = "cpgdirection", outcome = "y", covariates = "cov1",
    blocks = "cID", chip = FALSE, B = 19, plots = FALSE,
    tables = FALSE, summary = FALSE, progress = FALSE, beep = FALSE,
    systems = "HPA axis & glucocorticoid signalling", outdir = tempfile()))
  qB <- fB$corrections[fB$corrections$field == "probe QC", ]
  expect_match(qB$issue, "^0 CpG\\(s\\) mapping into the analysed system")
  expect_match(qB$action, "1 excluded across the whole submitted file")
})
