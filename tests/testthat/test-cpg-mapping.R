## Spec 6-11, 18, 20, 25: CpG -> gene mapping, co-effects, replicates, filter.

test_that("spec 20: only methylation sites are analysed, and the drop is reported", {
  cols <- c("cg00645497_BC21_T1_CD38",
            "nv-GRCh38-chr19-3118944-3118944-A-C_TC11_T1_GNA11",
            "cg25140571_BC21_T1_OXTR")
  m <- dmsa:::.frame_measurements(cols)
  expect_equal(nrow(m), 2L)
  expect_length(attr(m, "dropped_nonsite"), 1L)
  expect_match(attr(m, "dropped_nonsite"), "^nv-")
})

test_that("spec 20: replicate DESIGNS of one CpG collapse; co-effects do not", {
  cols <- c("cg26890182_BC11_T1_ACMSD", "cg26890182_BC12_T1_ACMSD",
            "cg26890182_BC13_T1_ACMSD", "cg25140571_BC21_T1_OXTR")
  m <- dmsa:::.frame_measurements(cols)                       # default collapse_site
  expect_equal(nrow(m), 2L)
  expect_setequal(m$cpg_id, c("cg26890182", "cg25140571"))
  expect_equal(attr(m, "replicate_cpgs"), "cg26890182")
  expect_length(m$replicates[[which(m$cpg_id == "cg26890182")]], 3L)
  ## keeping them is an explicit choice, never a silent one
  k <- dmsa:::.frame_measurements(cols, replicate_policy = "keep_correlated")
  expect_equal(nrow(k), 4L)
})

test_that("spec 10: the column suffix is provenance, never a restriction", {
  m <- dmsa:::.frame_measurements("cg25140571_BC21_T1_OXTR")
  expect_equal(m$cpg_id, "cg25140571")
  expect_equal(m$input_gene, "OXTR")
  ## the suffix says OXTR, but all four discovered targets are offered
  p <- dmsa:::.frame_cpg_gene_pairs(m, c("CAV3", "OXTR", "TRIM66", "HSA-MIR-548BA"),
                             pairs = fake_pairs_cg25140571())
  expect_setequal(p$target_gene, c("CAV3", "OXTR", "TRIM66", "HSA-MIR-548BA"))
})

test_that("spec 56: a BARE CpG discovers its targets, no suffix needed", {
  m <- dmsa:::.frame_measurements("cg1")
  p <- dmsa:::.frame_cpg_gene_pairs(m, c("GENE_A", "GENE_B"),
                             pairs = fake_pairs("cg1", c("GENE_A", "GENE_B"),
                                                d = c(1L, -1L)))
  expect_equal(nrow(p), 2L)                    # 1 measurement, 1 CpG, 2 pairs
  expect_equal(length(unique(p$measurement_id)), 1L)
  expect_setequal(p$target_gene, c("GENE_A", "GENE_B"))
  expect_true(all(p$is_coeffect_selected))
})

test_that("spec 57 + 25: selected genes FILTER; a resource cannot expand them", {
  m <- dmsa:::.frame_measurements("cg1")
  fx <- fake_pairs("cg1", c("A", "B", "C"), d = 1L)
  p <- dmsa:::.frame_cpg_gene_pairs(m, c("A", "C"), pairs = fx)
  expect_setequal(p$target_gene, c("A", "C"))
  expect_equal(attr(p, "genes_outside_reference"), "B")
})

test_that("spec 18/11: one CpG feeds several genes; nothing dedups it away", {
  m <- dmsa:::.frame_measurements(c("cg1_BC21_T1_A", "cg2_BC21_T1_C"))
  fx <- rbind(fake_pairs("cg1", c("A", "B"), d = 1L),
              fake_pairs("cg2", "C", d = -1L))
  p <- dmsa:::.frame_cpg_gene_pairs(m, c("A", "B", "C"), pairs = fx)
  expect_equal(nrow(p), 3L)
  expect_equal(sum(p$cpg_id == "cg1"), 2L)
  expect_equal(p$n_targets_selected[p$cpg_id == "cg2"], 1L)
  expect_false(p$is_coeffect_selected[p$cpg_id == "cg2"])
})

test_that("direction and abstention are PAIR-specific, not CpG-specific", {
  m <- dmsa:::.frame_measurements("cg25140571")
  p <- dmsa:::.frame_cpg_gene_pairs(m, c("CAV3", "HSA-MIR-548BA", "OXTR", "TRIM66"),
                             pairs = fake_pairs_cg25140571())
  expect_equal(sum(p$usable), 1L)
  expect_equal(p$target_gene[p$usable], "OXTR")
  expect_equal(p$best_direction[p$usable], -1L)
  expect_true(all(is.na(p$best_direction[!p$usable])))
})

test_that("spec 58: pair invariance - provenance may differ, the record may not", {
  m <- dmsa:::.frame_measurements("cg1")
  bare   <- fake_pairs("cg1", "A", d = -1L, mapping_primary = "EPICv2_manifest")
  suffix <- fake_pairs("cg1", "A", d = -1L, mapping_primary = "input_annotation")
  pb <- dmsa:::.frame_cpg_gene_pairs(m, "A", pairs = bare)
  ps <- dmsa:::.frame_cpg_gene_pairs(m, "A", pairs = suffix)
  key <- c("cpg_id", "target_gene", "best_direction", "usable")
  expect_equal(pb[, key], ps[, key])
  expect_false(identical(pb$mapping_primary, ps$mapping_primary))
})

test_that("spec 9/16/17: the ledger keys on the PAIR and explains exclusions", {
  m <- dmsa:::.frame_measurements("cg25140571_BC21_T1_OXTR")
  p <- dmsa:::.frame_cpg_gene_pairs(m, c("CAV3", "HSA-MIR-548BA", "OXTR", "TRIM66"),
                             pairs = fake_pairs_cg25140571())
  led <- dmsa:::.frame_pair_ledger(p, fake_reference(c("CAV3", "HSA-MIR-548BA",
                                                "OXTR", "TRIM66")))
  expect_equal(nrow(led), 4L)
  expect_equal(anyDuplicated(led$pair_id), 0L)
  expect_true(all(grepl("::", led$pair_id)))
  expect_equal(sum(led$used), 1L)
  ## every excluded pair says why
  expect_false(any(is.na(led$reason_not_used[!led$used])))
  expect_true(all(c("system", "module", "mapping_primary", "best_evidence")
                  %in% names(led)))
})

test_that("mapping needs cpg_gene_pairs(), and says so plainly when absent", {
  ## the condition is whether the FUNCTION resolves, not whether the package is
  ## installed: cpgdirection 2.3.0 is present but predates cpg_gene_pairs(),
  ## and a version check that only asked "is the package there?" would have
  ## claimed the new API was available.
  skip_if(!is.null(dmsa:::.cpgd("cpg_gene_pairs")),
          "cpg_gene_pairs() is available")
  m <- dmsa:::.frame_measurements("cg1")
  e <- expect_error(dmsa:::.frame_cpg_gene_pairs(m, "A"))
  expect_match(conditionMessage(e), "cpgdirection")
  expect_match(conditionMessage(e), "cpg_gene_pairs")
})

test_that("a pair table missing required columns is refused", {
  m <- dmsa:::.frame_measurements("cg1")
  bad <- data.frame(cpg_id = "cg1", target_gene = "A")
  expect_error(dmsa:::.frame_cpg_gene_pairs(m, "A", pairs = bad),
               "missing required column")
})

## ---------------------------------------------------------------------------
## LIVE CONTRACT: run against the real cpgdirection when it is installed.
## Everything above uses the injection seam so the architecture is tested
## everywhere; this block additionally proves the adapter's expectations match
## the real cpg_gene_pairs(). It skips where the function is absent, which is
## the only thing that legitimately cannot be tested without it.
## ---------------------------------------------------------------------------

test_that("LIVE: the adapter's argument names all exist on cpg_gene_pairs()", {
  fn <- dmsa:::.cpgd("cpg_gene_pairs")
  skip_if(is.null(fn), "cpg_gene_pairs() not available")
  used <- c("cpgs", "genes", "gene_mode", "annotation_mode", "tissue",
            "include_brain", "probe_qc", "verbose")
  expect_length(setdiff(used, names(formals(fn))), 0L)
})

test_that("LIVE: real pairs carry every column the adapter requires", {
  fn <- dmsa:::.cpgd("cpg_gene_pairs")
  skip_if(is.null(fn), "cpg_gene_pairs() not available")
  p <- fn(cpgs = "cg25140571", universal = FALSE, verbose = FALSE)
  expect_true(all(c("cpg_id", "target_gene", "best_direction", "usable")
                  %in% names(p)))
  expect_true(all(c("mapping_primary", "best_evidence", "abstain_reason",
                    "n_targets_for_cpg", "is_coeffect") %in% names(p)))
  expect_s3_class(p, "cpgd_pairs")
})

test_that("LIVE: co-effects and pair-specific direction survive the adapter", {
  fn <- dmsa:::.cpgd("cpg_gene_pairs")
  skip_if(is.null(fn), "cpg_gene_pairs() absent")
  ## cpgdirection's own answer for every target of cg25140571 in blood
  q <- as.data.frame(fn(cpgs = "cg25140571", tissue = "blood",
                        universal = FALSE, verbose = FALSE))
  key <- function(x, g) paste(x$usable[x$target_gene == g],
                              x$best_direction[x$target_gene == g])
  ## a second target whose answer differs from OXTR's (CAV3 with the
  ## cpgdirectionData 0.99.5 lookups; TRIM66, used before, was a gene on
  ## another chromosome and is no longer a target)
  other <- setdiff(q$target_gene, "OXTR")
  other <- other[vapply(other, function(g) key(q, g) != key(q, "OXTR"), logical(1))]
  skip_if(!length(other), "no target of cg25140571 answers differently from OXTR")
  g2 <- if ("CAV3" %in% other) "CAV3" else other[1]
  m <- dmsa:::.frame_measurements("cg25140571_BC21_T1_OXTR")
  p <- as.data.frame(dmsa:::.frame_cpg_gene_pairs(m, c("OXTR", g2), tissue = "blood"))
  ## one CpG, two selected target genes -> two pairs, flagged as a co-effect
  expect_equal(nrow(p), 2L)
  expect_true(all(p$is_coeffect_selected))
  ## direction is resolved WITHIN the pair, as cpgdirection resolves it
  for (g in c("OXTR", g2)) {
    expect_identical(p$usable[p$target_gene == g], q$usable[q$target_gene == g], info = g)
    expect_equal(p$best_direction[p$target_gene == g], q$best_direction[q$target_gene == g],
                 info = g)
  }
  expect_true(p$usable[p$target_gene == "OXTR"])
  expect_equal(p$best_direction[p$target_gene == "OXTR"], -1L)
})

test_that("LIVE: the gene filter is never exceeded", {
  skip_if(is.null(dmsa:::.cpgd("cpg_gene_pairs")), "cpg_gene_pairs() absent")
  ## cg25140571 also supports CAV3 and HSA-MIR-548BA; neither is selected here
  m <- dmsa:::.frame_measurements("cg25140571")
  p <- dmsa:::.frame_cpg_gene_pairs(m, "OXTR", tissue = "blood")
  expect_setequal(unique(p$target_gene), "OXTR")
})
