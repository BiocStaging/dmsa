test_that("the bundled map is present, versioned, and self-describing", {
  m <- dmsa_directions()
  meta <- attr(m, "dmsa_direction_map")
  expect_true(nrow(m) > 4e5)
  expect_setequal(names(m), c("probe", "gene", "d", "p_plus", "tier"))
  ## an abstention is not a direction: nothing unusable should have shipped
  expect_true(all(m$d %in% c(-1L, 1L)))
  expect_false(any(is.na(m$p_plus)))
  expect_true(all(m$p_plus >= 0 & m$p_plus <= 1))
  ## provenance is the point, not decoration
  expect_true(all(c("tissue", "version", "source", "source_doi") %in%
                    names(meta)))
  expect_identical(meta$tissue, "blood")
})

test_that("lookup by probe, by gene, and by pinned pair each behave", {
  by_gene <- dmsa_directions(genes = "AVP")
  expect_true(nrow(by_gene) > 0)
  expect_true(all(by_gene$gene == "AVP"))

  by_probe <- dmsa_directions(c("cg00052046", "not_a_real_probe"))
  expect_true(all(by_probe$probe == "cg00052046"))
  expect_equal(nrow(dmsa_directions("not_a_real_probe")), 0)

  ## same length probes + genes selects those specific pairs, not the union
  pinned <- dmsa_directions(c("cg00052046", "cg00176879"), c("AVP", "AVP"))
  expect_equal(nrow(pinned), 2)
  expect_true(all(pinned$gene == "AVP"))
})

test_that("dmsa_align takes probe IDs alone and records which map it used", {
  al <- suppressMessages(
    dmsa_align(c("cg00052046", "cg00176879", "cg00308631")))
  expect_s3_class(al, "dmsa_alignment")
  expect_true(nrow(al) >= 3)
  expect_true(all(al$d %in% c(-1, 1)))
  expect_match(attr(al, "direction_map"), "bundled blood map v")

  ## pinning genes restricts to those pairs
  pinned <- suppressMessages(
    dmsa_align(c("cg00052046", "cg00176879"), genes = c("AVP", "AVP")))
  expect_equal(nrow(pinned), 2)
  expect_true(all(pinned$gene == "AVP"))

  ## a probe set with no calls is an error that names the alternative
  expect_error(suppressMessages(dmsa_align(c("nope1", "nope2"))),
               "cpgdirection")
})

test_that("a user-supplied table is untouched by the bundled map", {
  dd <- data.frame(cpg = c("p1", "p2"), d = c(-1, -1), p_plus = c(.2, .2))
  al <- dmsa_align(dd, genes = c("CRH", "NR3C1"), level = "system",
                   polarity = data.frame(gene = c("CRH", "NR3C1"),
                                         w_g = c(1, -1)))
  expect_equal(al$s, c(-1, 1))
  ## and it says so, rather than claiming provenance it does not have
  expect_identical(attr(al, "direction_map"), "user-supplied direction table")
  ## genes remains required for a table
  expect_error(dmsa_align(dd), "`genes` is required")
})

test_that("both bundled layers load, self-describe, and share one schema", {
  for (tis in c("blood", "epithelium")) {
    m <- dmsa_directions(tissue = tis)
    meta <- attr(m, "dmsa_direction_map")
    expect_setequal(names(m), c("probe", "gene", "d", "p_plus", "tier"))
    expect_true(all(m$d %in% c(-1L, 1L)))
    expect_false(any(is.na(m$p_plus)))
    expect_true(all(m$p_plus >= 0 & m$p_plus <= 1))
    expect_identical(meta$tissue, tis)
    expect_identical(meta$build, "hg19")
    expect_true(nzchar(meta$version) && nzchar(meta$source_doi))
    ## the metadata must describe the object it is attached to, not a
    ## sibling: a mislabelled layer is exactly the silent failure that
    ## shipping two layers is meant to prevent
    expect_identical(meta$n_pairs, nrow(m))
    expect_identical(meta$n_probes, length(unique(m$probe)))
  }
  expect_true(nrow(dmsa_directions(tissue = "epithelium")) > 6e5)
})

test_that("no bundled pair links a probe to a gene on another chromosome", {
  ## Map 2026.08 came from catalogue lookups that put many genes on the wrong
  ## chromosome: 84,117 blood and 15,313 epithelium rows paired a probe with a
  ## gene elsewhere in the genome. 2026.11 is rebuilt from the corrected
  ## lookups (cpgdirectionData 0.99.5). Checked here against the curated Alpha
  ## cascade, which records each probe's chromosome and so each reference
  ## gene's; on 2026.08 this test finds 167 blood and 19 epithelium rows.
  cas <- utils::read.csv(gzfile(system.file("extdata", "alpha_cascade_2026c.csv.gz",
                                            package = "dmsa")),
                         stringsAsFactors = FALSE)
  gt <- table(cas$gene, cas$chr)
  gchr <- colnames(gt)[max.col(gt, ties.method = "first")]
  for (tis in c("blood", "epithelium")) {
    m <- dmsa_directions(tissue = tis)
    pc <- cas$chr[match(m$probe, cas$cpg)]
    gc <- gchr[match(m$gene, rownames(gt))]
    ok <- !is.na(pc) & !is.na(gc)
    expect_gt(sum(ok), 1000)
    expect_identical(sum(pc[ok] != gc[ok]), 0L, label = tis)
    expect_identical(attr(m, "dmsa_direction_map")$version, "2026.11")
  }
})

test_that("asking for a layer that does not ship is an error, not a guess", {
  expect_error(dmsa_directions(tissue = "brain"))
  expect_error(dmsa_directions(tissue = "buccal"))
})

test_that("the cache keeps the layers apart", {
  ## order matters: a cache keyed carelessly would hand back whichever
  ## layer was read first, under either name
  e1 <- dmsa_directions(tissue = "epithelium")
  b1 <- dmsa_directions(tissue = "blood")
  e2 <- dmsa_directions(tissue = "epithelium")
  expect_identical(attr(b1, "dmsa_direction_map")$tissue, "blood")
  expect_identical(attr(e2, "dmsa_direction_map")$tissue, "epithelium")
  expect_identical(nrow(e1), nrow(e2))
  expect_false(identical(nrow(b1), nrow(e1)))
})

test_that("the two layers are substantively different maps", {
  ## This is the whole reason two layers ship. If the epithelial calls were
  ## a copy of the blood calls, bundling both would cost 1.6 MB and buy
  ## nothing; if they differ, then running a buccal sample against blood
  ## calls silently mis-signs a large share of its probes.
  b <- dmsa_directions(tissue = "blood")
  e <- dmsa_directions(tissue = "epithelium")
  kb <- paste(b$probe, b$gene, sep = "\r")
  ke <- paste(e$probe, e$gene, sep = "\r")
  shared <- intersect(kb, ke)
  expect_true(length(shared) > 1e5)
  flip <- b$d[match(shared, kb)] != e$d[match(shared, ke)]
  expect_true(mean(flip) > 0.2)

  ## and they are complementary, not nested: each layer calls genes the
  ## other declines to call, so the wrong layer is not merely noisier - it
  ## can have nothing at all to say about the genes in the set
  expect_true(length(setdiff(b$probe, e$probe)) > 5e4)
  expect_true(length(setdiff(e$probe, b$probe)) > 5e4)
})

test_that("dmsa_align reads the layer it was asked for and stamps it", {
  probes <- head(unique(dmsa_directions(tissue = "epithelium")$probe), 20)
  al <- suppressMessages(dmsa_align(probes, tissue = "epithelium"))
  expect_s3_class(al, "dmsa_alignment")
  expect_match(attr(al, "direction_map"), "bundled epithelium map v")
  expect_true(nrow(al) > 0)

  ## a probe called in only one layer must not leak across
  epi_only <- setdiff(dmsa_directions(tissue = "epithelium")$probe,
                      dmsa_directions(tissue = "blood")$probe)
  expect_error(suppressMessages(dmsa_align(head(epi_only, 5))),
               "cpgdirection")
  expect_silent(suppressMessages(
    dmsa_align(head(epi_only, 5), tissue = "epithelium")))
})

test_that("confidence is recovered from every layer a direction table offers", {
  ## Regression: a cpgdirection tissue = "all" result fills best_confidence
  ## only for its catalogue layers. SMR rows arrive with none, got an NA
  ## weight, and then evaporated inside dmsa_triangulate() several steps
  ## later - a silent 61% loss on a real Alpha system, with no error anywhere.
  base <- data.frame(cpg_id = paste0("cg", 1:6),
                     best_direction = c(1, -1, 1, -1, 1, -1),
                     stringsAsFactors = FALSE)
  g <- rep("OXTR", 6)

  ## p_plus wins when present
  d1 <- base; d1$p_plus <- c(.9, .1, .8, .2, .7, .3)
  expect_equal(dmsa_align(d1, genes = g)$p_plus, d1$p_plus)

  ## the single-tissue schema
  d2 <- base; d2$probability_plus1 <- c(.9, .1, .8, .2, .7, .3)
  expect_equal(dmsa_align(d2, genes = g)$p_plus, d2$probability_plus1)

  ## confidence in the CALLED direction converts back to P(+1)
  d3 <- base; d3$best_confidence <- c(.6, .6, .4, .4, 1, 1)
  expect_equal(dmsa_align(d3, genes = g)$p_plus, c(.8, .2, .7, .3, 1, 0))

  ## SMR rows carry no best_confidence; their tier's PUBLISHED accuracy is
  ## the right source, never the universal/distance probability, which would
  ## substitute a 0.60-0.65 prior for a 0.95-0.97 validated causal call
  d4 <- base
  d4$best_confidence <- c(.6, NA, NA, NA, NA, NA)
  d4$smr_tier <- c(NA, "S1", "S2", "S3", "S1", "S2")
  al <- dmsa_align(d4, genes = g)
  expect_false(anyNA(al$p_plus))
  expect_equal(round(al$p_plus, 3), c(.8, .04, .85, .295, .96, .15))

  ## a plain d-only table has no confidence COLUMN: certainty is the
  ## documented behaviour there, and it must not warn
  expect_silent(al5 <- dmsa_align(base, genes = g))
  expect_equal(al5$p_plus, c(1, 0, 1, 0, 1, 0))
})

test_that("a missing confidence VALUE is refused, not assumed certain", {
  ## Filling an empty confidence with certainty would hand the
  ## least-supported rows the most weight - exactly backwards.
  d <- data.frame(cpg_id = paste0("cg", 1:6), best_direction = rep(c(1, -1), 3),
                  best_confidence = c(.8, .8, NA, NA, .6, .6),
                  stringsAsFactors = FALSE)
  expect_warning(al <- dmsa_align(d, genes = rep("OXTR", 6)), "no confidence value")
  expect_equal(sum(al$usable), 4L)
  expect_equal(al$reason[3:4], c("no_confidence", "no_confidence"))
  expect_true(all(!is.na(al$p_plus[al$usable])))
})

test_that("a one-way direction layer is flagged, not silently aligned", {
  ## cpgdirection's distance curves require unanimity and the blood curve
  ## never exceeds 0.449, so distance_only cannot return +1 for any CpG at
  ## any distance. Aligning on a block that points one way by construction
  ## manufactures coherence that is not in the data.
  d <- data.frame(cpg_id = paste0("cg", 1:6), best_direction = rep(-1, 6),
                  best_confidence = rep(.3, 6),
                  best_evidence = c(rep("distance_only", 4), "measured", "measured"),
                  stringsAsFactors = FALSE)
  expect_warning(dmsa_align(d, genes = rep("OXTR", 6)), "ONE-WAY")

  ## and a generic all-one-sign set is flagged whatever its source
  dn <- data.frame(cpg_id = paste0("cg", 1:60), best_direction = rep(-1, 60),
                   p_plus = rep(.2, 60), stringsAsFactors = FALSE)
  expect_warning(dmsa_align(dn, genes = rep("G", 60)), "every one of 60")
})

test_that("no usable probe is dropped between alignment and the test", {
  ## usable and "reaches the test" must be the same set: an NA weight does
  ## not fail in dmsa_align(), it disappears later inside dmsa_triangulate().
  d <- data.frame(cpg_id = paste0("cg", 1:8), best_direction = rep(c(1, -1), 4),
                  best_confidence = c(.6, .5, NA, NA, .8, NA, .3, NA),
                  smr_tier = c(NA, NA, "S1", "S2", NA, "S1", NA, "S3"),
                  stringsAsFactors = FALSE)
  al <- dmsa_align(d, genes = rep("OXTR", 8), level = "gene")
  expect_false(anyNA(al$p_plus))
  expect_true(all(!is.na(al$p_s_plus[al$usable])))
  expect_equal(sum(al$usable), sum(al$usable & !is.na(al$p_s_plus)))
})

test_that("all-one-direction sets: evidence-backed calls get a once-per-set message naming where to check", {
  ## PI 2026-08-29: "If cpgdirection works, it works. This message confuses
  ## users... if the user CAN do something, tell the user what and where."
  probes <- sprintf("cg77%06d", 1:55)
  d_ev <- data.frame(cpg_id = probes, best_direction = -1L,
                     best_confidence = .8, best_evidence = "curated_egene",
                     direction_tier = "M", stringsAsFactors = FALSE)
  g <- rep("GENEX", 55)

  ## evidence-backed: a message (not a warning), pointing at the exact files
  m1 <- capture.output(a1 <- dmsa_align(d_ev, genes = g), type = "message")
  hit <- grep("all 55 direction calls", m1, value = TRUE)
  expect_length(hit, 1)
  expect_match(hit, "analysis_set.csv", fixed = TRUE)
  expect_match(hit, "direction_tier", fixed = TRUE)
  expect_match(hit, "frame\\$map")
  expect_match(hit, "Shown once per set", fixed = TRUE)
  expect_no_warning(dmsa_align(d_ev, genes = g))

  ## and it is throttled: the same set aligned again stays silent
  m2 <- capture.output(a2 <- dmsa_align(d_ev, genes = g), type = "message")
  expect_length(grep("direction calls", m2), 0)

  ## the throttle is per SET: a different probe block speaks again
  probes3 <- sprintf("cg88%06d", 1:55)
  d_ev3 <- d_ev; d_ev3$cpg_id <- probes3
  m3 <- capture.output(a3 <- dmsa_align(d_ev3, genes = g), type = "message")
  expect_length(grep("all 55 direction calls", m3), 1)
})

test_that("all-one-direction sets WITHOUT evidence tiers keep an actionable warning", {
  probes <- sprintf("cg99%06d", 1:52)
  d_raw <- data.frame(cpg_id = probes, d = 1L, stringsAsFactors = FALSE)
  g <- rep("GENEY", 52)
  w <- tryCatch(withCallingHandlers(
    { dmsa_align(d_raw, genes = g); NULL },
    warning = function(w) { stop(conditionMessage(w)) }),
    error = function(e) conditionMessage(e))
  expect_match(w, "every one of 52")
  expect_match(w, "NO evidence tiers")
  expect_match(w, "table\\(map\\$best_direction\\)")
  ## throttled on repeat too
  expect_no_warning(dmsa_align(d_raw, genes = g))

  ## a mixed-direction set triggers neither branch
  d_mix <- data.frame(cpg_id = sprintf("cg66%06d", 1:60),
                      d = rep(c(1L, -1L), 30), stringsAsFactors = FALSE)
  m <- capture.output(
    expect_no_warning(dmsa_align(d_mix, genes = rep("GENEZ", 60))),
    type = "message")
  expect_length(grep("direction calls", m), 0)
})
