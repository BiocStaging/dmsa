## Fixture builder standing in for cpgdirection::cpg_gene_pairs().
##
## The column set and semantics mirror the real cpgd_pairs contract: one row
## per CpG x gene pair, direction resolved WITHIN the pair, abstention
## pair-specific, mapping provenance kept separate from direction evidence.
## Injecting this through `.frame_cpg_gene_pairs(pairs = ...)` lets the whole
## architecture be tested where cpgdirection is not installed.
fake_pairs <- function(cpg, genes, d = NA_integer_, usable = NA,
                       mapping_primary = "EPICv2_manifest",
                       evidence = NULL, confidence = NA_real_,
                       abstain = NULL) {
  n <- length(genes)
  d <- rep_len(d, n)
  if (all(is.na(usable))) usable <- !is.na(d)
  usable <- rep_len(usable, n)
  if (is.null(evidence))
    evidence <- ifelse(usable, "catalogue_single", "no_evidence")
  if (is.null(abstain))
    abstain <- ifelse(usable, NA_character_, "no usable direction evidence")
  data.frame(
    cpg_id = rep_len(cpg, n), target_gene = genes,
    mapping_sources = rep_len(mapping_primary, n),
    mapping_primary = rep_len(mapping_primary, n),
    mapping_strength = rep_len("annotation", n),
    best_direction = as.integer(d),
    best_evidence = rep_len(evidence, n),
    best_confidence = rep_len(confidence, n),
    direction_tier = ifelse(usable, "A", NA_character_),
    probability_plus1 = ifelse(is.na(d), NA_real_, ifelse(d > 0, .8, .2)),
    usable = usable, abstain_reason = abstain,
    n_targets_for_cpg = n, is_coeffect = n > 1L,
    stringsAsFactors = FALSE)
}

## A cg25140571 result as cpgdirection returned it before the
## cpgdirectionData 0.99.5 lookups: four discovered targets, exactly one
## usable. TRIM66 was a gene on another chromosome and is no longer a target;
## the fixture is kept because the adapter's handling does not depend on
## which genes they are.
fake_pairs_cg25140571 <- function()
  fake_pairs("cg25140571",
             c("CAV3", "HSA-MIR-548BA", "OXTR", "TRIM66"),
             d = c(NA, NA, -1L, NA),
             usable = c(FALSE, FALSE, TRUE, FALSE),
             mapping_primary = c("EPICv2_manifest", "EPICv2_manifest",
                                 "EPICv2_manifest", "tissue_lookup"),
             evidence = c("distance_tissue_conflict", "no_evidence",
                          "catalogue_single", "no_evidence"))

fake_reference <- function(genes, system_id = "1", system = "Sys A",
                           module_id = "1.1", module = "Mod A1")
  data.frame(gene = genes, system_id = system_id, system = system,
             module_id = module_id, module = module, stringsAsFactors = FALSE)
