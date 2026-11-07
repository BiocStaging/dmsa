## ===========================================================================
## CpG -> GENE MAPPING  (spec sections 1, 6-11, 20, 25)
##
## THE INVARIANT THIS FILE EXISTS TO ENFORCE
##
##   BIOLOGICAL QUESTION        system > module > gene      (dmsa_reference)
##   MEASUREMENT AVAILABILITY   the user's actual CpGs      (their data)
##   TARGET DISCOVERY+DIRECTION cpgdirection::cpg_gene_pairs()
##   ANALYTICAL OBJECT          measurement x cpg_id x target_gene
##
## DMSA never pre-selects CpGs. It reads whatever methylation the user
## supplies, canonicalises each measurement to a CpG id, and asks cpgdirection
## which genes that CpG supports and in which direction - separately for every
## CpG x gene pair. The selected genes are a FILTER over that discovery, never
## a requirement and never an expander: a direction resource may not add a
## biological gene the user did not select (spec 25).
##
## ONE CpG MAPS TO MANY GENES. That is the biology, not a duplicate to clean
## up: cg25140571 supports CAV3, HSA-MIR-548BA, OXTR and TRIM66, and only the
## OXTR pair carries a usable direction. Deduplicating by cpg_id or by data
## column destroys exactly the co-effects the analysis is built on, so the
## uniqueness key here is the PAIR, never the CpG (spec 11, 18).
##
## Co-effects (one CpG -> many genes) and technical replicates (many probe
## designs -> one CpG) are DIFFERENT problems and are never solved by one
## duplicated() call (spec 20).
## ===========================================================================

## Canonical CpG id from a user's measurement column.
## EPIC v2 columns arrive as cg########_<design>_<time>_<GENE>, e.g.
## "cg00645497_BC21_T1_CD38". The CpG is the leading cg id; the trailing gene
## is INPUT PROVENANCE only - it never restricts the CpG to that gene, because
## cpgdirection's augment mode can legitimately return further targets
## (spec 10).
.map_canon_cpg <- function(x) sub("^(cg[0-9]+).*$", "\\1", as.character(x))

.map_input_gene <- function(x) {
  x <- as.character(x)
  g <- sub("^cg[0-9]+_[^_]+_[^_]+_(.+)$", "\\1", x)
  ifelse(g == x, NA_character_, g)
}

## Is this measurement a methylation SITE? PI ruling: DMSA analyses cg probes
## only. EPIC v2 also ships variant probes (nv-GRCh38-chr19-...-A-C_TC11_T1_GNA11,
## and rs/ch designs); they are not methylation sites and are excluded BY RULE,
## with a reported count - never silently.
.map_is_site <- function(x) grepl("^cg[0-9]+", as.character(x))

## Build the measurement table from the user's methylation column names.
##
## Returns one row per RETAINED measurement with
##   measurement_id  under "collapse_site" the canonical CpG id (the policy
##                   collapses designs, so the site is the measurement);
##                   under "keep_correlated" the exact user data column.
##                   `replicates` always lists the physical column(s).
##   cpg_id          the canonical CpG
##   input_gene      gene parsed from the column suffix, or NA
##   n_designs       how many columns measure this same CpG
## and carries, as attributes, what was dropped and why.
##
## replicate_policy (spec 20):
##   "collapse_site"    [default] several probe designs of ONE CpG are one
##                      measurement; the replicates are averaged before
##                      analysis, because they are repeated reads of the same
##                      site and counting them separately inflates that site.
##   "keep_correlated"  keep every design as its own measurement and accept
##                      that they are correlated.
.frame_measurements <- function(columns,
                                replicate_policy = c("collapse_site",
                                                     "keep_correlated")) {
  replicate_policy <- match.arg(replicate_policy)
  columns <- as.character(columns)
  site <- .map_is_site(columns)
  dropped_nonsite <- columns[!site]
  keep <- columns[site]
  if (!length(keep))
    stop("no methylation-site (cg) columns found among ", length(columns),
         " supplied measurement(s). DMSA analyses CpG methylation sites; ",
         "variant/control probes (nv-, rs-, ch-) are not sites.",
         call. = FALSE)

  m <- data.frame(measurement_id = keep,
                  cpg_id = .map_canon_cpg(keep),
                  input_gene = .map_input_gene(keep),
                  stringsAsFactors = FALSE)
  n <- table(m$cpg_id)
  m$n_designs <- as.integer(n[m$cpg_id])
  reps <- names(n)[n > 1L]

  if (replicate_policy == "collapse_site") {
    ## one row per canonical CpG; the members travel with it so the caller can
    ## average the right columns and the ledger can say which they were.
    ## Under this policy measurement_id IS the canonical CpG - for every row,
    ## replicated or not. Rewriting only when some OTHER CpG happened to be
    ## replicated made the same input column mint different pair_ids across
    ## datasets.
    grp <- split(m$measurement_id, m$cpg_id)
    m <- m[!duplicated(m$cpg_id), , drop = FALSE]
    m$replicates <- I(unname(grp[m$cpg_id]))
    m$measurement_id <- m$cpg_id
  } else {
    m$replicates <- I(as.list(m$measurement_id))
  }
  rownames(m) <- NULL
  attr(m, "dropped_nonsite") <- dropped_nonsite
  attr(m, "replicate_cpgs")  <- reps
  attr(m, "replicate_policy") <- replicate_policy
  m
}

## The cpgdirection adapter (spec 6, 7).
##
## `pairs` is the INJECTION SEAM. Production passes NULL and the real
## cpgdirection::cpg_gene_pairs() is resolved at call time; tests pass a
## fixture table so the architecture is exercised in builds where cpgdirection
## is not installed - the alternative, skip_if_not_installed(), would hide the
## core of the package from its own build report.
##
## Resolution is by name rather than `::` on purpose: cpgdirection ships from
## GitHub/Zenodo, and Bioconductor forbids declaring a dependency its build
## system cannot install.
##
## THE SWITCH, when cpgdirection is accepted into Bioconductor:
##   1. DESCRIPTION:  Imports: ..., cpgdirection
##   2. NAMESPACE:    importFrom(cpgdirection, cpg_gene_pairs)
##   3. here:         replace the .cpgd("cpg_gene_pairs") lookup and its
##                    is.null() guard with a direct cpg_gene_pairs(...) call.
## Nothing else in DMSA changes, because nothing else calls cpgdirection.
##
## DMSA depends on cpgdirection, NOT on cpgdirectionData. The planned layout is
##   cpgdirectionData  ExperimentHub package holding the large layers
##   cpgdirection      thin software package; retrieves them from the Hub
##   dmsa              Imports cpgdirection
## so the data package is cpgdirection's dependency to declare, not DMSA's.
## Wiring DMSA straight to cpgdirectionData would bypass the direction ladder
## that resolves a pair, which is the whole reason cpgdirection exists.
.frame_cpg_gene_pairs <- function(measurements, selected_genes,
                                  tissue = "blood",
                                  annotation_mode = c("augment", "strict"),
                                  include_brain = FALSE,
                                  probe_qc = c("exclude", "flag", "ignore"),
                                  pairs = NULL) {
  annotation_mode <- match.arg(annotation_mode)
  probe_qc <- match.arg(probe_qc)
  cpgs <- unique(measurements$cpg_id)

  .live_fn <- NULL
  if (is.null(pairs)) {
    fn <- .cpgd("cpg_gene_pairs")
    if (is.null(fn))
      stop("CpG-to-gene mapping needs the cpgdirection package, which is not ",
           "installed.\nDMSA reads which genes a CpG supports, and in which ",
           "direction, from cpgdirection::cpg_gene_pairs().\nInstall it, or ",
           "pass a pair table directly for testing.", call. = FALSE)
    ## cpgdirection announces its QC exclusion count over the WHOLE submitted
    ## file. dmsa_frame() reports the same fact scoped to the analysed
    ## systems (a file-wide count beside a two-system frame reads as the
    ## frame's own loss - PI, 2026-08-29), so the file-scope line is muffled
    ## here and the frame prints the reconciled, scoped note instead.
    pairs <- withCallingHandlers(
      fn(cpgs = cpgs, genes = selected_genes, gene_mode = "filter",
         annotation_mode = annotation_mode, tissue = tissue,
         include_brain = include_brain, probe_qc = probe_qc,
         verbose = FALSE),
      message = function(m) invokeRestart("muffleMessage"))
    .live_fn <- fn
  }
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)

  req <- c("cpg_id", "target_gene", "best_direction", "usable")
  miss <- setdiff(req, names(pairs))
  if (length(miss))
    stop("the CpG-gene pair table is missing required column(s): ",
         paste(miss, collapse = ", "),
         "\ncpg_gene_pairs() is expected to return one row per CpG x gene ",
         "pair with at least: ", paste(req, collapse = ", "), call. = FALSE)
  for (v in c("cpg_id", "target_gene")) pairs[[v]] <- as.character(pairs[[v]])

  ## spec 25: the selection is authoritative. gene_mode = "filter" should have
  ## done this already, but a mapping resource must never be able to enlarge
  ## the biological question, so the intersection is enforced here too and what
  ## it removed is recorded rather than dropped in silence.
  outside <- setdiff(unique(pairs$target_gene), selected_genes)
  .cpg_before <- unique(pairs$cpg_id)
  pairs <- pairs[pairs$target_gene %in% selected_genes, , drop = FALSE]
  ## a CpG whose EVERY discovered target lies outside the reference is a
  ## fourth fate, distinct from unmapped and QC-excluded: it was mapped fine,
  ## just not to the biology being asked about. Under gene_mode = "filter"
  ## cpgdirection filters INTERNALLY, so such a CpG returns no row at all and
  ## is invisible to a before/after diff - it is recovered as the residual:
  ## submitted - surviving - unmapped - QC-excluded. Both routes are taken so
  ## injected pair tables (which carry no attributes) reconcile too.
  cpgs_outside <- setdiff(.cpg_before, unique(pairs$cpg_id))
  .surv <- unique(pairs$cpg_id)
  .attr_unm <- attr(pairs, "unmapped_cpgs")
  .attr_qcx <- attr(pairs, "qc_excluded_cpgs")
  if (!is.null(attr(pairs, "n_submitted")))
    cpgs_outside <- union(cpgs_outside,
                          setdiff(cpgs, c(.surv, .attr_unm, .attr_qcx)))

  ## merge pairs back onto the PHYSICAL measurements (spec 9): several designs
  ## may share one canonical CpG, and one CpG may support several genes, so
  ## this join is deliberately many-to-many.
  out <- merge(measurements[, setdiff(names(measurements), "replicates"),
                            drop = FALSE],
               pairs, by = "cpg_id", all.x = FALSE, all.y = FALSE)
  if (!nrow(out))
    stop("no CpG x gene pair survived: none of the ", length(cpgs),
         " canonical CpG(s) in your data maps to any of the ",
         length(selected_genes), " selected gene(s).", call. = FALSE)

  ## co-effect accounting is computed over the SELECTED genes, so it answers
  ## "how many of the genes I am analysing does this CpG serve", not "how many
  ## genes exist in the annotation"
  tg <- tapply(out$target_gene, out$cpg_id, function(z) length(unique(z)))
  out$n_targets_selected <- as.integer(tg[out$cpg_id])
  out$is_coeffect_selected <- out$n_targets_selected > 1L
  out <- out[order(out$cpg_id, out$target_gene), , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "genes_outside_reference") <- outside
  attr(out, "cpgs_outside_reference") <- cpgs_outside
  attr(out, "n_cpgs_submitted") <- length(cpgs)
  attr(out, "n_cpgs_mapped") <- length(unique(out$cpg_id))
  ## rule 4 of the integration brief: probe QC removes CpGs SILENTLY by
  ## default (~3.85% of the array) and unmapped CpGs return no row at all.
  ## Both lists ride along so the frame can reconcile counts and report them
  ## rather than letting the loss read as a bug - or worse, go unread.
  for (at in c("qc_excluded_cpgs", "unmapped_cpgs", "source_counts",
               "tissue", "gene_mode", "probe_qc", "qc_excluded_pairs"))
    if (!is.null(attr(pairs, at))) attr(out, at) <- attr(pairs, at)
  ## WHERE the QC-excluded CpGs would have mapped (PI, 2026-08-29): under
  ## probe_qc = "exclude" a masked CpG returns no row, so the frame cannot
  ## say how many of the excluded CpGs even belonged to the systems being
  ## analysed - the count reads as the frame's own loss when it is mostly
  ## the rest of the file's. One extra call over JUST the excluded CpGs, in
  ## flag mode, recovers their target genes; the analysis set is untouched.
  .qcx2 <- attr(out, "qc_excluded_cpgs")
  if (is.null(attr(out, "qc_excluded_pairs")) && length(.qcx2) &&
      !is.null(.live_fn)) {
    qp <- tryCatch(withCallingHandlers(
      .live_fn(cpgs = .qcx2, genes = selected_genes, gene_mode = "filter",
               annotation_mode = annotation_mode, tissue = tissue,
               include_brain = include_brain, probe_qc = "flag",
               verbose = FALSE),
      message = function(m) invokeRestart("muffleMessage")),
      error = function(e) NULL)
    if (!is.null(qp)) {
      qp <- as.data.frame(qp, stringsAsFactors = FALSE)
      if (all(c("cpg_id", "target_gene") %in% names(qp)))
        attr(out, "qc_excluded_pairs") <-
          unique(qp[qp$target_gene %in% selected_genes,
                    c("cpg_id", "target_gene"), drop = FALSE])
    }
  }
  out
}

## Attach the biological reference (system, module) and mint the pair key.
## pair_id = measurement_id::system_id::gene is collision-safe for the 2026c
## codebook, where every gene sits in exactly one system and one module.
.frame_pair_ledger <- function(pairs, reference) {
  ref <- if (inherits(reference, "dmsa_reference")) reference$systems else
    as.data.frame(reference, stringsAsFactors = FALSE)
  keep <- intersect(c("gene", "system_id", "system", "module_id", "module"),
                    names(ref))
  led <- merge(pairs, unique(ref[, keep, drop = FALSE]),
               by.x = "target_gene", by.y = "gene", all.x = TRUE)
  ## cpgdirection ships its OWN `pair_id` ("<cpg>|<gene>"), which is the stable
  ## key of its published column contract. Minting DMSA's key over the top of
  ## it destroyed that provenance silently: two different identifiers, one name.
  ## Keep theirs under its source name and mint ours beside it.
  if ("pair_id" %in% names(led)) {
    led$cpgd_pair_id <- as.character(led$pair_id)
    led$pair_id <- NULL
  }
  led$pair_id <- paste(led$measurement_id, led$system_id, led$target_gene,
                       sep = "::")
  if (anyDuplicated(led$pair_id))
    warning("pair_id is not unique; a gene appears in more than one module of ",
            "the same system in this reference")
  led$used <- as.logical(led$usable) %in% TRUE
  .why <- if ("abstain_reason" %in% names(led))
    as.character(led$abstain_reason) else rep(NA_character_, nrow(led))
  .why[is.na(.why) | !nzchar(.why)] <- "no usable pair-specific direction"
  led$reason_not_used <- ifelse(led$used, NA_character_, .why)
  ord <- c("pair_id", "cpgd_pair_id",
           "measurement_id", "cpg_id", "target_gene", "system_id",
           "system", "module_id", "module", "mapping_sources", "mapping_primary",
           "mapping_strength", "best_direction", "best_evidence",
           "best_accuracy", "best_confidence", "direction_tier", "probability_plus1",
           "n_targets_selected", "is_coeffect_selected", "used",
           "reason_not_used")
  led <- led[, c(intersect(ord, names(led)),
                 setdiff(names(led), ord)), drop = FALSE]
  led <- led[order(led$system_id, led$target_gene, led$cpg_id), , drop = FALSE]
  rownames(led) <- NULL
  ## merge() drops attributes; the provenance counts belong on the ledger too
  for (at in c("genes_outside_reference", "cpgs_outside_reference",
               "n_cpgs_submitted", "n_cpgs_mapped",
               "qc_excluded_cpgs", "unmapped_cpgs", "source_counts",
               "tissue", "gene_mode", "probe_qc", "qc_excluded_pairs"))
    if (!is.null(attr(pairs, at))) attr(led, at) <- attr(pairs, at)
  led
}

## ---------------------------------------------------------------------------
## Spec 1: the runtime that dmsa_frame() calls when direction comes from
## CpG x gene PAIRS (cpgdirection, or an injected pair table) rather than from
## the bundled per-column snapshot.
##
## Order (the spec's): measurements -> pair discovery -> reference filter ->
## ledger -> the engine's map. The ENGINE is untouched: what changes is where
## `mp` comes from, and mp keeps exactly the columns the engine has always
## consumed (probe, column, gene, system_id, system, best_direction, p_plus)
## plus provenance the report may show.
##
## p_plus is resolved by dmsa_align()'s OWN documented rule, by handing it
## the raw pair columns once - not by re-implementing it here: best_accuracy
## when cpg_gene_pairs() returns it (cpgdirection's November 2026 builds, E17), otherwise the
## ladder p_plus -> probability_plus1 -> best_confidence -> smr-tier
## midpoints -> unusable.
.frame_pairs_runtime <- function(available_columns, reference,
                                 tissue = "blood",
                                 replicate_policy = "collapse_site",
                                 probe_qc = "exclude",
                                 pairs = NULL) {
  ref_sys <- if (inherits(reference, "dmsa_reference")) reference$systems
             else as.data.frame(reference, stringsAsFactors = FALSE)
  need <- c("gene", "system_id", "system")
  miss <- setdiff(need, names(ref_sys))
  if (length(miss))
    stop("the biological reference is missing column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)

  ## 1. measurements: only cg site columns are candidates; everything else in
  ##    `data` (outcomes, covariates, ids) is simply not a measurement
  cand <- available_columns[.map_is_site(available_columns)]
  if (!length(cand))
    stop("no CpG methylation-site (cg...) columns found among the ",
         length(available_columns), " available column(s). DMSA reads which ",
         "genes a CpG supports from the CpG ids themselves, so the ",
         "methylation columns must be named by their cg identifiers.",
         call. = FALSE)
  m <- .frame_measurements(cand, replicate_policy = replicate_policy)

  ## 2. + 3. pair discovery, with the reference's genes as the FILTER
  p <- .frame_cpg_gene_pairs(m, selected_genes = unique(ref_sys$gene),
                             tissue = tissue, probe_qc = probe_qc,
                             pairs = pairs)

  ## 4. the ledger: every pair, why used or not
  led <- .frame_pair_ledger(p, reference)

  ## 5. resolve p_plus through the engine's own ladder, once, on every pair
  ##    that carries a direction (usable or not - the ledger reports both)
  dir_tab <- data.frame(cpg_id = led$cpg_id,
                        best_direction = led$best_direction,
                        stringsAsFactors = FALSE)
  for (cn in c("best_accuracy", "probability_plus1", "best_confidence", "smr_tier"))
    if (cn %in% names(led)) dir_tab[[cn]] <- led[[cn]]
  al0 <- suppressWarnings(withCallingHandlers(
    dmsa_align(dir_tab, genes = led$target_gene, level = "gene"),
    message = function(m) invokeRestart("muffleMessage")))
  led$p_plus <- al0$p_plus[match(paste(led$cpg_id, led$target_gene),
                                 paste(al0$probe, al0$gene))]
  ## a pair whose confidence could not be recovered is NOT usable at the
  ## frame level, whatever `usable` said: with an NA weight it would
  ## evaporate inside the engine instead of being reported here
  no_conf <- led$used & !is.finite(led$p_plus)
  if (any(no_conf)) {
    led$used[no_conf] <- FALSE
    led$reason_not_used[no_conf] <- "direction has no recoverable confidence"
  }

  ## 6. the engine's map: one row per USABLE pair
  u <- led[led$used, , drop = FALSE]
  if (!nrow(u))
    stop("no usable CpG x gene pair remains: ", nrow(led), " pair(s) were ",
         "discovered for ", length(unique(led$cpg_id)), " CpG(s), but none ",
         "carries a usable pair-specific direction. tables in the pair ",
         "ledger say why, pair by pair.", call. = FALSE)
  mp <- data.frame(probe = u$cpg_id,
                   column = u$measurement_id,
                   gene = u$target_gene,
                   system_id = as.character(u$system_id),
                   system = u$system,
                   best_direction = u$best_direction,
                   p_plus = u$p_plus,
                   stringsAsFactors = FALSE)
  for (cn in c("module_id", "module", "best_evidence", "direction_tier",
               "best_rung", "best_accuracy",
               "pair_id", "cpgd_pair_id", "is_coeffect_selected",
               "n_targets_selected"))
    if (cn %in% names(u)) mp[[cn]] <- u[[cn]]

  ## spec 13/14: pair-level meaning of full / confidence. "full" is every
  ## usable pair. "confidence" is the strict subset whose direction rests on
  ## the top of the evidence ladder: tier M (measured in the target tissue)
  ## or S1 (top SMR). It is a SENSITIVITY analysis - it never redefines
  ## membership, and a pair in "confidence" is always in "full".
  hi <- if ("direction_tier" %in% names(mp))
    mp$direction_tier %in% c("M", "S1") else rep(FALSE, nrow(mp))
  list(full = mp,
       confidence = mp[hi, , drop = FALSE],
       ledger = led,
       measurements = m)
}

## ---------------------------------------------------------------------------
## CpG-set selection by NAME PATTERN (PI request, 2026-08-29): one methylation
## file often carries several biological groups in its column names - in
## Alpha's child file, "_maternal" and "_paternal" mark the parents' probes
## and the child's carry no group tag. Which group a question is about is a
## SELECTION, and it must be one line, not a hand-built column vector.

#' Select CpG measurement columns by name pattern
#'
#' Returns the CpG-site column names of \code{x} that survive an
#' include/exclude filter on the NAMES. Patterns are fixed substrings
#' (not regular expressions), matched case-sensitively; a column survives
#' \code{include} if it contains ANY of the include patterns, and is then
#' removed by \code{exclude} if it contains ANY of the exclude patterns.
#' Only CpG-site columns (names starting \code{cg...}) are considered -
#' outcomes and covariates are never touched. The same filter is available
#' directly inside \code{dmsa_frame()} as \code{cpgs_include} /
#' \code{cpgs_exclude}.
#'
#' @param x A data.frame, a matrix, or a character vector of column names.
#' @param include Optional character vector of fixed substrings; keep only
#'   CpG columns containing at least one. \code{NULL} keeps all.
#' @param exclude Optional character vector of fixed substrings; drop CpG
#'   columns containing at least one. Applied after \code{include}.
#' @return Character vector of surviving CpG column names.
#' @examples
#' cols <- c("cg00247334_BC21_T4_OXTR",
#'           "cg00247334_BC21_maternal_T1_OXTR",
#'           "cg00247334_BC21_paternal_T1_OXTR",
#'           "EPDS_Total_T4")
#' ## the child's probes: "not the parents'"
#' dmsa_cpg_columns(cols, exclude = c("_maternal", "_paternal"))
#' ## the mother's probes
#' dmsa_cpg_columns(cols, include = "_maternal")
#' @export
dmsa_cpg_columns <- function(x, include = NULL, exclude = NULL) {
  cols <- if (is.character(x) && is.null(dim(x))) x
          else if (is.matrix(x)) colnames(x)
          else names(as.data.frame(x))
  cols <- as.character(cols)
  sites <- cols[.map_is_site(cols)]
  .cpg_name_filter(sites, include, exclude)
}

.cpg_name_filter <- function(sites, include = NULL, exclude = NULL) {
  keep <- sites
  if (!is.null(include)) {
    include <- as.character(include)
    hit <- Reduce(`|`, lapply(include, function(p)
      grepl(p, keep, fixed = TRUE)))
    keep <- keep[hit]
  }
  if (!is.null(exclude) && length(keep)) {
    exclude <- as.character(exclude)
    out <- Reduce(`|`, lapply(exclude, function(p)
      grepl(p, keep, fixed = TRUE)))
    keep <- keep[!out]
  }
  keep
}
