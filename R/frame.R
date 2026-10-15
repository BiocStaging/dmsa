# ============================================================================
# THE TWO-CALL INTERFACE, CALL 1: dmsa_frame()
#
# Declare everything once - data, levels, outcomes, covariates, blocks,
# moderation, map - and get back a validated "dmsa_frame" that dmsa_report()
# runs without further questions. The frame is where every silent data
# problem becomes a visible, documented correction: dmsa_frame() TEST-DRIVES
# the declaration (coercion audit, rank audit, block audit, outcome audit,
# B = 49 pilot) and records every action in frame$corrections.
#
# Moderation layout (Tsachi, 14 Aug 2026):
#   moderation = FALSE by default; mod = ""; mod2 = "".
#   moderation = TRUE + mod            -> tested term is frame x mod
#   moderation = TRUE + mod + mod2     -> frame x mod x mod2 (moderators
#     mean-centered, all lower-order terms in the model, the highest-order
#     term is the test). Interactions ride the COMPOSITE lens - the one lens
#     that carries to products - so moderation runs are composite-only and say
#     so in the report.
#
# cpg_map: "full" (DEFAULT) keeps every CpG-gene pair that carries a usable
# direction call; "confidence" is an OPT-IN sensitivity that keeps only the
# pairs meeting the high-confidence bar (tier A or SMR S1). "full" is the
# default because a usable direction is what DMSA needs to align a probe -
# the confidence bar is a stricter subset of that, not a different
# definition of membership, and it must never change which genes belong to
# a system. BOTH alignments are built regardless; any gene or
# system whose pooled direction flips between the two maps is analysed under
# the chosen map but warned about, listed in frame$map_conflicts, and repeated
# in the report's Design notes.
# ============================================================================

## The four Project Alpha builds each document a DIFFERENT covariate contract in
## their own Covariates sheet. A single hard-coded set silently degraded three of
## the four to sex alone - and on the before/after builds that drops `time`,
## which is the entire design. Detect the build from the columns present and
## return its documented contract; refuse to guess when nothing matches.
.frame_contracts <- list(
  parent_T1 = list(
    need = c("age_at_array_T1", "Epi_T1", "Fib_T1", "ctrlSV3_T1", "ctrlSV5_T1"),
    fixed = c("sex_c", "age_at_array_T1", "Epi_T1", "Fib_T1",
              "ctrlSV3_T1", "ctrlSV5_T1"),
    block = c("cID"), chip = "chip_T1",
    label = "build 1: parent T1"),
  child_T4 = list(
    need = c("birth_week", "birth_weight", "Epi_T4"),
    fixed = c("sex_c", "birth_week", "birth_weight", "Epi_T4"),
    block = c("cID"), chip = "chip_T4",
    label = "build 4: child T4 (leanest set - n is small)"),
  before_after = list(
    need = c("time", "Epi_array", "Fib_array", "ctrlSV5_array"),
    fixed = c("time", "sex_c", "Epi_array", "Fib_array", "ctrlSV5_array"),
    block = c("ID", "cID"), chip = NULL,
    label = "build 2/3: before-after (repeated measures; chip is NEVER included)")
)

.frame_which_build <- function(nms) {
  hit <- vapply(.frame_contracts, function(z) all(z$need %in% nms), logical(1))
  if (!any(hit)) return(NULL)
  names(.frame_contracts)[which(hit)[1L]]
}

.frame_contract <- function(nms = NULL) {
  if (is.null(nms)) return(.frame_contracts$parent_T1$fixed)
  b <- .frame_which_build(nms)
  if (is.null(b))
    stop("covariates = \"contract\" could not identify which Project Alpha ",
         "build this data frame is: none of the documented contracts had its ",
         "columns present. Expected one of - parent T1 (age_at_array_T1, ",
         "Epi_T1, Fib_T1, ctrlSV3_T1, ctrlSV5_T1), child T4 (birth_week, ",
         "birth_weight, Epi_T4), or before-after (time, Epi_array, Fib_array, ",
         "ctrlSV5_array). Pass the covariates explicitly instead of ",
         "\"contract\".", call. = FALSE)
  .frame_contracts[[b]]$fixed
}

## BUILD-AWARE COLUMN NAMES. The same CpG is named differently in every
## Project Alpha workbook: cg..._TC21_T1_AVP in the parent build,
## cg..._TC21_AVP in the two before-after builds, and cg..._TC21_T4_AVP plus
## cg..._TC21_maternal_T1_AVP / _paternal_T1_AVP in the child build. The
## bundled direction map carries the parent-T1 name only, so before 1.16.0 a
## child or before-after frame matched ZERO probes and stopped with "no mapped
## methylation columns found" - three of the four builds could not be analysed
## at all. The selection cascade lists every variant explicitly; use it.
.frame_col_variants <- function(cs) {
  if (is.null(cs) || is.null(cs$cascade)) return(NULL)
  cd <- cs$cascade
  vs <- grep("^col_", names(cd), value = TRUE)
  if (!length(vs) || !"col_parent_T1" %in% vs) return(NULL)
  cd[, vs, drop = FALSE]
}

## Returns the column vector to use plus a note when it had to be translated.
## `build` is the contract-detected build; it breaks ties, because the child
## workbook carries THREE matching sets (the child's own T4 array and both
## parents' T1 arrays) and picking by match count alone would be arbitrary.
.frame_resolve_columns <- function(mp, avail, cs, build = NULL) {
  keep0 <- list(column = mp$column, note = NULL)
  tab <- .frame_col_variants(cs)
  if (is.null(tab) || !length(avail)) return(keep0)
  hit0 <- sum(mp$column %in% avail)
  i   <- match(mp$column, tab$col_parent_T1)
  if (all(is.na(i))) return(keep0)
  vs  <- names(tab)
  cnt <- vapply(vs, function(v) {
    x <- tab[[v]][i]; sum(!is.na(x) & x %in% avail)
  }, integer(1))
  if (!any(cnt > 0L)) return(keep0)
  best <- vs[which.max(cnt)]
  pref <- c(parent_T1 = "col_parent_T1", child_T4 = "col_child_T4",
            before_after = "col_long")
  if (!is.null(build) && build %in% names(pref) && pref[[build]] %in% vs &&
      cnt[[pref[[build]]]] > 0L) best <- pref[[build]]
  if (cnt[[best]] <= hit0) return(keep0)
  new  <- tab[[best]][i]
  bad  <- is.na(new)
  new[bad] <- mp$column[bad]
  other <- vs[vs != best & cnt[vs] > 0L]
  list(column = new, note = sprintf(
    paste0("the direction map names probes in the parent-T1 form, this data ",
           "frame uses the '%s' form -> translated, %d probe(s) matched%s"),
    sub("^col_", "", best), cnt[[best]],
    if (length(other))
      paste0(" (this workbook also carries: ",
             paste(sub("^col_", "", other), collapse = ", "),
             " - pass `methylation` as the column names you want to override ",
             "the default)")
    else ""))
}

## SHARED PROBES. 229 probes in the Alpha panel sit inside two genes and
## appear as a column under each name (IGF2+INS, GABRA5+GABRB3,
## HLA-DPA1+HLA-DPB1, PEG3+ZIM2 and others). Two such genes in the same
## level-local family are NOT two independent units - they are the same
## measurement twice - and if both survive they read as two findings. Say so
## at frame time, while the design can still be changed.
.frame_shared_probes <- function(map, cs = NULL) {
  need <- c("probe", "gene", "system_id")
  if (!all(need %in% names(map))) return(NULL)
  cd <- if (is.null(cs)) NULL else cs$cascade
  ## The cascade states sharing EXPLICITLY in `shared_with`. Inferring it from
  ## duplicated probe names in the analysed map does not work: an Alpha column
  ## name embeds its gene, so one probe under two genes is two differently
  ## named columns, and after the map is de-duplicated on `column` the
  ## duplication is invisible.
  if (is.null(cd) || !"shared_with" %in% names(cd)) return(NULL)
  sw <- ifelse(is.na(cd$shared_with), "", as.character(cd$shared_with))
  ok <- nzchar(sw)
  if (!any(ok)) return(NULL)
  cd <- cd[ok, , drop = FALSE]; sw <- sw[ok]
  pk <- if ("cpg" %in% names(cd)) cd$cpg else cd$probe_id
  if (is.null(pk)) return(NULL)

  ## "IGF2; INS" is two partners, not one gene called "IGF2; INS"
  parts <- strsplit(sw, ";", fixed = TRUE)
  n <- lengths(parts)
  long <- data.frame(sid   = rep(as.character(cd$system_id), n),
                     probe = rep(as.character(pk), n),
                     g1    = rep(as.character(cd$gene), n),
                     g2    = trimws(unlist(parts, use.names = FALSE)),
                     stringsAsFactors = FALSE)
  long <- long[nzchar(long$g2) & long$g1 != long$g2, , drop = FALSE]
  if (!nrow(long)) return(NULL)

  ## Only sharing that actually reached the analysis counts: the probe must have
  ## entered under BOTH gene names, in a system that was selected.
  ana  <- unique(data.frame(sid   = as.character(map$system_id),
                            gene  = as.character(map$gene),
                            probe = as.character(map$probe),
                            stringsAsFactors = FALSE))
  keyp <- paste(ana$sid, ana$gene, ana$probe)
  long <- long[paste(long$sid, long$g1, long$probe) %in% keyp &
               paste(long$sid, long$g2, long$probe) %in% keyp, , drop = FALSE]
  if (!nrow(long)) return(NULL)

  syl  <- if ("system" %in% names(map))
            stats::setNames(as.character(map$system),
                            as.character(map$system_id)) else NULL
  pair <- paste0(ifelse(is.null(syl), long$sid, syl[long$sid]), "|",
                 pmin(long$g1, long$g2), " + ", pmax(long$g1, long$g2))
  u <- unique(data.frame(pair = pair, probe = long$probe,
                         stringsAsFactors = FALSE))
  tb <- table(u$pair)
  out <- data.frame(pair = names(tb), n_shared = as.integer(tb),
                    stringsAsFactors = FALSE)
  out[order(-out$n_shared, out$pair), , drop = FALSE]
}



## OUTCOME KIND. The Results prose used one template - "Higher X was associated
## with ..." - which assumes a CONTINUOUS outcome. It produced "Higher sex was
## associated with ..." and "Higher time was associated with ...", both
## meaningless. Three cases have to be told apart, and by the DATA, not by the
## column name:
##   continuous - the template is correct
##   binary     - a between-group contrast (sex)
##   wave       - a WITHIN-PERSON change (time, in the before/after builds),
##                which is what makes those builds worth having
## `wave` requires two things: a before/after contract, and an outcome that
## actually varies inside the finest permutation block. Testing only "binary +
## has blocks" would misclassify sex in the parent build, where a family block
## (cID) legitimately contains both parents.
## Spec 50: resolve `outcome_type` into one declared family per outcome.
## Accepts a single value (applies to every outcome) or a vector naming every
## outcome. Unnamed multi-element input is refused: positional matching of
## families to outcomes is exactly the kind of silent mis-pairing this fixes.
.frame_otypes <- function(value, outcome) {
  ok <- c("gaussian", "logistic", "multinomial")
  ## the caller resolves the missing-argument default BEFORE this function, so
  ## a literal c("gaussian","logistic","multinomial") arriving here is a real
  ## (unnamed, positional) user input and falls through to the refusal below.
  if (!is.character(value) || !length(value) || anyNA(value))
    stop("`outcome_type` must be a character vector naming one family per ",
         "outcome. Allowed values: ", paste(ok, collapse = ", "), ".",
         call. = FALSE)
  bad <- setdiff(value, ok)
  if (length(bad))
    stop("`outcome_type` value(s) not recognised: ",
         paste(unique(bad), collapse = ", "), ".\nAllowed: ",
         paste(ok, collapse = ", "), ".", call. = FALSE)
  nm <- names(value)
  if (length(value) == 1L && (is.null(nm) || !nzchar(nm)))
    return(stats::setNames(rep(value, length(outcome)), outcome))
  if (is.null(nm) || !all(nzchar(nm)))
    stop("`outcome_type` has ", length(value), " values for ", length(outcome),
         " outcome(s) but does not name them.\nEach family must be attached to ",
         "its outcome by name, e.g.\n  outcome_type = c(",
         paste(sprintf("%s = \"gaussian\"", outcome), collapse = ", "), ")\n",
         "Positional matching is not accepted here: a family silently paired ",
         "with the wrong outcome is not detectable in the output.",
         call. = FALSE)
  extra <- setdiff(nm, outcome)
  if (length(extra))
    stop("`outcome_type` names column(s) that are not outcomes: ",
         paste(extra, collapse = ", "), ".\nOutcomes: ",
         paste(outcome, collapse = ", "), ".", call. = FALSE)
  gap <- setdiff(outcome, nm)
  if (length(gap))
    stop("`outcome_type` does not declare a family for outcome(s): ",
         paste(gap, collapse = ", "), ".\nWhen `outcome_type` is named it must ",
         "name every outcome - an undeclared one would fall back to gaussian ",
         "without saying so.", call. = FALSE)
  if (anyDuplicated(nm))
    stop("`outcome_type` names an outcome more than once: ",
         paste(unique(nm[duplicated(nm)]), collapse = ", "), call. = FALSE)
  value[outcome]
}

.frame_outcome_kind <- function(base, oc, blockv, build) {
  v <- base[[oc]]
  ## a two-level outcome is binary WHATEVER its storage type: a 0/1 numeric,
  ## a logical, a factor, or a character column all describe the same group
  ## contrast, and the report's direction sentence must say so for all of
  ## them (a factor-coded pills column used to fall through to "continuous"
  ## and produce "Higher pills_past_T1 was associated ..." for a 0/1
  ## variable - PI, 2026-08-29)
  u <- if (is.numeric(v)) sort(unique(v[is.finite(v)])) else
    sort(unique(as.character(v[!is.na(v) & as.character(v) != ""])))
  if (length(u) != 2L) return(list(kind = "continuous", lo = NA, hi = NA))
  within <- FALSE
  if (identical(build, "before_after") && length(blockv)) {
    id <- blockv[[1]]
    if (id %in% names(base)) {
      k <- tapply(v, base[[id]], function(z) length(unique(z[is.finite(z)])))
      within <- isTRUE(mean(k > 1, na.rm = TRUE) > 0.5)
    }
  }
  list(kind = if (within) "wave" else "binary", lo = u[1], hi = u[2])
}

## non-empty and not NA - `chip` can legitimately be NA (the before/after
## contract forbids chip entirely), and nzchar(NA) is NA, not FALSE.
.nz <- function(x) !is.null(x) && length(x) >= 1L && !is.na(x[1]) && nzchar(x[1])

.frame_note <- function(cor, field, issue, action)
  rbind(cor, data.frame(field = field, issue = issue, action = action,
                        stringsAsFactors = FALSE))

## ---- the bundled Alpha map -------------------------------------------------
.frame_read_map <- function(map = "alpha") {
  if (is.data.frame(map)) return(as.data.frame(map))
  f <- system.file("extdata", "coverage_v4_full.csv", package = "dmsa")
  if (!nzchar(f) || !file.exists(f))
    stop("bundled Alpha map not found; pass `map` as a data.frame with ",
         "columns probe, column, gene, system_id, system, best_direction, ",
         "p_plus (and optionally best_tier, smr_tier, module)", call. = FALSE)
  utils::read.csv(f, stringsAsFactors = FALSE)
}

.frame_maps <- function(map = "alpha") {
  mp <- .frame_read_map(map)
  need <- c("probe", "column", "gene", "system_id", "system",
            "best_direction", "p_plus")
  miss <- setdiff(need, names(mp))
  if (length(miss))
    stop("map is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  ## spec 11: NO dedup by column. One physical column mapping to several
  ## genes is the design (co-effects), not a duplicate; the bundled snapshot
  ## happens to be one row per column, so for it this changes nothing, but a
  ## user map may legitimately carry pairs.
  usable <- is.finite(mp$best_direction)
  hi <- usable
  if (!is.null(mp$best_tier) || !is.null(mp$smr_tier)) {
    ta <- if (!is.null(mp$best_tier)) mp$best_tier %in% "A" else FALSE
    s1 <- if (!is.null(mp$smr_tier))  mp$smr_tier  %in% "S1" else FALSE
    hi <- usable & (ta | s1)
  }
  list(full = mp[usable, , drop = FALSE],
       confidence = mp[hi, , drop = FALSE])
}

.frame_modules <- function(sets = "alpha") {
  ## The selection cascade is the authority on module membership when one is
  ## available: it is the same file the user selects systems from, so the family
  ## a unit is corrected within cannot drift from the family it was declared in.
  ## modules_alpha.csv is kept only as a fallback for installs without a cascade.
  cs <- .frame_sets(sets)
  if (!is.null(cs)) {
    mm <- unique(cs$cascade[, c("system_id", "module_id", "module", "gene")])
    return(mm[order(.cas_num(mm$system_id), .cas_modnum(mm$module_id),
                    mm$gene), , drop = FALSE])
  }
  f <- system.file("extdata", "modules_alpha.csv", package = "dmsa")
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  utils::read.csv(f, stringsAsFactors = FALSE)
}

## Load a cascade for use inside dmsa_frame(), tolerantly: a frame must still be
## constructible when no cascade is installed or a user passes cascade = NULL.
.frame_sets <- function(sets = "alpha") {
  if (is.null(sets) || isFALSE(sets)) return(NULL)
  if (inherits(sets, "dmsa_sets")) return(sets)
  if (inherits(sets, "dmsa_selection"))
    return(structure(list(cascade = sets$cascade,
                          systems = sets$systems,
                          modules = sets$modules,
                          source = sets$source, name = sets$name,
                          columns = grep("^col_", names(sets$cascade),
                                         value = TRUE)),
                     class = "dmsa_sets"))
  out <- try(dmsa_sets(sets), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

## `systems =` accepts the cascade's short names ("hpa"), ids, or full names.
## Resolution runs through the cascade first so that the handle a user typed is
## the handle the cascade documents; the old exact-match behaviour is the
## fallback, so existing scripts keep working.
.frame_pick_systems <- function(systems, sys_tab, cs) {
  if (is.null(systems)) return(sys_tab)
  if (!is.null(cs)) {
    ## Only route through the cascade when the direction map is describing the
    ## same panel: system NAMES must overlap. Otherwise a numeric `systems = 2`
    ## against an unrelated map would silently pick up the cascade's system 2.
    same_panel <- any(sys_tab$system %in% cs$systems$system)
    got <- if (same_panel) try(.cas_resolve_systems(cs, systems), silent = TRUE)
           else structure("skip", class = "try-error")
    if (!inherits(got, "try-error") && nrow(got)) {
      keep <- sys_tab$system %in% got$system
      if (any(keep)) {
        out <- sys_tab[keep, , drop = FALSE]
        attr(out, "short") <- got$system_short[match(out$system, got$system)]
        attr(out, "absent") <- got$system_short[!got$system %in% out$system]
        return(out)
      }
      stop("`systems` resolved to ", paste(got$system_short, collapse = ", "),
           " but no mapped probe of those system(s) is present in the data",
           call. = FALSE)
    }
  }
  keep <- if (is.numeric(systems)) sys_tab$system_id %in% systems else
    sys_tab$system %in% systems |
    as.character(sys_tab$system_id) %in% as.character(systems)
  if (!any(keep))
    stop("`systems` matched nothing; available: ",
         paste(sprintf("%s '%s'", sys_tab$system_id, sys_tab$system),
               collapse = "; "),
         if (!is.null(cs)) paste0("\nshort names: ",
           paste(cs$systems$system_short, collapse = ", ")) else "",
         call. = FALSE)
  sys_tab[keep, , drop = FALSE]
}

## pooled direction of a unit under one map: sign of sum d * (2 p_plus - 1)
.frame_pooled_sign <- function(mp, by) {
  m <- mp$best_direction * (2 * pmin(pmax(mp$p_plus, 0), 1) - 1)
  m[!is.finite(m)] <- 0
  s <- tapply(m, mp[[by]], function(v) sign(sum(v)))
  s[is.na(s)] <- 0
  s
}

.frame_conflicts <- function(maps, systems_keep) {
  out <- data.frame(level = character(), unit = character(),
                    sign_confidence = integer(), sign_full = integer(),
                    stringsAsFactors = FALSE)
  for (lv in c("gene", "system")) {
    a <- .frame_pooled_sign(maps$confidence[maps$confidence$system_id %in%
                                              systems_keep, ], lv)
    b <- .frame_pooled_sign(maps$full[maps$full$system_id %in%
                                        systems_keep, ], lv)
    u <- intersect(names(a), names(b))
    bad <- u[a[u] != 0 & b[u] != 0 & a[u] != b[u]]
    if (length(bad))
      out <- rbind(out, data.frame(level = lv, unit = bad,
                                   sign_confidence = as.integer(a[bad]),
                                   sign_full = as.integer(b[bad]),
                                   stringsAsFactors = FALSE))
  }
  out
}

## ---- beta -> M when needed -------------------------------------------------
.frame_mvalues <- function(BE, eps = 1e-4) {
  rng <- range(BE, na.rm = TRUE)
  if (rng[1] >= 0 && rng[2] <= 1) {
    BE <- pmin(pmax(BE, eps), 1 - eps)
    list(M = log2(BE / (1 - BE)), converted = TRUE)
  } else list(M = BE, converted = FALSE)
}

#' Declare a DMSA analysis frame
#'
#' Call 1 of the two-call interface. Declares the data, the levels
#' (system / module / gene / probe), the outcomes, the covariates and blocks,
#' optional moderation, and the direction map - validates all of it, runs the
#' test drive (with \code{autofix = TRUE} correcting what is safely
#' correctable, LOUDLY), and returns a \code{dmsa_frame} for
#' \code{\link{dmsa_report}}.
#'
#' @param data data.frame holding outcomes, covariates, ids and (by default)
#'   the Alpha methylation columns.
#' @param methylation NULL (auto-detect the bundled map's columns in
#'   \code{data}), a character vector of column names, or a numeric matrix
#'   whose colnames are probe/column ids.
#' @param map "alpha" for the bundled Project Alpha map, or a data.frame.
#'   Used by the BUNDLED direction source only; passing a table forces
#'   \code{direction_source = "bundled"}.
#' @param reference The biological reference (system > module > gene) that
#'   defines membership on the pair path: \code{"alpha"} (default, the 2026c
#'   codebook - 30 systems, 187 modules, 1,282 genes) or a
#'   \code{dmsa_reference} object. The reference never selects CpGs: which
#'   CpGs exist is read from the data, which genes they support from
#'   cpgdirection, and the reference then says which of those genes belong to
#'   the biology being asked about.
#' @param direction_source Where CpG-to-gene direction comes from.
#'   \code{"auto"} (default) uses \code{cpgdirection::cpg_gene_pairs()} -
#'   one row per CpG x gene pair, direction resolved within the pair,
#'   pair-specific abstention - when cpgdirection is installed, and otherwise
#'   falls back to the bundled per-column Alpha snapshot with a recorded
#'   note. \code{"cpgdirection"} requires the package; \code{"bundled"}
#'   forces the snapshot. On the pair path the frame keeps the full pair
#'   ledger (\code{$pair_ledger}; see \code{dmsa_coverage()}). For testing,
#'   \code{options(dmsa.pair_table = <pair data.frame>)} injects a pair table
#'   in place of the live \code{cpg_gene_pairs()} call - the regression suite
#'   uses this seam so the pair path is exercised on machines where
#'   cpgdirection is not installed.
#' @param tissue Tissue context for pair discovery: \code{"blood"} (default),
#'   \code{"nasal_epithelium"}, or \code{"solid_tissue"}. Pair path only.
#' @param replicate_policy \code{"collapse_site"} (default): several probe
#'   designs of one CpG are ONE measurement, averaged on the beta scale
#'   before analysis (a missing value in any replicate keeps the site visible
#'   to \code{missing_methylation} rather than quietly narrowing its basis);
#'   \code{"keep_correlated"} keeps every design as its own measurement.
#'   Pair path only.
#' @param cpgs_include,cpgs_exclude Optional CpG-set selection by NAME
#'   pattern, for methylation files that carry several biological groups in
#'   one table (in Alpha's child file, \code{"_maternal"} and
#'   \code{"_paternal"} mark the parents' probes; the child's carry no
#'   tag). Fixed substrings, not regular expressions: a CpG column survives
#'   \code{cpgs_include} if its name contains ANY of the patterns
#'   (\code{NULL} = all), and is then removed by \code{cpgs_exclude} if it
#'   contains any of those. Analysing the child is one line:
#'   \code{cpgs_exclude = c("_maternal", "_paternal")}. Only CpG-site
#'   columns are filtered - outcomes and covariates are never touched - and
#'   the drop is recorded in the corrections table. See also
#'   \code{dmsa_cpg_columns()} for the same filter as a standalone helper.
#' @param system,module,gene,probe logical: which levels to analyse.
#' @param sample_id Optional column of \code{data} holding sample
#'   identifiers. When the methylation matrix carries rownames, rows are
#'   matched by identifier and reordered to match \code{data}; any
#'   mismatch or duplication is an error rather than a silent
#'   realignment.
#' @param sample_align "id" (default) matches a methylation matrix to
#'   \code{data} by sample identifier where identifiers exist;
#'   "positional" declares that the rows are already in the same order.
#'   With no identifiers on either side a shuffle is undetectable, so
#'   alignment falls back to positional and is recorded in the
#'   corrections table.
#' @param missing_methylation What to do when the methylation block contains
#'   missing or non-finite values. \code{"error"} (default) stops and reports
#'   how many probes and samples are affected; \code{"drop_probes"} removes the
#'   affected probes and records it; \code{"common_complete_rows"} keeps every
#'   probe and analyses the samples complete across all of them. There is no
#'   silent option: downstream a probe with any non-finite value is dropped
#'   from its set, and that must be a declared choice.
#' @param methylation_scale \code{"auto"} (default) treats values inside
#'   [0, 1] as beta and converts them to M-values; \code{"beta"} declares them
#'   beta, and values outside [0, 1] are then an error rather than being
#'   clamped; \code{"M"} declares M-values and converts nothing.
#' @param systems Optional subset of systems to analyse. Accepts the selection
#'   cascade's short names - \code{systems = c("hpa", "oxytocin")} - and also
#'   system ids or full names. Matching is case-insensitive and a unique prefix
#'   of a short name is enough. \code{NULL} (default) analyses every covered
#'   system. Everything below the system - modules, genes, probes - is taken in
#'   full; use \code{dmsa_select()} and pass the result as \code{sets} to
#'   narrow a lower level. Run \code{dmsa_systems()} to print the short names.
#' @param sets The selection cascade declaring system > module > gene > probe.
#'   \code{"alpha"} (default) uses the bundled, module-audited Project Alpha
#'   2026c cascade; a path, \code{data.frame}, \code{dmsa_sets} or
#'   \code{dmsa_selection} may be supplied instead, and \code{NULL} disables it
#'   (module membership then falls back to the legacy bundled module map). See
#'   \code{dmsa_sets_template()} for the schema.
#' @param ... Not for new arguments. The former \code{outcome} argument is
#'   accepted here so existing scripts keep working (it means exactly
#'   \code{outcomes}); anything else is an error. A side effect of this slot:
#'   every argument after \code{data} and \code{methylation} must be given
#'   by its full name.
#' @param outcomes Which column(s) of \code{data} the methylation is tested
#'   against - the phenotype side of the model. With several, each is tested
#'   in turn with the others entering as covariates (mutual adjustment). If
#'   your question runs the other way - these columns PREDICT the
#'   methylation - name them with \code{predictors =} instead and skip
#'   \code{frame_role} entirely.
#' @param predictors The natural spelling when METHYLATION IS YOUR OUTCOME:
#'   which column(s) of \code{data} predict the methylation frame. Passing
#'   \code{predictors =} sets \code{frame_role = "outcome"} for you (an
#'   explicit \code{frame_role = "predictor"} alongside it is a
#'   contradiction and errors). Same columns as \code{outcome}, opposite
#'   reading; give exactly ONE of outcome / outcomes / predictors.
#' @param covariates "contract" for the Alpha student contract, or a character
#'   vector. chip is added as a fixed factor automatically when present.
#' @param blocks The exchangeability blocks: the column(s) of \code{data}
#'   whose rows must travel TOGETHER under permutation because they are not
#'   independent - the participant id when a person contributes more than one
#'   row, the couple id for partners, the family id for relatives. There is
#'   NO default (PI ruling, 2026-09-02): if you declare none, every row is
#'   permuted independently and \code{dmsa_frame()} says so and asks you to
#'   confirm that is intended; a block column you DID name but that is absent
#'   from \code{data} is a hard error (spec 43), because dropping it would
#'   silently change the null model. \code{covariates = "contract"} carries
#'   the Alpha build's block structure and is adopted when you give none.
#'   TWO factors and \code{dmsa_frame()} works out the structure itself and
#'   says so in lme4 notation: NESTED (e.g. couples on the same chip,
#'   \code{(1 | chip/cID)}) blocks on the nested pair, with the coarse
#'   factor's variance component left to \code{chip =}; CROSSED (partners of
#'   one couple on DIFFERENT chips, \code{(1 | cID) + (1 | chip)}) routes the
#'   finest dependence factor to the blocks and the other into the
#'   \code{(1 | .)} random intercept. Aliased factors are detected; three or
#'   more are refused (the engine fits one block structure and one
#'   intercept).
#' @param random_effects Deprecated spelling of \code{blocks} (it always named
#'   the exchangeability blocks, not lme4 random effects). Still accepted,
#'   with a message; giving both is an error.
#' @param chip how to handle the array/chip batch factor. \code{TRUE} (default)
#'   enters it as a fixed effect, using \code{chip_f} if present and otherwise
#'   building it from a \code{chip}/\code{chip_T*} column. \code{FALSE} leaves
#'   it out entirely - on the Alpha parent build the factor has 61 levels on
#'   about 400 arrays, so it costs 60 degrees of freedom that the control-probe
#'   surrogate variables already partly absorb, and 7 of those 61 chips hold a
#'   single array - fitted exactly by their own dummy (leverage 1, residual 0),
#'   so they are counted in \code{n} while contributing nothing to any estimate.
#'   \code{"pool"} keeps the batch adjustment but merges every level holding
#'   fewer than 3 arrays into one \code{chip_small} level, which recovers the
#'   degrees of freedom and removes the leverage-1 rows. Any other character
#'   string names the column to use. The before/after contracts forbid chip and
#'   override this argument.
#' @param chip_effect how the chip factor enters the model, once \code{chip}
#'   has named it. \code{"random"} (default, and what the Alpha covariate
#'   contract specifies) takes chip out of the fixed design and fits it as a
#'   one-way random intercept, \code{(1 | chip)}, by REML quasi-demeaning: a
#'   median of about 16 effective degrees of freedom on the parent build, and
#'   zero on the probes where the chip variance is estimated as nil.
#'   \code{"fixed"} restores the pre-1.18.0 behaviour and enters chip as a
#'   fixed factor, at one degree of freedom per chip beyond the first.
#'   \code{"none"} drops the batch adjustment altogether, whatever
#'   \code{chip} says. Validity does not rest on the variance component being
#'   right: the permutation null is computed under the same transform, so the
#'   random-intercept fit moves power, not type-I error.
#' @param outcome_levels Optional labels for the two levels of a two-level
#'   outcome, used only in the written report so a sentence can say
#'   "taking pills [1] ... than not taking pills [0]" instead of falling
#'   back to \code{outcome = 0} / \code{outcome = 1}. Either a length-2
#'   character vector - honoured when the frame has ONE outcome - or a
#'   named list with one such pair per outcome, e.g.
#'   \code{list(pills_past_T1 = c("not taking pills", "taking pills"))},
#'   which works for any number of outcomes. Each pair is safest given as
#'   a NAMED vector matched by coded value, e.g. \code{c("0" = "not
#'   taking pills", "1" = "taking pills")} - immune to flipping whatever
#'   the coding (1/2, 0/1, "a"/"b"); an unnamed pair is taken in
#'   sorted-value order (low, high), and a name that is not one of the
#'   variable's two actual levels is an immediate error naming them. \code{NULL} (default) uses the coded-value fallback.
#'   Purely cosmetic: no level name reaches the model, the permutation or
#'   the family-wise correction. Spell it the way you declared the
#'   variable: with \code{outcomes = } use this argument; with
#'   \code{predictors = } use \code{predictor_levels} - same meaning,
#'   your vocabulary.
#' @param predictor_levels The same level labels, for a frame declared with
#'   \code{predictors = } (\code{frame_role = "outcome"}): someone who
#'   wrote \code{predictors = "pills_past_T1"} thinks of pills as a
#'   predictor, so the labels are asked for in those words. Exactly one of
#'   \code{predictor_levels}/\code{outcome_levels} may be given, and each
#'   requires its matching declaration.
#' @param moderation,mod,mod2 moderation switch and moderator column(s); see
#'   the layout note above.
#' @param outcome_label,predictor_labels,covariate_labels,mod_label,mod2_label
#'   Optional display NAMES used in figures and the written report instead of
#'   raw column names, spelled the way you declared each variable:
#'   \code{predictor_labels} names the \code{predictors = } column(s),
#'   \code{outcome_label} the \code{outcomes = } column(s) (the two are the
#'   same column - give one, matching your declaration),
#'   \code{covariate_labels} the covariates, and \code{mod_label}/
#'   \code{mod2_label} are ONE name each for the moderator column(s) (a
#'   moderator's LEVEL labels go in \code{mod_levels}). Either one per
#'   column, in order, or a named vector keyed by column name. Purely
#'   cosmetic: no label reaches the model, the permutation or the
#'   family-wise correction.
#' @param mod_levels Labels for the two LEVELS of a two-level moderator,
#'   used in the moderation figures and prose. Safest as a named vector
#'   matched by coded value - \code{mod_levels = c("1" = "Husband",
#'   "2" = "Wife")} - which cannot flip whatever the coding; an unnamed
#'   pair is taken in sorted-value order (low, high).
#' @param type "linear", "non-linear" (linear/quadratic/threshold arms,
#'   ACAT-combined), or "exponential" (adds an exp arm).
#' @param outcome_type The declared family of each outcome: "gaussian",
#'   "logistic", or "multinomial". A single value applies to every outcome. With
#'   several outcomes of different kinds, name them -
#'   \code{outcome_type = c(epds = "gaussian", depressed = "logistic")} - and
#'   every outcome must be named. An unnamed vector of several values is
#'   refused: a family paired positionally with the wrong outcome is not
#'   visible anywhere in the output. The stored value is the bare string for a
#'   single-outcome frame and a vector named by outcome otherwise.
#' @param frame_role Which way the SECONDARY models face. Most callers never
#'   need this argument: passing the phenotype columns as \code{predictors =}
#'   sets it for you. The main DMSA path is
#'   always \code{methylation ~ outcome + covariates}: a directional set
#'   statistic on methylation is what DMSA estimates, and no argument turns that
#'   around. \code{"predictor"} (default) faces the moderated model the same way
#'   as the shape scan, \code{outcome ~ tone score x moderator};
#'   \code{"outcome"} declares that the named columns predict the frame, giving
#'   \code{tone score ~ outcome x moderator}. Because the shape scan cannot be
#'   turned around, \code{frame_role = "outcome"} requires
#'   \code{type = "linear"}. The resolved orientation of every path is stored on
#'   the frame as \code{model_orientation} and printed.
#' @param cpg_map "full" (default) or "confidence"; see the note above.
#'   "full" analyses every CpG-gene pair with a usable direction call.
#'   "confidence" is an opt-in sensitivity restricted to high-confidence
#'   calls (tier A or SMR S1); it is a strict subset of "full" and never
#'   redefines system, module or gene membership.
#' @param B,alpha,seed,correction inference settings (correction is applied
#'   within level-local families only).
#' @param fdr Level of the naming rule (PI ruling 2026-10-15): a unit is named
#'   when its dependence-calibrated Benjamini-Hochberg q-value
#'   (\code{q_unit_dep}, see \code{\link{dmsa_triangulate}}) is at or below
#'   \code{fdr} within its level-local family. Default .10. \code{alpha}
#'   remains the level of the family-wise badges and of the per-lens
#'   pre-registered any-lens rule, which is reported beside the naming rule.
#' @param progress Show a 0-100% progress bar while \code{dmsa_report()} runs.
#'   Defaults to \code{interactive()}: a bar writes carriage returns, which make
#'   a piped log unreadable, so it is off under \code{Rscript} and in tests.
#' @param beep Sound a completion signal when \code{dmsa_report()} finishes.
#'   \code{TRUE} uses \pkg{beepr} sound 8; a number picks another sound;
#'   \code{FALSE} is silent. Defaults to \code{interactive()}, so a run never
#'   makes noise during \code{R CMD check} or a testthat run. \pkg{beepr} is a
#'   Suggests: without it this is silently a no-op.
#' @param gene_models Draw real exon/intron structure (with the TSS) under
#'   the probes in each locus panel. The default \code{"auto"} fetches the
#'   model from Ensembl for NAMED genes only - a handful of small requests,
#'   announced in the report log - and degrades to the bare coordinate axis
#'   with a message when offline. \code{TRUE} forces the fetch,
#'   \code{FALSE} never goes to the network, and a table built once with
#'   \code{\link{dmsa_gene_model}} can be passed to stay offline with the
#'   exons still drawn.
#' @param palette Colour palette for every figure the report draws. One of
#'   the \pkg{viridisLite} palettes - \code{"viridis"} (default; the
#'   colour-blind-safe standard), \code{"magma"}, \code{"plasma"},
#'   \code{"inferno"}, \code{"cividis"} (designed for deuteranopia),
#'   \code{"mako"}, \code{"rocket"}, or \code{"turbo"} (vivid, not
#'   colour-blind-safe - avoid for publication). All are perceptually
#'   uniform and print safely in greyscale except \code{"turbo"}. If
#'   \pkg{viridisLite} is not installed, any value falls back to base R's
#'   \code{hcl.colors(..., "viridis")}, so the figures always draw. The
#'   per-system accent colours of the overview panels come from a separate
#'   fixed scheme and are not affected.
#' @param plots,tables,summary,plot_type,table_type,outdir output
#'   styling for \code{dmsa_report}.
#' @param autofix TRUE: the test drive corrects what is safely correctable and
#' @param weighting Character. Probe weighting engine within a unit: \code{"combined"} (default) fuses the flat and reliability statistics on one shared permutation stream, \code{"flat"} weights every usable aligned probe equally, \code{"reliability"} weights each probe by its item-rest correlation with the rest of its unit. Weights are computed from methylation alone, so the permutation null is unaffected.
#' @param w_floor Numeric. Lower bound applied to reliability weights before normalisation, so a single poorly behaved probe cannot be driven to zero influence. Ignored when \code{weighting = "flat"}.
#'   records it; FALSE: strict errors instead.
#' @return object of class \code{dmsa_frame}.
#' @examples
#' set.seed(42)
#' map <- data.frame(gene = rep(paste0("G", 1:3), each = 3), system_id = 1L,
#'                   system = "Sim system")
#' map$probe <- sprintf("cg%07d", seq_len(nrow(map)))
#' map$column <- paste0(map$probe, "_", map$gene)
#' map$best_direction <- rep(c(-1, 1), length.out = nrow(map))
#' map$p_plus <- ifelse(map$best_direction > 0, 0.9, 0.1)
#' d <- data.frame(out1 = rnorm(60), cov1 = rnorm(60),
#'                 cID = rep(1:30, each = 2))
#' sig <- 0.5 * outer(d$out1, map$best_direction * (map$gene == "G1"))
#' d[map$column] <- plogis(matrix(rnorm(60 * nrow(map)), 60) + sig)
#' dmsa_frame(d, map = map, outcomes = "out1", covariates = "cov1",
#'            blocks = "cID", B = 99, outdir = tempfile("dmsa_ex"))
#' @export
dmsa_frame <- function(data, methylation = NULL, map = "alpha",
                       reference = "alpha",
                       direction_source = c("auto", "cpgdirection", "bundled"),
                       tissue = "blood",
                       replicate_policy = c("collapse_site", "keep_correlated"),
                       cpgs_include = NULL, cpgs_exclude = NULL,
                       sample_id = NULL,
                       sample_align = c("id", "positional"),
                       missing_methylation = c("error", "drop_probes",
                                               "common_complete_rows"),
                       methylation_scale = c("auto", "beta", "M"),
                       system = TRUE, module = FALSE, gene = TRUE,
                       probe = TRUE, systems = NULL, sets = "alpha",
                       ...,
                       outcomes = NULL, predictors = NULL,
                       covariates = "contract", blocks = NULL,
                       random_effects = NULL,
                       chip = TRUE,
                       chip_effect = c("random", "fixed", "none"),
                       outcome_levels = NULL, predictor_levels = NULL,
                       moderation = FALSE, mod = "", mod2 = "",
                       outcome_label = NULL, predictor_labels = NULL,
                       mod_label = NULL, mod2_label = NULL,
                       covariate_labels = NULL, mod_levels = NULL,
                       type = c("linear", "non-linear", "exponential"),
                       outcome_type = c("gaussian", "logistic", "multinomial"),
                       frame_role = c("predictor", "outcome"),
                       cpg_map = c("full", "confidence"),
                       B = 1999, alpha = 0.05, fdr = 0.10, seed = 1,
                       correction = c("maxT", "minP"),
                       weighting = c("combined", "reliability", "flat"),
                       w_floor = 1.5,
                       palette = "viridis", plots = TRUE, tables = TRUE,
                       summary = TRUE, gene_models = "auto",
                       progress = interactive(), beep = interactive(),
                       plot_type = c("png", "pdf"),
                       table_type = c("html", "docx", "rtf"),
                       outdir = "dmsa_output", autofix = TRUE) {
  .fr_explicit <- !missing(frame_role)
  if (!is.numeric(fdr) || length(fdr) != 1L || !is.finite(fdr) ||
      fdr <= 0 || fdr >= 1)
    stop("`fdr` must be a single number between 0 and 1 (default .10)",
         call. = FALSE)
  chip_effect <- match.arg(chip_effect)
  type <- match.arg(type)
  frame_role <- match.arg(frame_role); cpg_map <- match.arg(cpg_map)

  ## ---- the phenotype columns: outcomes OR predictors (PI, 2026-08-29) -----
  ## One question, two spellings. `outcomes` = the columns the methylation is
  ## tested against; `predictors` = the same columns when the question runs
  ## the other way (methylation is the outcome), and naming them that way IS
  ## the orientation declaration - frame_role follows automatically. The
  ## former `outcome` argument is accepted through `...` so every existing
  ## script keeps working, but it is no longer part of the interface.
  .dots <- list(...)
  .legacy <- if ("outcome" %in% names(.dots)) .dots[["outcome"]] else NULL
  .dots[["outcome"]] <- NULL
  if (length(.dots))
    stop("unknown argument(s) to dmsa_frame(): ",
         paste(names(.dots), collapse = ", "),
         "\n(arguments after the first two must be given by their full ",
         "names)", call. = FALSE)
  .given <- c(outcomes = !is.null(outcomes),
              predictors = !is.null(predictors),
              outcome = !is.null(.legacy))
  if (sum(.given) > 1L)
    stop("give the phenotype columns ONCE: `",
         paste(names(.given)[.given], collapse = "` and `"),
         "` name the same columns (outcomes for outcome ~ frame, predictors ",
         "for frame ~ predictors). Pick one.", call. = FALSE)
  if (!is.null(predictors)) {
    if (.fr_explicit && frame_role == "predictor")
      stop("`predictors = ` says these columns PREDICT the methylation ",
           "frame, but frame_role = \"predictor\" says the frame predicts ",
           "them. Drop one of the two: `predictors` alone already sets the ",
           "orientation.", call. = FALSE)
    outcome <- predictors
    frame_role <- "outcome"
  } else outcome <- if (!is.null(outcomes)) outcomes else .legacy
  ## level labels follow the SPELLING the user declared the variable with
  ## (PI, 2026-08-29): someone who wrote `predictors = "pills_past_T1"`
  ## thinks of pills as a predictor and reaches for `predictor_levels`;
  ## `outcome_levels` is the same thing for `outcomes = `. One argument, in
  ## the user's own vocabulary - never both.
  if (!is.null(predictor_levels)) {
    if (!is.null(outcome_levels))
      stop("give the level labels ONCE: `predictor_levels` and ",
           "`outcome_levels` label the same variable. Use predictor_levels ",
           "with `predictors = `, outcome_levels with `outcomes = `.",
           call. = FALSE)
    if (is.null(predictors))
      stop("`predictor_levels` labels the variable given as `predictors = `,",
           " but the phenotype was given as `", 
           if (!is.null(outcomes)) "outcomes" else "outcome",
           " = `. Use `outcome_levels = ` to label it.", call. = FALSE)
    outcome_levels <- predictor_levels
  }
  ## display-label spellings follow the same principle. `predictor_labels`
  ## names the `predictors = ` column(s) - the PI's own call
  ## (predictor_labels = "Used birth-control pills when met") is the
  ## intended use. Covariate display names have their own argument,
  ## `covariate_labels`. The old behaviour (predictor_labels silently
  ## labelling the covariates) is gone: with `outcomes = ` it now errors
  ## and points to the right argument.
  if (!is.null(predictor_labels) && is.null(predictors))
    stop("`predictor_labels` names the `predictors = ` column(s), but the ",
         "phenotype was given as `outcomes = `. Use `outcome_label = ` for ",
         "it; covariate display names go in `covariate_labels = `.",
         call. = FALSE)
  if (!is.null(predictor_labels) && !is.null(outcome_label))
    stop("give the display name ONCE: `predictor_labels` and ",
         "`outcome_label` name the same column(s).", call. = FALSE)
  if (!is.null(predictor_labels)) outcome_label <- predictor_labels
  ## mod_label is ONE display name for the moderator column; a length-2
  ## vector is almost always intended as LEVEL labels - point there instead
  ## of silently keeping only the first entry
  if (!is.null(mod_label) && length(mod_label) > 1L)
    stop("`mod_label` is one display NAME for the moderator column. ",
         "Level labels (e.g. Husband/Wife for a 1/2 sex code) go in ",
         "`mod_levels = c(\"1\" = \"Husband\", \"2\" = \"Wife\")`.",
         call. = FALSE)
  if (is.null(outcome) || !length(outcome))
    stop("name the phenotype column(s) the methylation is tested against: ",
         "`outcomes =`, or - when the columns PREDICT the methylation ",
         "frame - `predictors =`.", call. = FALSE)
  correction <- match.arg(correction); plot_type <- match.arg(plot_type)
  weighting <- match.arg(weighting)
  ## spec 42 (PI ruling 2026-09-02): NO default exchangeability blocks. The
  ## user declares them with `blocks = ` when participants contribute more
  ## than one row; `random_effects = ` is kept as a deprecated alias of the
  ## same argument (it always named the blocks, not lme4 random effects).
  if (!is.null(random_effects)) {
    if (!is.null(blocks))
      stop("give the exchangeability blocks ONCE: `blocks = ` ",
           "(`random_effects = ` is the deprecated spelling of the same ",
           "argument)", call. = FALSE)
    message("dmsa_frame(): `random_effects = ` is deprecated - it names the ",
            "exchangeability BLOCKS (rows that travel together under ",
            "permutation); use `blocks = `. Treated as blocks = ",
            paste0("c(", paste0('"', random_effects, '"', collapse = ", "), ")"))
    blocks <- random_effects
  }
  random_effects <- blocks        # internal name below, unchanged
  .blocks_declared <- !is.null(blocks)
  table_type <- match.arg(table_type)
  sample_align <- match.arg(sample_align)
  missing_methylation <- match.arg(missing_methylation)
  methylation_scale <- match.arg(methylation_scale)
  direction_source <- match.arg(direction_source)
  replicate_policy <- match.arg(replicate_policy)

  ## ---- spec 21-24: analysis-LEVEL switches vs biological SELECTORS --------
  ## `system` is a TRUE/FALSE switch (run the system level?); `systems` picks
  ## which biological systems to analyse. Passing a system NAME to the singular
  ## argument used to bind the switch, leave `systems = NULL` (= every system),
  ## and run a whole-panel analysis the user never asked for. That must stop.
  .lvl_switch <- function(value, nm, plural, selects) {
    if (is.character(value) || is.factor(value))
      stop("`", nm, "` does NOT select ", selects, ".\n",
           "It is a TRUE/FALSE switch controlling whether the ",
           toupper(nm), " ANALYSIS LEVEL is run.\n\n",
           "You supplied:\n  ", nm, " = ",
           paste(sprintf('"%s"', as.character(value)), collapse = ", "),
           "\n\nTo select ", selects, " use the PLURAL argument:\n",
           "  ", plural, " = ",
           paste(sprintf('"%s"', as.character(value)), collapse = ", "), "\n\n",
           "No analysis was run.", call. = FALSE)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop("`", nm, "` must be a single TRUE or FALSE (it switches the ",
           toupper(nm), " analysis level on or off); got ",
           if (length(value) != 1L) paste0("length ", length(value), " ")
           else "", class(value)[1],
           if (length(value) == 1L && is.na(value)) " NA" else "",
           ". No analysis was run.", call. = FALSE)
    invisible(TRUE)
  }
  .lvl_switch(system, "system", "systems", "biological systems")
  .lvl_switch(module, "module", "modules", "biological modules")
  .lvl_switch(gene,   "gene",   "genes",   "genes")
  .lvl_switch(probe,  "probe",  "probes",  "probes")
  ## gene_models is NOT in the strict-logical list: it also accepts "auto"
  ## (fetch exon models for named genes, degrade offline) and a prebuilt
  ## dmsa_gene_model() table
  if (!(is.logical(gene_models) && length(gene_models) == 1L &&
        !is.na(gene_models)) &&
      !identical(gene_models, "auto") && !is.data.frame(gene_models))
    stop("`gene_models` must be TRUE, FALSE, \"auto\", or a table from ",
         "dmsa_gene_model()", call. = FALSE)
  for (.nm in c("moderation", "plots", "tables", "summary",
                "progress", "beep", "autofix")) {
    .v <- get(.nm)
    if (!is.logical(.v) || length(.v) != 1L || is.na(.v))
      stop("`", .nm, "` must be a single TRUE or FALSE; got ",
           if (length(.v) != 1L) paste0("length ", length(.v), " ") else "",
           class(.v)[1], ". No analysis was run.", call. = FALSE)
  }

  ## `systems = "all"` is the readable spelling of "every system" (NULL keeps
  ## working). It cannot be mixed with named systems - that request is
  ## contradictory and silently resolving it either way would be a guess.
  if (!is.null(systems) && is.character(systems) &&
      any(tolower(systems) == "all")) {
    if (length(systems) > 1L)
      stop("`systems = \"all\"` cannot be combined with named systems; ",
           "you supplied: ", paste(systems, collapse = ", "),
           ". Use systems = \"all\" alone, or name the systems you want.",
           call. = FALSE)
    systems <- NULL
  }
  ## dmsa_import() hands the frame both halves in one object: unwrap it so
  ## `dmsa_frame(dmsa_import(x, ...), outcome = ...)` is the whole bridge
  ## from any preprocessing pipeline to the analysis.
  if (inherits(data, "dmsa_import")) {
    if (is.null(methylation)) methylation <- data$methylation
    data <- data$data
  }
  data <- as.data.frame(data)
  cor <- data.frame(field = character(), issue = character(),
                    action = character(), stringsAsFactors = FALSE)

  ## ---- moderation layout validation ---------------------------------------
  if (isTRUE(moderation) && !nzchar(mod))
    stop("moderation = TRUE requires `mod` to name a column ",
         "(mod2 is optional and only together with mod)", call. = FALSE)
  if (nzchar(mod2) && !nzchar(mod))
    stop("`mod2` was set without `mod`; set `mod` first", call. = FALSE)
  if (!isTRUE(moderation) && (nzchar(mod) || nzchar(mod2))) {
    ## say it LIVE, not only in the corrections table: a user who typed
    ## mod = "sex" meant to test a moderation and would otherwise discover
    ## pages later that none was run (PI, 2026-08-29)
    message("dmsa_frame(): `mod = ` was given but moderation = FALSE - no ",
            "moderation is tested. Set moderation = TRUE to test it.")
    cor <- .frame_note(cor, "moderation",
                       "mod/mod2 given but moderation = FALSE",
                       "moderators ignored; set moderation = TRUE to test them")
  }
  for (mv in c(mod, mod2)) if (nzchar(mv) && !mv %in% names(data))
    stop("moderator column '", mv, "' is not in `data`", call. = FALSE)


  ## ---- display labels -----------------------------------------------------
  ## Column names are the user's own and often unreadable in a figure
  ## ("Attachment_Anxiety_General_T1"). Labels are cosmetic only: nothing in the
  ## model, the permutation or the correction ever sees them, so a mislabelled
  ## run is still a correct run - but every panel and sentence can read in
  ## plain language.
  .lab_pairs <- function(cols, labs, what) {
    if (is.null(labs) || !length(labs)) return(character(0))
    labs <- as.character(labs)
    if (!is.null(names(labs)) && all(nzchar(names(labs)))) {
      unknown <- setdiff(names(labs), cols)
      if (length(unknown))
        warning("`", what, "` names not in ", what, ": ",
                paste(unknown, collapse = ", "), " - ignored", call. = FALSE)
      labs <- labs[intersect(names(labs), cols)]
      return(labs)
    }
    if (length(labs) != length(cols))
      stop("`", what, "` has ", length(labs), " label(s) for ", length(cols),
           " column(s); supply one per column or a named vector",
           call. = FALSE)
    stats::setNames(labs, cols)
  }

  ## ---- outcomes -----------------------------------------------------------
  miss <- setdiff(outcome, names(data))
  if (length(miss))
    stop("outcome column(s) not in `data`: ", paste(miss, collapse = ", "),
         call. = FALSE)

  ## ---- spec 50: THE OUTCOME FAMILY IS PER OUTCOME, NOT PER FRAME ----------
  ## `outcome` is a vector, `outcome_type` was a single scalar, and the binary
  ## autofix below overwrote that scalar. So analysing
  ##   outcome = c("epds_total", "depressed")
  ## let the second, two-level outcome flip the family for the FIRST one, and
  ## the frame then reported a continuous score as logistic - with the switch
  ## recorded against the wrong column. A family belongs to an outcome.
  outcome_type <- .frame_otypes(if (missing(outcome_type)) "gaussian"
                                else outcome_type, outcome)

  ## ---- spec 51: THE MODEL ORIENTATION IS DECLARED AND CONSISTENT ---------
  ## Three models can be fitted from one frame, and they do not all face the
  ## same way:
  ##   main DMSA path   methylation ~ outcome + covariates   (always)
  ##   shape scan       outcome     ~ tone score + covariates
  ##   moderation       whichever `frame_role` declares
  ## The main path's direction is the DEFINITION of the DMSA estimand - a
  ## directional set statistic on methylation - so `frame_role` neither does
  ## nor can turn it around, and documenting it as "outcome ~ frame" was simply
  ## wrong. What `frame_role` governs is the secondary models. Until now only
  ## moderation read it: `frame_role = "outcome"` flipped the moderated model
  ## to `tone ~ outcome x moderator` while the shape scan accompanying it in
  ## the same report stayed `outcome ~ tone` - two answers to two different
  ## questions, printed side by side as though they agreed.
  .orient <- list(
    main = paste0("methylation ~ ", paste(outcome, collapse = " + "),
                  " + covariates"),
    shape = if (identical(type, "linear")) NA_character_
            else paste0(paste(outcome, collapse = " + "),
                        " ~ tone score + covariates"),
    moderation = if (!isTRUE(moderation)) NA_character_
                 else if (identical(frame_role, "outcome"))
                   paste0("tone score ~ ", paste(outcome, collapse = " + "),
                          " x moderator + covariates")
                 else paste0(paste(outcome, collapse = " + "),
                             " ~ tone score x moderator + covariates"))
  if (identical(frame_role, "outcome") && !identical(type, "linear"))
    stop("frame_role = \"outcome\" cannot be honoured together with type = \"",
         type, "\".\nThe shape scan that `type` switches on fits\n  ",
         .orient$shape, "\nwhich is the opposite orientation to the one ",
         "frame_role = \"outcome\" declares, and it is not switchable: the ",
         "curvature arms (Sasabuchi, Fieller) are defined on the outcome as ",
         "the response.\nSo the report would carry a moderated model facing ",
         "one way and a shape scan facing the other.\nEither declare ",
         "frame_role = \"predictor\", or keep frame_role = \"outcome\" and set ",
         "type = \"linear\".\nNo analysis was run.", call. = FALSE)
  if (identical(frame_role, "outcome") && !isTRUE(moderation))
    cor <- .frame_note(cor, "frame_role", "declared \"outcome\"",
                       paste0("no effect on this frame: only the moderated ",
                              "model reads it, and moderation = FALSE. The ",
                              "main path is ", .orient$main,
                              " either way - that is the DMSA estimand."))

  ## A nominal outcome cannot enter the main path. DMSA regresses methylation ON
  ## the outcome and tests ONE coefficient, so a k-level nominal variable coerced
  ## to numeric is silently analysed as an ordered score - "mode 5 is five times
  ## mode 1". That is not a defensible model and it must not run.
  if (any(outcome_type == "multinomial")) {
    .mn <- names(outcome_type)[outcome_type == "multinomial"]
    stop("outcome_type = 'multinomial' (outcome ",
         paste(.mn, collapse = ", "),
         ") is not supported by dmsa_frame(): the ",
         "main path regresses methylation ON the outcome and tests a single ",
         "1-df coefficient, so a nominal outcome would be coerced to a numeric ",
         "score and analysed as if its levels were equally spaced and ordered. ",
         "Use dmsa_outcome(), which implements a rank-based scheme for ",
         "multinomial and ordinal outcomes, or recode to a binary contrast ",
         "(e.g. one level vs the rest) and pass it as a two-level outcome.",
         call. = FALSE)
  }

  ## ---- spec 46: SAMPLE ALIGNMENT ------------------------------------------
  ## A methylation matrix and a phenotype table are two separate objects. If
  ## their rows are not the same samples in the same order, every coefficient
  ## in the report is computed from mismatched people - and nothing about the
  ## output looks wrong. Positional equality is therefore never ASSUMED when
  ## identifiers exist to check it against.
  if (is.matrix(methylation)) {
    if (nrow(methylation) != nrow(data))
      stop("methylation has ", nrow(methylation), " row(s) but `data` has ",
           nrow(data), ". They must describe the same samples.\n",
           "No analysis was run.", call. = FALSE)
    rn <- rownames(methylation)
    idcol <- sample_id
    if (sample_align == "id") {
      ## auto-detect only when the user did not name a column, and say so
      if (is.null(idcol) && !is.null(rn)) {
        ## numeric id columns (1001, 1002, ...) are as common as character
        ## ones; compare everything as character so they can be detected too
        cand <- names(data)[vapply(data, function(v)
          is.character(v) || is.factor(v) ||
            (is.numeric(v) && !anyNA(v)), TRUE)]
        hit <- cand[vapply(cand, function(v)
          setequal(as.character(data[[v]]), rn), TRUE)]
        if (length(hit) == 1L) {
          idcol <- hit
          message("dmsa_frame(): matching methylation rows by `", idcol,
                  "` (auto-detected from rownames). Pass sample_id = to be ",
                  "explicit.")
        }
      }
      if (!is.null(idcol)) {
        if (!idcol %in% names(data))
          stop("`sample_id` names a column that is not in `data`: ", idcol,
               call. = FALSE)
        if (is.null(rn))
          stop("`sample_id` was supplied but the methylation matrix has no ",
               "rownames to match against.\nEither add sample identifiers as ",
               "rownames, or set sample_align = \"positional\" if the rows are ",
               "known to be in the same order.\nNo analysis was run.",
               call. = FALSE)
        ids <- as.character(data[[idcol]])
        if (anyDuplicated(ids))
          stop("`", idcol, "` has duplicated identifier(s) (e.g. ",
               paste(utils::head(unique(ids[duplicated(ids)]), 3),
                     collapse = ", "), "); rows cannot be matched ",
               "unambiguously.\nNo analysis was run.", call. = FALSE)
        if (anyDuplicated(rn))
          stop("the methylation matrix has duplicated rownames (e.g. ",
               paste(utils::head(unique(rn[duplicated(rn)]), 3),
                     collapse = ", "), "); rows cannot be matched ",
               "unambiguously.\nNo analysis was run.", call. = FALSE)
        miss_d <- setdiff(ids, rn); miss_m <- setdiff(rn, ids)
        if (length(miss_d) || length(miss_m))
          stop("sample identifiers do not match between `data` and the ",
               "methylation matrix.\n",
               length(miss_d), " in data but not in methylation",
               if (length(miss_d)) paste0(" (e.g. ",
                 paste(utils::head(miss_d, 5), collapse = ", "), ")") else "",
               "\n", length(miss_m), " in methylation but not in data",
               if (length(miss_m)) paste0(" (e.g. ",
                 paste(utils::head(miss_m, 5), collapse = ", "), ")") else "",
               "\nNo analysis was run.", call. = FALSE)
        ord <- match(ids, rn)
        if (!identical(ord, seq_along(ord))) {
          methylation <- methylation[ord, , drop = FALSE]
          cor <- .frame_note(cor, idcol, "methylation rows in a different order",
                             "reordered to match `data` by sample identifier")
        }
      } else if (!is.null(rn)) {
        ## the METHYLATION side carries identifiers but no unique matching
        ## column was found in `data`. Recording "no identifiers on either
        ## side" here would be false, and positional alignment against a
        ## matrix that HAS rownames is exactly the unverifiable state spec 46
        ## exists to prevent - so this stops rather than proceeding.
        stop("the methylation matrix has rownames, but no column of `data` ",
             "matches them",
             if (length(hit) > 1L) paste0(" uniquely (", length(hit),
               " columns match: ", paste(utils::head(hit, 3), collapse = ", "),
               ")") else "",
             ".\nRows cannot be verified to describe the same samples in ",
             "the same order.\nName the identifier column with sample_id =, ",
             "or declare sample_align = \"positional\" if the order is known ",
             "to be right.\nNo analysis was run.", call. = FALSE)
      } else {
        ## no identifiers on either side: shuffling cannot be detected, so
        ## positional is the only possible behaviour. It is allowed, but it is
        ## recorded rather than assumed silently.
        cor <- .frame_note(cor, "sample alignment",
                           "no sample identifiers on either side",
                           paste("rows aligned POSITIONALLY and unverified;",
                                 "supply sample_id = and matrix rownames to",
                                 "have this checked"))
      }
    } else {
      cor <- .frame_note(cor, "sample alignment", "sample_align = positional",
                         "rows aligned by position at the user's request")
    }
  }

  ## ---- maps and coverage --------------------------------------------------
  ## Spec 1/3: WHERE direction comes from is now explicit. "cpgdirection"
  ## builds the map from CpG x gene PAIRS (cpgdirection::cpg_gene_pairs),
  ## the biological reference supplying membership; "bundled" is the interim
  ## per-column Alpha snapshot; "auto" (default) uses cpgdirection when it is
  ## installed and otherwise falls back to the snapshot WITH A RECORDED NOTE.
  ## A user-supplied `map` table is a bundled-style map and forces "bundled".
  cs <- .frame_sets(sets)
  .src <- direction_source
  if (is.data.frame(map) && .src != "bundled") {
    if (.src == "cpgdirection")
      stop("`map` was supplied as a table AND direction_source = ",
           "\"cpgdirection\" was requested. A user map is a bundled-style ",
           "per-column snapshot; pick one source.", call. = FALSE)
    .src <- "bundled"
  }
  if (.src == "auto")
    .src <- if (!is.null(.cpgd("cpg_gene_pairs"))) "cpgdirection" else "bundled"
  if (.src == "cpgdirection" && is.null(.cpgd("cpg_gene_pairs")) &&
      is.null(getOption("dmsa.pair_table", NULL)))
    stop("direction_source = \"cpgdirection\" but the cpgdirection package ",
         "(with cpg_gene_pairs) is not installed.\n",
         "BiocManager::install(\"teindor/cpgdirection\") installs the full ",
         "build; direction_source = \"bundled\" uses the interim snapshot.",
         call. = FALSE)
  if (direction_source == "auto" && .src == "bundled" &&
      !is.data.frame(map))
    cor <- .frame_note(cor, "direction source",
                       "cpgdirection is not installed",
                       paste("falling back to the bundled per-column Alpha",
                             "snapshot; install cpgdirection for pair-level",
                             "direction resolution"))

  pair_rt <- NULL
  if (.src == "cpgdirection") {
    ## the biological reference (spec 3): membership is system > module > gene
    .ref <- if (inherits(reference, "dmsa_reference")) reference
            else if (identical(reference, "alpha")) alpha_reference()
            else stop("`reference` must be \"alpha\" or a dmsa_reference ",
                      "object (see dmsa_reference(), dmsa_reference_csv())",
                      call. = FALSE)
    .avail0 <- if (is.matrix(methylation)) colnames(methylation)
               else if (is.character(methylation) && length(methylation) >= 1)
                 intersect(names(data), methylation)
               else names(data)
    if (!is.null(cpgs_include) || !is.null(cpgs_exclude)) {
      .sites0 <- .avail0[.map_is_site(.avail0)]
      .keep0 <- .cpg_name_filter(.sites0, cpgs_include, cpgs_exclude)
      if (!length(.keep0))
        stop("cpgs_include/cpgs_exclude removed every CpG-site column (",
             length(.sites0), " candidate(s)). Check the patterns - they are ",
             "matched as fixed substrings of the column names.",
             call. = FALSE)
      cor <- .frame_note(cor, "cpg selection",
                         sprintf("include: %s | exclude: %s",
                                 if (is.null(cpgs_include)) "-" else
                                   paste(cpgs_include, collapse = ", "),
                                 if (is.null(cpgs_exclude)) "-" else
                                   paste(cpgs_exclude, collapse = ", ")),
                         sprintf("%d of %d CpG column(s) kept",
                                 length(.keep0), length(.sites0)))
      .avail0 <- c(setdiff(.avail0, .sites0), .keep0)
    }
    pair_rt <- .frame_pairs_runtime(.avail0, .ref, tissue = tissue,
                                    replicate_policy = replicate_policy,
                                    pairs = getOption("dmsa.pair_table", NULL))
    maps <- pair_rt[c("full", "confidence")]
    ## rule 4 (integration brief): what probe QC excluded and what simply has
    ## no mapping is RECONCILED here, not left for the user to discover as a
    ## count mismatch. Both lists live on the ledger's attributes in full.
    ## the probe-QC note is deferred until AFTER the `systems =` subset: a
    ## file-wide exclusion count printed beside a two-system frame reads as
    ## the frame's own loss (PI, 2026-08-29). See the scoped note below.
    .orf <- attr(pair_rt$ledger, "cpgs_outside_reference")
    if (length(.orf))
      cor <- .frame_note(cor, "outside reference",
                         sprintf(paste0("%d CpG(s) map only to genes outside ",
                                        "the biological reference"),
                                 length(.orf)),
                         paste0("mapped fine, just not to the biology asked ",
                                "about (e.g. ",
                                paste(utils::head(.orf, 5), collapse = ", "),
                                "); not analysed"))
    .unm <- attr(pair_rt$ledger, "unmapped_cpgs")
    if (length(.unm))
      cor <- .frame_note(cor, "unmapped CpGs",
                         sprintf("%d CpG(s) map to no supported target",
                                 length(.unm)),
                         paste0("no CpG x gene pair exists for them (e.g. ",
                                paste(utils::head(.unm, 5), collapse = ", "),
                                "); they are not analysed and not errors"))
    cor <- .frame_note(cor, "direction source", "cpgdirection pairs",
                       sprintf(paste0("%d usable pair(s) from %d CpG(s) x %d ",
                                      "gene(s); tissue %s; ledger kept on ",
                                      "the frame"),
                               nrow(pair_rt$full),
                               length(unique(pair_rt$full$probe)),
                               length(unique(pair_rt$full$gene)), tissue))
  } else {
    maps <- .frame_maps(map)
  }

  ## translate the map's column names to whichever build this data frame is,
  ## BEFORE the intersection that decides whether any probe was found at all
  ## spec 47: a character `methylation` selects columns however MANY it names.
  ## The old guard was length > 1, so methylation = "cg12345678" fell through
  ## to "use every column of `data`" - a single-probe request silently became
  ## a whole-panel analysis. Any non-NULL character selection is honoured.
  .avail <- if (is.matrix(methylation)) colnames(methylation)
            else if (is.character(methylation) && length(methylation) >= 1)
              intersect(names(data), methylation)
            else names(data)
  if (is.null(pair_rt)) {
    ## bundled snapshot: build-specific column-name translation applies
    .res <- .frame_resolve_columns(maps[[cpg_map]], .avail, cs,
              if (is.character(methylation) && length(methylation) >= 1) NULL
              else .frame_which_build(names(data)))
    if (!is.null(.res$note)) {
      for (.k in names(maps))
        maps[[.k]]$column <- .frame_resolve_columns(maps[[.k]], .avail, cs,
          if (is.character(methylation) && length(methylation) >= 1) NULL
          else .frame_which_build(names(data)))$column
      cor <- .frame_note(cor, "column names", "build-specific CpG naming",
                         .res$note)
    }
  }
  ## pair path: the map was BUILT from the available columns, so no
  ## translation and no membership filter is needed - the columns are the
  ## user's own measurement names (or, under collapse_site, the canonical
  ## CpG ids whose replicate columns are averaged below).

  chosen_map <- maps[[cpg_map]]
  if (is.null(pair_rt)) {
    in_data <- chosen_map$column %in% names(data)
    if (is.character(methylation) && length(methylation) >= 1)
      in_data <- in_data & chosen_map$column %in% methylation
    if (is.matrix(methylation)) {
      keyed <- chosen_map$column %in% colnames(methylation) |
               chosen_map$probe  %in% colnames(methylation)
      in_data <- keyed
    }
    mp <- chosen_map[in_data, , drop = FALSE]
  } else mp <- chosen_map
  if (is.null(pair_rt) &&
      (!is.null(cpgs_include) || !is.null(cpgs_exclude))) {
    ## bundled path: the same name filter, applied to the mapped columns
    ## AFTER build translation, so the patterns match the names actually
    ## being analysed
    .keepb <- .cpg_name_filter(mp$column, cpgs_include, cpgs_exclude)
    .nb <- nrow(mp)
    mp <- mp[mp$column %in% .keepb, , drop = FALSE]
    cor <- .frame_note(cor, "cpg selection",
                       sprintf("include: %s | exclude: %s",
                               if (is.null(cpgs_include)) "-" else
                                 paste(cpgs_include, collapse = ", "),
                               if (is.null(cpgs_exclude)) "-" else
                                 paste(cpgs_exclude, collapse = ", ")),
                       sprintf("%d of %d mapped column(s) kept", nrow(mp), .nb))
  }
  if (!nrow(mp))
    stop("no mapped methylation columns found (cpg_map = '", cpg_map,
         "'); ",
         if (!is.null(cpgs_include) || !is.null(cpgs_exclude))
           "the cpgs_include/cpgs_exclude patterns removed every mapped column - check them first. Otherwise "
         else "",
         if (is.null(pair_rt))
           "check `methylation`/`data` column names against the map"
         else paste0("no usable pair rests on tier M or S1 evidence here. ",
                     "cpg_map = 'full' analyses every usable pair; the pair ",
                     "ledger says what each pair's direction rests on."),
         call. = FALSE)
  ## (untestable-gene tracking happens AFTER the systems subset below - a
  ## 1.1.0 bug counted untestable genes across ALL systems in the data)

  sys_tab <- unique(mp[, c("system_id", "system")])
  if (!is.null(systems)) {
    sys_tab <- .frame_pick_systems(systems, sys_tab, cs)
    gone <- attr(sys_tab, "absent")
    if (length(gone))
      cor <- .frame_note(cor, "systems", "selected system carries no probe here",
                         paste("no mapped probe in the data for:",
                               paste(gone, collapse = ", ")))
    mp <- mp[mp$system_id %in% sys_tab$system_id, , drop = FALSE]
  }
  mp <- mp[order(mp$system_id, mp$gene, mp$probe), , drop = FALSE]

  ## probe-QC accounting, SCOPED TO THE FRAME (PI ruling, 2026-08-29): the
  ## headline number is how many excluded CpGs would have mapped into the
  ## systems actually analysed; the whole-file count rides along as context.
  ## The mapping of the excluded CpGs comes from one extra flag-mode call
  ## over just those CpGs (see .frame_cpg_gene_pairs); when it is not
  ## available (injected pair tables, or the call failed) the note falls
  ## back to the file-wide count WITH its scope stated, never implied.
  if (!is.null(pair_rt)) {
    .qcx <- attr(pair_rt$ledger, "qc_excluded_cpgs")
    if (length(.qcx)) {
      .qcp <- attr(pair_rt$ledger, "qc_excluded_pairs")
      .act <- paste0("probe QC mask: cross-hybridizing/degenerate mapping/",
                     "SNP-contaminated; probe_qc = \"flag\" in cpgdirection ",
                     "recovers them marked")
      if (!is.null(.qcp) && nrow(.qcp)) {
        ## the scope: with an explicit `systems =`, the selected systems'
        ## reference genes (even genes with no surviving probe - the excluded
        ## CpG could have been that probe); with no selection, the whole
        ## reference, because everything in it is being analysed
        .selg <- if (is.null(systems)) unique(.ref$systems$gene) else
          unique(.ref$systems$gene[
            .ref$systems$system_id %in% sys_tab$system_id |
            .ref$systems$system %in% sys_tab$system])
        .inx <- unique(.qcp$cpg_id[.qcp$target_gene %in% .selg])
        .iss <- sprintf(paste0("%d CpG(s) mapping into the analysed ",
                               "system(s) excluded as unreliable probes"),
                        length(.inx))
        if (length(.inx))
          .act <- paste0(.act, " (e.g. ",
                         paste(utils::head(.inx, 5), collapse = ", "), ")")
        .act <- sprintf("%s; %d excluded across the whole submitted file",
                        .act, length(.qcx))
        message("dmsa_frame(): ", .iss, " (", length(.qcx),
                " across the whole file). probe_qc = \"flag\" keeps them, ",
                "marked.")
      } else {
        .iss <- sprintf(paste0("%d CpG(s) excluded as unreliable probes ",
                               "across the whole submitted file (their ",
                               "mapping is unknown, so this count is NOT ",
                               "scoped to the analysed systems)"),
                        length(.qcx))
        .act <- paste0(.act, " (e.g. ",
                       paste(utils::head(.qcx, 5), collapse = ", "), ")")
      }
      cor <- .frame_note(cor, "probe QC", .iss, .act)
    }
  }

  ## genes present in the raw map for the CHOSEN systems but with ZERO probes
  ## surviving the chosen direction bar are UNTESTABLE - tracked, reported,
  ## never silently dropped
  ## The raw map is re-read from disk, so it still carries the parent-T1
  ## column names - the build translation applied to `maps` above is NOT in it.
  ## Testing raw$column against names(data) therefore matched NOTHING in the
  ## child and before-after builds, `ut` came out empty, and every one of those
  ## reports claimed "Coverage: 0 gene(s) ... untestable" while the parent build
  ## on the SAME two systems and the SAME 18 genes correctly reported 36.
  ## Translate raw the same way before the intersection.
  if (is.null(pair_rt)) {
    raw <- .frame_read_map(map)
    raw$column <- .frame_resolve_columns(
      raw, .avail, cs,
      if (is.character(methylation) && length(methylation) >= 1) NULL
      else .frame_which_build(names(data)))$column
    raw_in <- if (is.matrix(methylation))
                raw$column %in% colnames(methylation) |
                raw$probe  %in% colnames(methylation)
              else if (is.character(methylation) && length(methylation) >= 1)
                raw$column %in% names(data) & raw$column %in% methylation
              else raw$column %in% names(data)
    raw <- raw[raw$system_id %in% sys_tab$system_id & raw_in, , drop = FALSE]
    ut <- setdiff(unique(raw$gene), unique(mp$gene))
    untestable <- unique(raw[raw$gene %in% ut, c("system_id", "system", "gene")])
  } else {
    ## spec 15: on the pair path the REFERENCE defines membership, so an
    ## untestable gene is a reference gene of a chosen system with no usable
    ## pair here - whether unmeasured, abstained, or QC-excluded. The ledger
    ## carries the pair-by-pair reason; this table is the per-gene summary.
    .rs <- if (inherits(.ref, "dmsa_reference")) .ref$systems else .ref
    .rs <- .rs[as.character(.rs$system_id) %in%
                 as.character(sys_tab$system_id), , drop = FALSE]
    ut <- setdiff(unique(.rs$gene), unique(mp$gene))
    untestable <- unique(.rs[.rs$gene %in% ut,
                             c("system_id", "system", "gene"), drop = FALSE])
  }

  ## modules requested -> module map must cover the chosen systems
  mods <- NULL
  if (isTRUE(module)) {
    mm <- .frame_modules(sets)
    if (!is.null(mm) && !any(mm$gene %in% mp$gene)) mm <- NULL
    have <- if (is.null(mm)) integer(0) else
      intersect(unique(mm$system_id), sys_tab$system_id)
    lack <- setdiff(sys_tab$system_id, have)
    if (!length(have))
      stop("module = TRUE but no module map covers the chosen system(s) ",
           paste(lack, collapse = ", "),
           "; the bundled module map covers systems ",
           paste(sort(unique(mm$system_id)), collapse = ", "), call. = FALSE)
    if (length(lack)) {
      msgt <- paste("module level runs only for system(s)",
                    paste(have, collapse = ", "), "- no module map for",
                    paste(lack, collapse = ", "))
      if (!autofix) stop(msgt, call. = FALSE)
      cor <- .frame_note(cor, "module", "module map incomplete", msgt)
    }
    mods <- mm[mm$system_id %in% have, , drop = FALSE]
  }

  ## ---- covariates + autofixable derived columns ---------------------------
  if (identical(covariates, "contract")) {
    .bld <- .frame_which_build(names(data))
    covariates <- .frame_contract(names(data))
    if (!is.null(.bld)) {
      ct <- .frame_contracts[[.bld]]
      cor <- .frame_note(cor, "contract", ct$label,
                         paste0("fixed: ", paste(ct$fixed, collapse = ", "),
                                " | blocks: ", paste(ct$block, collapse = ", "),
                                if (is.null(ct$chip)) " | chip: NEVER included"
                                else paste0(" | chip: ", ct$chip)))
      ## the contract the user asked for DECLARES the blocks; adopt them when
      ## the user gave none (an explicit `blocks = ` always wins)
      if (is.null(random_effects) && all(ct$block %in% names(data))) {
        random_effects <- ct$block
        .blocks_declared <- TRUE
        cor <- .frame_note(cor, "blocks",
                           paste0("set by the ", ct$label, " contract"),
                           paste0("blocks = ", paste(ct$block, collapse = ", "),
                                  " (covariates = \"contract\" carries the block structure)"))
        message("dmsa_frame(): blocks = ", paste(ct$block, collapse = ", "),
                " taken from the ", ct$label, " contract")
      }
      ## the before/after builds forbid chip; honour that
      if (is.null(ct$chip)) attr(covariates, "no_chip") <- TRUE
    }
  }
  ## sex_c is DERIVED, not a column of any workbook. Deriving it only under
  ## covariates = "contract" meant a caller who typed the contract out by hand
  ## silently lost the sex adjustment - same intent, different model, different
  ## p-values. Derive it whenever it is asked for.
  if ("sex_c" %in% covariates && !"sex_c" %in% names(data) &&
      "sex" %in% names(data)) {
    sx <- suppressWarnings(as.numeric(data$sex))
    data$sex_c <- ifelse(sx == 1, -0.5, +0.5)
    cor <- .frame_note(cor, "sex_c", "absent",
                       "built from `sex` (1 -> -0.5, else +0.5)")
  }
  ## `chip` is an argument, not only an auto-detection: on the Alpha parent
  ## build chip_T1 has 61 levels on ~400 arrays (7 of them singletons), which
  ## spends 60 df on a nuisance factor that ctrlSV3/ctrlSV5 already partly
  ## absorb. Whether to pay that is a judgement the analyst has to be able to
  ## make and report, so chip = FALSE has to be sayable.
  ## ---- random_effects: mixed-model-style routing (PI, 2026-08-29) ---------
  ## Users state their random effects the way they think about them - as
  ## grouping factors, lme4-style. With TWO factors, dmsa_frame() works out
  ## the structure itself and routes each factor to the machinery that
  ## implements it in a permutation engine:
  ##   the FINEST dependence factor  -> the exchangeability blocks (rows that
  ##                                    must travel together under permutation)
  ##   the other (batch-like) factor -> the (1 | .) random intercept, i.e. the
  ##                                    chip machinery (REML quasi-demeaning)
  ## The decision is announced in lme4 notation. Crossing the factors into one
  ## block interaction - the old reading - made every group a singleton
  ## whenever partners sat on different chips, and stopped the run.
  random_effects <- unique(as.character(
    random_effects[!is.na(random_effects) & nzchar(random_effects)]))
  if (length(random_effects) > 2L)
    stop("blocks lists ", length(random_effects), " grouping ",
         "factors. dmsa supports two: one exchangeability block (the rows ",
         "that stay together under permutation) and one (1 | .) random ",
         "intercept. Fold or drop the rest.", call. = FALSE)
  if (length(random_effects) == 2L && all(random_effects %in% names(data))) {
    .f1n <- random_effects[1]; .f2n <- random_effects[2]
    .f1 <- factor(data[[.f1n]]); .f2 <- factor(data[[.f2n]])
    .n12 <- all(rowSums(table(.f1, .f2) > 0) == 1L)   # f1 nested in f2
    .n21 <- all(rowSums(table(.f2, .f1) > 0) == 1L)   # f2 nested in f1
    if (.n12 && .n21) {
      ## the two columns induce the SAME grouping: keeping both is harmless
      ## (their interaction IS that grouping), so only say so
      cor <- .frame_note(cor, "blocks",
                         sprintf("`%s` and `%s` are aliased (identical grouping)",
                                 .f1n, .f2n),
                         "one grouping, stated twice - treated as one block")
    } else if (.n12 || .n21) {
      ## NESTED, e.g. (1 | chip/cID) with couples on the SAME chip, or the
      ## before/after contract's (ID within cID). The interaction of nested
      ## factors IS the fine factor, so the existing blocks already carry the
      ## dependence exactly - nothing needs rerouting, and nothing changes.
      ## The coarse factor's own (1 | .) intercept, when it is a chip, comes
      ## through `chip =` (whose default auto-detects chip columns).
      .fine <- if (.n12) .f1n else .f2n
      .coarse <- setdiff(c(.f1n, .f2n), .fine)
      cor <- .frame_note(cor, "blocks",
                         sprintf("(1 | %s/%s): %s nested in %s",
                                 .coarse, .fine, .fine, .coarse),
                         sprintf(paste0("blocks on the nested pair equal ",
                           "blocks on `%s` alone - rows of one %s travel ",
                           "together under permutation; a (1 | %s) variance ",
                           "component is the `chip =` argument's job"),
                           .fine, .fine, .coarse))
    } else {
      ## CROSSED/INDEPENDENT, e.g. partners of one couple on DIFFERENT
      ## chips: their interaction is all singletons and can never permute,
      ## so this case used to stop the run. Route instead: the finest
      ## dependence factor is the exchangeability block, the other enters as
      ## the (1 | .) random intercept through the chip machinery.
      .fine <- if (nlevels(.f1) >= nlevels(.f2)) .f1n else .f2n
      .coarse <- setdiff(c(.f1n, .f2n), .fine)
      if (isTRUE(attr(covariates, "no_chip")))
        stop("blocks = c(\"", .f1n, "\", \"", .f2n, "\") are ",
             "CROSSED, which needs `", .coarse, "` as a (1 | .) random ",
             "intercept - but this build's contract forbids a chip/batch ",
             "term. Use one exchangeability factor, or drop the contract.",
             call. = FALSE)
      if (is.character(chip) && length(chip) == 1L && nzchar(chip) &&
          !identical(chip, .coarse) && !identical(chip, "pool"))
        stop("blocks = c(\"", .f1n, "\", \"", .f2n, "\") routes `",
             .coarse, "` to the (1 | .) random intercept, but chip = \"",
             chip, "\" already claims that slot. The engine fits ONE random ",
             "intercept: drop one of the two.", call. = FALSE)
      if (identical(chip_effect, "none"))
        stop("blocks lists `", .coarse, "` as a random effect, but ",
             "chip_effect = \"none\" says to drop it. Pick one.",
             call. = FALSE)
      cor <- .frame_note(cor, "blocks",
                         sprintf("(1 | %s) + (1 | %s): crossed/independent",
                                 .f1n, .f2n),
                         sprintf(paste0("`%s` (finest dependence) is the ",
                           "exchangeability block - its rows travel together ",
                           "under permutation; `%s` enters as the (1 | %s) ",
                           "random intercept (REML quasi-demeaning). Crossing ",
                           "them into one block would make every group a ",
                           "singleton."),
                           .fine, .coarse, .coarse))
      chip <- .coarse
      if (!identical(chip_effect, "fixed")) chip_effect <- "random"
      random_effects <- .fine
    }
  }

  chip_req <- chip
  chip <- NULL
  if (isTRUE(attr(covariates, "no_chip"))) {
    cor <- .frame_note(cor, "chip", "forbidden by this build's contract",
                       "not entered: the before/after builds list chip as NEVER include")
  } else if (identical(chip_req, FALSE)) {
    cor <- .frame_note(cor, "chip", "chip = FALSE",
                       paste0("not entered as a fixed effect - residual chip ",
                              "effects are left to ctrlSV3/ctrlSV5; say so in ",
                              "the write-up"))
  } else if (is.character(chip_req) && length(chip_req) == 1L && nzchar(chip_req)) {
    pool <- identical(chip_req, "pool") && !("pool" %in% names(data))
    src <- if (pool) {
      if ("chip_f" %in% names(data)) "chip_f" else
        grep("^chip(_T[0-9]+)?$", names(data), value = TRUE)[1]
    } else chip_req
    if (is.na(src) || !src %in% names(data))
      stop("chip = \"", chip_req, "\" is not a column of `data`", call. = FALSE)
    cf <- factor(data[[src]])
    if (pool) {
      ## A chip level holding a single array is fitted exactly by its own dummy:
      ## leverage 1, residual 0. That array is counted in n and contributes
      ## nothing to any estimate. On the Alpha parent build 7 of 61 chips are
      ## singletons, so 7 of 396 arrays are structurally uninformative while
      ## still being reported. Pooling the rare levels keeps the batch
      ## adjustment, recovers the degrees of freedom, and removes the
      ## leverage-1 rows.
      tb <- table(cf)
      rare <- names(tb)[tb < 3L]
      if (length(rare)) {
        lv <- as.character(cf); lv[lv %in% rare] <- "chip_small"
        cf <- factor(lv)
        cor <- .frame_note(cor, "chip_f", "chip = \"pool\"",
          sprintf(paste0("%d chip level(s) holding fewer than 3 arrays pooled ",
                         "into 'chip_small' (%d array(s)); %d level(s) remain, ",
                         "%d df instead of %d"),
                  length(rare), sum(tb[rare]), nlevels(cf), nlevels(cf) - 1L,
                  length(tb) - 1L))
      } else {
        cor <- .frame_note(cor, "chip_f", "chip = \"pool\"",
                           "no chip level held fewer than 3 arrays - nothing pooled")
      }
    } else {
      cor <- .frame_note(cor, "chip_f", paste0("chip = \"", chip_req, "\""),
                         paste0("built as factor(", src, "), fixed effect"))
    }
    data$chip_f <- cf; chip <- "chip_f"
  } else
  if ("chip_f" %in% names(data)) chip <- "chip_f"
  else {
    chip_col <- grep("^chip(_T[0-9]+)?$", names(data), value = TRUE)[1]
    if (!is.na(chip_col)) {
      data$chip_f <- factor(data[[chip_col]]); chip <- "chip_f"
      cor <- .frame_note(cor, "chip_f", "absent",
                         paste0("built as factor(", chip_col, ")"))
    }
  }

  ## The Alpha contract enters chip as a RANDOM intercept, `(1 | chip_T1)`, not
  ## as a fixed factor. Through 1.17.0 dmsa used a fixed factor: 61 levels and
  ## 60 df on the parent build, with 7 singleton chips fitted perfectly. Under
  ## `chip_effect = "random"` (the default, and what the contract specifies)
  ## chip is removed from the fixed design and handled by a quasi-demeaning GLS
  ## transform instead - a median of ~16 effective df on these data, and zero on
  ## the probes where the chip variance is estimated as nil.
  chip_random <- ""
  if (identical(chip_effect, "none")) {
    if (.nz(chip)) cor <- .frame_note(cor, "chip", "chip_effect = 'none'",
      "chip was NOT entered in any form")
    chip <- ""
  } else if (identical(chip_effect, "random") && .nz(chip)) {
    chip_random <- chip; chip <- ""
    cor <- .frame_note(cor, "chip", "entered as a RANDOM intercept",
      paste0("(1 | chip) per the Alpha covariate contract - fitted by REML ",
             "quasi-demeaning, not as a ", nlevels(factor(data$chip_f)),
             "-level fixed factor"))
  } else if (identical(chip_effect, "fixed") && .nz(chip)) {
    cor <- .frame_note(cor, "chip", "entered as a FIXED factor",
      paste0("deviation from the Alpha contract, which specifies (1 | chip); ",
             "costs ", nlevels(factor(data$chip_f)) - 1L, " df"))
  }
  miss <- setdiff(covariates, names(data))
  if (length(miss)) {
    if (!autofix)
      stop("covariate(s) not in `data`: ", paste(miss, collapse = ", "),
           call. = FALSE)
    cor <- .frame_note(cor, paste(miss, collapse = ","),
                       "covariate absent from data", "dropped from the model")
    ## a mistyped covariate name is indistinguishable from an intended one;
    ## dropping it silently produces an UNADJUSTED model that still looks
    ## adjusted in the write-up. Make it a warning, not a table footnote.
    warning("covariate(s) not found in `data` and dropped from the model: ",
            paste(miss, collapse = ", "),
            " - the reported model is NOT adjusted for them", call. = FALSE)
    covariates <- setdiff(covariates, miss)
  }

  ## A contract covariate derived from an outcome must never enter the model.
  ## Without this the frame builds sex_c from `sex`, then discovers the design
  ## is rank-deficient, then drops it - three messages that leave the user to
  ## work out for themselves that the covariate WAS the outcome.
  if (length(covariates)) {
    drop_oc <- character(0)
    for (oc0 in outcome) {
      hits <- grep(paste0("^", oc0, "(_c|_z|_f)?$"), covariates, value = TRUE)
      drop_oc <- c(drop_oc, hits)
    }
    drop_oc <- unique(setdiff(drop_oc, character(0)))
    if (length(drop_oc)) {
      covariates <- setdiff(covariates, drop_oc)
      cor <- .frame_note(cor, paste(drop_oc, collapse = ","),
                         "derived from an outcome column",
                         paste0("excluded before fitting: a covariate built ",
                                "from ", paste(outcome, collapse = "/"),
                                " cannot adjust a model that tests it"))
    }
    ## a column given BOTH as a covariate and as a moderator (PI, 2026-08-29:
    ## "sex" in covariates AND mod = "sex") would enter the moderated model
    ## twice - once raw, once as the scaled moderator - a perfect
    ## collinearity the fit would resolve silently. Say it out loud, drop
    ## the covariate copy, and tell the user how to undo it.
    .dup_mod <- intersect(covariates, Filter(nzchar, c(mod, mod2)))
    if (length(.dup_mod) && isTRUE(moderation)) {
      covariates <- setdiff(covariates, .dup_mod)
      message("dmsa_frame(): ", paste0("'", .dup_mod, "'", collapse = ", "),
              " was given both as a covariate and as a moderator. It is ",
              "dropped from the covariates - the moderated model already ",
              "carries its main effect. If that was not intended, re-run ",
              "dmsa_frame() with the intended setting.")
      cor <- .frame_note(cor, paste(.dup_mod, collapse = ","),
                         "given as covariate AND moderator",
                         paste0("dropped from the covariates; the moderated ",
                                "model carries its main effect once"))
    } else if (length(.dup_mod)) {
      ## mod named but moderation = FALSE: the moderator is IGNORED (a note
      ## already records that), so the covariate copy is the one that runs -
      ## still worth a live line, because the user plainly meant to moderate
      message("dmsa_frame(): '", paste(.dup_mod, collapse = "', '"),
              "' is a covariate AND `mod = `, but moderation = FALSE, so ",
              "no moderation is tested - it stays a plain covariate. Set ",
              "moderation = TRUE to test the interaction.")
    }
  }

  ## ---- TEST DRIVE 1: numeric coercion with per-column loss ---------------
  ## spec 45: CATEGORICAL COVARIATES MUST STAY CATEGORICAL. Blanket
  ## as.numeric(as.character(x)) had two failure modes, both silent:
  ##   - levels like "control"/"treated" became all-NA, so every row dropped;
  ##   - levels like "1"/"2"/"3" became 1,2,3 - a categorical adjusted as a
  ##     CONTINUOUS variable with equal spacing, which is a wrong adjustment
  ##     that no error ever announced.
  ## A column is only coerced when every non-missing value genuinely parses as
  ## a number (numbers stored as text). Otherwise it is kept as a factor and
  ## model.matrix() builds the contrasts, which is what the fitting path
  ## already expects. The outcome keeps the old treatment: its family is
  ## declared through `outcome_type`, not inferred here.
  catvars <- unique(c(covariates, c(mod, mod2)[nzchar(c(mod, mod2))]))
  numvars <- unique(c(outcome, covariates,
                      c(mod, mod2)[nzchar(c(mod, mod2))]))
  for (v in numvars) {
    old <- data[[v]]
    if (is.numeric(old)) next
    ## A logical column is already a two-level numeric in disguise; sending it
    ## through as.numeric(as.character(x)) turns "TRUE" into NA and deletes
    ## every row. TRUE -> 1, FALSE -> 0 is the standard, lossless reading.
    if (is.logical(old)) {
      data[[v]] <- as.numeric(old)
      cor <- .frame_note(cor, v, "logical column",
                         "converted TRUE/FALSE to 1/0")
      next
    }
    new <- suppressWarnings(as.numeric(as.character(old)))
    lost <- sum(is.na(new) & !is.na(old))
    ## An explicit factor/ordered column is a DECLARATION of categorical
    ## intent: honour it even when the labels happen to look numeric
    ## (levels "1","2","3" are site codes, not a continuous variable).
    ## Only bare character columns are parse-tested.
    if (v %in% catvars && (is.factor(old) || lost > 0)) {
      f <- droplevels(as.factor(old))
      nl <- nlevels(f)
      if (nl < 2L)
        stop("covariate '", v, "' is categorical with ", nl,
             " usable level(s); it cannot adjust anything. ",
             "Drop it from `covariates` or supply a variable that varies.",
             call. = FALSE)
      data[[v]] <- f
      cor <- .frame_note(cor, v, "categorical",
                         sprintf(paste0("kept as a factor with %d levels (%s); ",
                                        "model.matrix() builds the contrasts"),
                                 nl, paste(utils::head(levels(f), 4),
                                           collapse = ", ")))
      next
    }
    if (!autofix && lost > 0)
      stop("column '", v, "' is non-numeric and coercion would lose ",
           lost, " values", call. = FALSE)
    data[[v]] <- new
    cor <- .frame_note(cor, v, "non-numeric",
                       sprintf("coerced to numeric%s",
                               if (lost > 0) sprintf(" (%d values -> NA)", lost)
                               else ""))
  }

  ## complete cases over model variables + blocks + methylation columns
  ## spec 43: a REQUESTED exchangeability block may never silently disappear.
  ## Dropping it does not merely lose a covariate - it changes the NULL MODEL
  ## from within-block to unrestricted permutation, so every p-value afterwards
  ## answers a question the user did not ask. This is a hard error even under
  ## autofix, which is allowed to repair nuisances, never the null.
  blockv <- intersect(random_effects, names(data))
  if (!length(random_effects)) {
    ## spec 42: say it, live, and ask whether it is intended
    message("dmsa_frame(): no exchangeability blocks declared (blocks = NULL) - ",
            "the permutation null treats every row as independent. Make sure ",
            "this is intended: if any participant contributes more than one ",
            "row (couples, family members, repeated measures), pass ",
            "blocks = \"<id column>\"; if every row is a different person, ",
            "nothing to do.")
    cor <- .frame_note(cor, "blocks", "none declared",
                       "rows permuted independently - intended? pass blocks = <id column> if participants repeat")
  }
  if (length(blockv) < length(random_effects)) {
    gone <- setdiff(random_effects, blockv)
    stop("requested exchangeability block column(s) not found in `data`: ",
         paste(gone, collapse = ", "),
         "\nContinuing without them would change the null model to ",
         "unrestricted permutation.\nNo analysis was run.", call. = FALSE)
  }
  ## chip_random must be RETAINED even though it is not in the fixed design:
  ## it is the grouping vector for the random intercept. Dropping it here made
  ## ri_group NULL and the random intercept silently never applied.
  needcols <- unique(c(numvars, blockv, if (.nz(chip)) chip,
                       if (.nz(chip_random)) chip_random))
  cc <- stats::complete.cases(data[, needcols, drop = FALSE])
  loss <- vapply(needcols, function(v) sum(is.na(data[[v]])), integer(1))
  worst <- needcols[which.max(loss)]
  if (any(loss > 0))
    cor <- .frame_note(cor, "complete cases",
                       sprintf("%d of %d rows dropped", sum(!cc), length(cc)),
                       sprintf("costliest column: %s (%d missing)",
                               worst, max(loss)))
  base <- data[cc, , drop = FALSE]
  if (nrow(base) < 30)
    stop("only ", nrow(base), " complete rows remain - too few", call. = FALSE)
  ## spec 43 (second half): blocks that survive as singletons are not blocks.
  ## An interaction such as cID x chip can leave every group of size 1, in which
  ## case within-block permutation has nothing to permute and silently degenerates
  ## to the identity - or, worse, to unrestricted permutation downstream.
  if (length(blockv)) {
    .bk <- interaction(base[, blockv, drop = FALSE], drop = TRUE)
    .sz <- table(.bk)
    if (all(.sz < 2L))
    {
      ## Say not only WHAT went wrong but WHICH column caused it: with several
      ## block columns, crossing them is usually the mistake (e.g. couples on
      ## DIFFERENT chips make cID x chip all-singleton), and one of the
      ## columns alone typically permutes fine - while a chip/batch column
      ## belongs in `chip =`, where it enters as a random intercept instead.
      .hint <- ""
      if (length(blockv) > 1L) {
        .solo <- vapply(blockv, function(v)
          sum(table(base[[v]]) >= 2L), integer(1))
        .good <- names(.solo)[.solo >= 2L]
        .batchy <- blockv[grepl("chip|slide|plate|batch|array", blockv,
                                ignore.case = TRUE)]
        .hint <- paste0(
          if (length(.good))
            paste0("\nNote: ", paste(sprintf("`%s` alone gives %d permutable group(s)",
                   .good, .solo[.good]), collapse = "; "),
                   " - crossing the columns is what created the singletons.")
          else "",
          if (length(.batchy))
            paste0("\n", paste0("`", .batchy, "`", collapse = " and "),
                   " looks like a batch/chip variable: it belongs in ",
                   "`chip = ` (random intercept), not in the ",
                   "exchangeability block.")
          else "")
      }
      stop("The requested exchangeability block has no permutable groups.\n",
           "Block: ", paste(blockv, collapse = " x "), " -> ", length(.sz),
           " group(s), all of size 1 on the ", nrow(base),
           " analysable rows.\n",
           "Continuing would change the null model to unrestricted permutation.",
           .hint,
           "\nNo analysis was run.", call. = FALSE)
    }
    if (sum(.sz >= 2L) < 2L)
      cor <- .frame_note(cor, paste(blockv, collapse = " x "),
                         sprintf("only %d block(s) of size >= 2", sum(.sz >= 2L)),
                         "permutation is nearly degenerate - interpret with care")
  }
  if (.nz(chip)) base[[chip]] <- droplevels(factor(base[[chip]]))
  for (v in covariates) if (is.factor(base[[v]]))
    base[[v]] <- droplevels(base[[v]])

  ## methylation matrix on the complete-case rows
  if (!is.null(pair_rt)) {
    ## PAIR PATH. mp has one row per usable CpG x gene pair, several of which
    ## can share one physical measurement (co-effects) - and, under
    ## collapse_site, one measurement can average several replicate probe
    ## designs of the same CpG (spec 20). Build one column per MEASUREMENT
    ## first, then expand to one column per PAIR ROW, so everything
    ## downstream (missing-value policy, M-values, the per-system sets) works
    ## index-parallel with mp exactly as it always has.
    mm <- pair_rt$measurements
    mm <- mm[mm$measurement_id %in% mp$column, , drop = FALSE]
    src <- if (is.matrix(methylation)) methylation[cc, , drop = FALSE]
           else as.matrix(base[, intersect(unique(unlist(mm$replicates)),
                                           names(base)), drop = FALSE])
    mode(src) <- "numeric"
    BEm <- vapply(seq_len(nrow(mm)), function(i) {
      ph <- intersect(mm$replicates[[i]], colnames(src))
      if (!length(ph)) return(rep(NA_real_, nrow(src)))
      if (length(ph) == 1L) src[, ph]
      ## replicate designs of ONE CpG average on the beta scale; na.rm stays
      ## FALSE so a missing replicate keeps the site visible to the declared
      ## missing_methylation policy instead of quietly narrowing its basis
      else rowMeans(src[, ph, drop = FALSE], na.rm = FALSE)
    }, numeric(nrow(src)))
    colnames(BEm) <- mm$measurement_id
    .navg <- sum(vapply(mm$replicates, length, 1L) > 1L)
    if (.navg > 0)
      cor <- .frame_note(cor, "replicates",
                         sprintf("%d CpG(s) measured by several probe designs",
                                 .navg),
                         "replicate designs averaged into one measurement per site (replicate_policy = collapse_site)")
    BE <- BEm[, mp$column, drop = FALSE]
    colnames(BE) <- mp$column
  } else if (is.matrix(methylation)) {
    ## the map-matching step accepts each probe by its column name OR its bare
    ## probe id, per row - so the extraction key must be chosen per row too. An
    ## all-or-nothing key silently produced all-NA columns for whichever naming
    ## convention was in the minority in a mixed matrix.
    key <- ifelse(mp$column %in% colnames(methylation), mp$column, mp$probe)
    BE <- methylation[cc, match(key, colnames(methylation)), drop = FALSE]
    ## every downstream consumer (sets$columns, the reporters, the aligner)
    ## keys on the map's `column` names. A matrix keyed by bare probe ids was
    ## accepted here but then indexed out of bounds in the pilot run - so
    ## normalise the extracted block to the canonical column names at once.
    colnames(BE) <- mp$column
  } else BE <- as.matrix(base[, mp$column, drop = FALSE])
  mode(BE) <- "numeric"
  ## ---- spec 48: MISSING METHYLATION NEEDS A DECLARED POLICY ---------------
  ## Downstream, a probe carrying any non-finite value is dropped from its set.
  ## That is a defensible engine rule, but it must not be what happens by
  ## DEFAULT and out of sight: one missing value in one sample deletes a whole
  ## CpG, and a gene can lose all of its evidence without the reader being
  ## asked. The decision is taken here, by the caller, and recorded.
  .na_cell <- !is.finite(BE)
  if (any(.na_cell)) {
    .nna <- sum(.na_cell)
    .np <- sum(colSums(.na_cell) > 0)
    .nr <- sum(rowSums(.na_cell) > 0)
    if (missing_methylation == "error")
      stop("`methylation` carries ", .nna, " missing or non-finite value(s), ",
           "affecting ", .np, " probe column(s) and ", .nr, " sample(s).\n",
           "Left alone, every affected probe is dropped from its set, so a ",
           "gene can lose all of its evidence without the report saying so.\n",
           "Declare a policy:\n",
           "  missing_methylation = \"drop_probes\"          ",
           "drop the affected probes; the correction log records which.\n",
           "  missing_methylation = \"common_complete_rows\"  ",
           "keep the probes; analyse the samples complete across all of them.",
           "\nNo analysis was run.", call. = FALSE)
    if (missing_methylation == "drop_probes") {
      okc <- colSums(.na_cell) == 0
      .dropped <- mp$column[!okc]
      cor <- .frame_note(cor, "methylation",
                         sprintf("%d missing value(s) in %d probe column(s)",
                                 .nna, .np),
                         sprintf(paste0("missing_methylation = \"drop_probes\"",
                                        ": %d probe(s) dropped (%s), %d kept"),
                                 .np,
                                 paste(utils::head(.dropped, 8), collapse = ", "),
                                 sum(okc)))
      mp <- mp[okc, , drop = FALSE]; BE <- BE[, okc, drop = FALSE]
    } else {
      keep_r <- rowSums(.na_cell) == 0
      if (sum(keep_r) < 30)
        stop("missing_methylation = \"common_complete_rows\" leaves only ",
             sum(keep_r), " sample(s) complete across all ", ncol(BE),
             " probe(s) - too few to analyse.\nUse missing_methylation = ",
             "\"drop_probes\" to keep the samples and lose the ", .np,
             " affected probe(s) instead.\nNo analysis was run.",
             call. = FALSE)
      BE <- BE[keep_r, , drop = FALSE]
      base <- base[keep_r, , drop = FALSE]
      ## the analysable row set just changed: factor levels and the
      ## exchangeability blocks must be re-derived on it, not inherited.
      if (.nz(chip)) base[[chip]] <- droplevels(factor(base[[chip]]))
      for (v in covariates) if (is.factor(base[[v]]))
        base[[v]] <- droplevels(base[[v]])
      if (length(blockv)) {
        .sz2 <- table(interaction(base[, blockv, drop = FALSE], drop = TRUE))
        if (all(.sz2 < 2L))
          stop("missing_methylation = \"common_complete_rows\" left the ",
               "exchangeability block with no permutable groups.\n",
               "Block: ", paste(blockv, collapse = " x "), " -> ",
               length(.sz2), " group(s), all of size 1 on the ", nrow(base),
               " retained rows.\nContinuing would change the null model to ",
               "unrestricted permutation.\nNo analysis was run.",
               call. = FALSE)
      }
      cor <- .frame_note(cor, "methylation",
                         sprintf("%d missing value(s) affecting %d sample(s)",
                                 .nna, .nr),
                         sprintf(paste0("missing_methylation = \"common_",
                                        "complete_rows\": %d sample(s) ",
                                        "dropped, %d kept, all %d probe(s) ",
                                        "retained"),
                                 sum(!keep_r), sum(keep_r), ncol(BE)))
    }
  }
  ## Dropping EVERY probe column is a different situation from dropping some,
  ## and it has one usual cause: the methylation set being analysed was not
  ## measured on the same people as the covariates. Saying so beats letting the
  ## frame continue with an empty map and fail later inside the aligner.
  if (!nrow(mp))
    stop("every mapped probe column had missing values on the ", nrow(base),
         " complete-case row(s), so no probe is left to analyse. This usually ",
         "means the methylation set and the covariates were measured on ",
         "different people - e.g. analysing one family member's array with ",
         "another's cell-composition and batch covariates. Check that the ",
         "covariates belong to the same array as `methylation`.",
         call. = FALSE)

  ## ---- spec 49: THE METHYLATION SCALE IS DECLARED, NOT GUESSED ------------
  ## "auto" keeps the historical behaviour (values inside [0,1] are treated as
  ## beta). Declaring "beta" makes it checkable: values outside [0,1] are then
  ## an error rather than being clamped to the boundary and called beta.
  if (methylation_scale == "beta") {
    rng <- range(BE, na.rm = TRUE)
    if (rng[1] < 0 || rng[2] > 1)
      stop("methylation_scale = \"beta\" was declared, but values run from ",
           signif(rng[1], 4), " to ", signif(rng[2], 4),
           ", outside [0, 1].\nThese are not beta values. Clamping them to ",
           "the boundary would report something else entirely as beta.\n",
           "Use methylation_scale = \"M\" if they are M-values, or fix the ",
           "input.\nNo analysis was run.", call. = FALSE)
  }
  mv <- if (methylation_scale == "M") list(M = BE, converted = FALSE)
        else .frame_mvalues(BE)
  if (methylation_scale == "beta" && !mv$converted)
    mv <- list(M = log2(pmin(pmax(BE, 1e-4), 1 - 1e-4) /
                        (1 - pmin(pmax(BE, 1e-4), 1 - 1e-4))), converted = TRUE)
  if (mv$converted)
    cor <- .frame_note(cor, "methylation", "beta values detected",
                       "converted to M-values (log2 b/(1-b), eps 1e-4)")
  M <- mv$M

  ## ---- TEST DRIVE 2: design rank audit ------------------------------------
  ff <- stats::as.formula(paste("~", paste(
    c(outcome, covariates, if (.nz(chip)) chip,
      c(mod, mod2)[nzchar(c(mod, mod2))]), collapse = " + ")))
  X <- stats::model.matrix(ff, base)
  nearconst <- colnames(X)[-1][apply(X[, -1, drop = FALSE], 2,
                                     function(z) stats::sd(z) < 1e-10)]
  qr_ <- qr(X)
  if (qr_$rank < ncol(X)) {
    aliased <- colnames(X)[qr_$pivot[(qr_$rank + 1):ncol(X)]]
    keepvars <- covariates
    for (a in aliased) {
      hit <- covariates[vapply(covariates, function(v) grepl(v, a, fixed = TRUE),
                               logical(1))]
      if (length(hit)) keepvars <- setdiff(keepvars, hit)
    }
    ## compositional note: does a dropped column track the sum of the others?
    comp <- ""
    Z <- X[, -1, drop = FALSE]
    for (a in intersect(aliased, colnames(Z))) {
      others <- setdiff(colnames(Z), a)
      s <- rowSums(Z[, others, drop = FALSE])
      cc_ <- if (stats::sd(s) > 0 && stats::sd(Z[, a]) > 0)
        suppressWarnings(stats::cor(Z[, a], s)) else NA_real_
      if (is.finite(cc_) && abs(cc_) > .999)
        comp <- " (compositional group: it mirrors the sum of the others - e.g. cell fractions)"
    }
    if (!autofix)
      stop("design is rank-deficient; aliased: ",
           paste(aliased, collapse = ", "), call. = FALSE)
    cor <- .frame_note(cor, paste(aliased, collapse = ","),
                       "aliased (rank-deficient design)",
                       paste0("dropped from the model", comp))
    covariates <- keepvars
  }
  if (length(nearconst)) {
    if (!autofix) stop("near-constant column(s): ",
                       paste(nearconst, collapse = ", "), call. = FALSE)
    hit <- covariates[vapply(covariates, function(v)
      any(grepl(v, nearconst, fixed = TRUE)), logical(1))]
    if (length(hit)) {
      covariates <- setdiff(covariates, hit)
      cor <- .frame_note(cor, paste(hit, collapse = ","), "near-constant",
                         "dropped from the model")
    }
  }
  Xs <- scale(X[, -1, drop = FALSE])
  Xs <- Xs[, apply(Xs, 2, function(z) all(is.finite(z))), drop = FALSE]
  kappa_ <- if (ncol(Xs)) kappa(crossprod(Xs) / (nrow(Xs) - 1)) else 1

  ## ---- TEST DRIVE 3: random-effects / block audit -------------------------
  blk <- if (length(blockv)) interaction(base[, blockv, drop = FALSE],
                                         drop = TRUE) else
    factor(seq_len(nrow(base)))
  bt <- table(table(blk))
  sizes <- as.integer(names(bt))
  ## permutation count within equal-size strata (couples + singletons lesson)
  logperm <- 0
  for (i in seq_along(bt)) logperm <- logperm + lfactorial(bt[i]) / log(2)
  singles <- if ("1" %in% names(bt)) bt[["1"]] else 0
  if (all(sizes == 1))
    cor <- .frame_note(cor, paste(blockv, collapse = ","),
                       "all blocks are singletons",
                       "reduces to unrestricted permutation")
  lost_strata <- sum(bt == 1 & sizes > 1)
  if (lost_strata > 0)
    cor <- .frame_note(cor, paste(blockv, collapse = ","),
                       sprintf("%d stratum/strata hold a single block", lost_strata),
                       "those blocks cannot move; counted as lost exchangeable units")
  ## block confounded with a fixed covariate?
  for (v in intersect(covariates, names(base))) {
    if (!is.numeric(base[[v]])) next
    wv <- tapply(base[[v]], blk, function(z) stats::var(z))
    if (all(!is.finite(wv) | wv < 1e-12) && stats::sd(base[[v]]) > 0 &&
        nlevels(blk) < nrow(base))
      cor <- .frame_note(cor, v, "constant within every block",
                         "block and covariate are confounded; permutation respects both")
  }
  if (logperm < 10) {
    fmt <- paste0("only ~2^%.1f distinct block permutations are ",
                  "available - too few for stable p-values. Merge ",
                  "blocks or drop the random term.")
    stop(sprintf(fmt, logperm), call. = FALSE)
  }

  ## ---- TEST DRIVE 4: outcome audit ----------------------------------------
  for (oc in outcome) {
    u <- unique(base[[oc]][is.finite(base[[oc]])])
    if (length(u) <= 1) stop("outcome '", oc, "' is constant", call. = FALSE)
    if (length(u) == 2 && outcome_type[[oc]] == "gaussian") {
      if (autofix) {
        ## spec 50: the switch lands on THIS outcome's family, not on a scalar
        ## shared with every other outcome in the same frame.
        outcome_type[[oc]] <- "logistic"
        ## "switched to logistic" reads as though a logistic regression is
        ## fitted. None is. Methylation is the RESPONSE and the outcome is a
        ## predictor of it, so a two-level outcome is a group contrast on a
        ## continuous response - and, in the before/after builds, a
        ## within-person one.
        .k <- .frame_outcome_kind(base, oc, blockv,
                                  .frame_which_build(names(data)))
        cor <- .frame_note(cor, oc, "two-level outcome",
          if (identical(.k$kind, "wave"))
            paste0("treated as a WITHIN-PERSON contrast between the two waves ",
                   "(", oc, " = ", .k$lo, " vs ", .k$hi, "); methylation is the ",
                   "response, so no logistic model is fitted")
          else
            paste0("treated as a group contrast (", oc, " = ", .k$lo, " vs ",
                   .k$hi, ") on a continuous response; no logistic model is ",
                   "fitted"))
      } else stop("outcome '", oc, "' is binary but its declared family is ",
                  "gaussian.\nDeclare it: outcome_type = c(",
                  paste(sprintf("%s = \"%s\"", outcome,
                                ifelse(outcome == oc, "logistic",
                                       outcome_type[outcome])),
                        collapse = ", "), ")", call. = FALSE)
    }
    if (length(u) > 2 && length(u) <= 6 && outcome_type[[oc]] == "gaussian" &&
        all(u == round(u)))
      cor <- .frame_note(cor, oc,
                         sprintf("%d distinct integer levels", length(u)),
                         "kept gaussian and analysed as an ORDERED score - if the levels are unordered this is wrong; recode to a binary contrast or use dmsa_outcome()")
  }

  ## ---- alignment per system under the chosen map --------------------------
  sets <- list()
  for (sid in sys_tab$system_id) {
    s <- mp[mp$system_id == sid, , drop = FALSE]
    ## A selected system with no surviving row used to reach dmsa_align() with a
    ## zero-row table and die inside it on "replacement has 1 row, data has 0",
    ## which says nothing about what the user did wrong.
    if (!nrow(s))
      stop("system '", sys_tab$system[match(sid, sys_tab$system_id)][1],
           "' (id ", sid, ") has no direction-called probe among the ",
           "methylation columns being analysed. ",
           nrow(mp), " mapped probe(s) remain overall, covering system id(s) ",
           paste(unique(mp$system_id), collapse = ", "),
           ". If `methylation` was given explicitly, it may not contain this ",
           "system's columns.",
           if (missing_methylation == "drop_probes")
             paste0(" Note missing_methylation = \"drop_probes\" is in ",
                    "effect: the corrections log says which probes were ",
                    "dropped for missing values - this system's probes may ",
                    "be among them.") else "",
           call. = FALSE)
    al <- dmsa_align(data.frame(cpg = s$probe, d = s$best_direction,
                                p_plus = s$p_plus),
                     genes = s$gene, level = "gene")
    sets[[as.character(sid)]] <-
      list(system_id = sid, system = s$system[1], map = s, alignment = al,
           columns = s$column, genes = unique(s$gene))
  }

  ## ---- map disagreement check (both maps, always) -------------------------
  conflicts <- .frame_conflicts(maps, sys_tab$system_id)
  if (nrow(conflicts))
    warning("direction map disagreement (confidence vs full) for: ",
            paste(sprintf("%s '%s'", conflicts$level, conflicts$unit),
                  collapse = "; "),
            " - analysed under cpg_map = '", cpg_map,
            "'; see frame$map_conflicts", call. = FALSE)

  shared_pairs <- tryCatch(.frame_shared_probes(mp, cs),
                           error = function(e) NULL)
  if (!is.null(shared_pairs) && nrow(shared_pairs))
    cor <- .frame_note(cor, "shared probes",
      sprintf("%d gene pair(s) share probe columns inside a system",
              nrow(shared_pairs)),
      paste0("the same measurement appears under both names, so they are NOT ",
             "independent units in their family: ",
             paste(sprintf("%s (%d)", shared_pairs$pair, shared_pairs$n_shared),
                   collapse = "; ")))

  ## resolve level labels BY VALUE when the user names them (PI, 2026-08-29:
  ## sex coded 1/2 with positional labels can silently flip in a human hand;
  ## c("1" = "Husband", "2" = "Wife") cannot). An unnamed pair keeps the
  ## documented positional rule: sorted level values, low then high.
  .lvl_resolve <- function(.ol, .col, .what) {
    .nm <- names(.ol)                 # as.character() drops names - keep them
    .ol <- as.character(.ol); names(.ol) <- .nm
    if (length(.ol) != 2L)
      stop("each `", .what, "` entry must be exactly two labels", call. = FALSE)
    if (is.null(names(.ol)) || !all(nzchar(names(.ol)))) return(unname(.ol))
    .k <- .frame_outcome_kind(base, .col, blockv,
                              .frame_which_build(names(data)))
    if (is.na(.k$lo))
      stop("'", .col, "' is not a two-level variable, so `", .what,
           "` does not apply to it", call. = FALSE)
    .want <- as.character(c(.k$lo, .k$hi))
    if (!setequal(names(.ol), .want))
      stop("`", .what, "` for '", .col, "' names value(s) ",
           paste(setdiff(names(.ol), .want), collapse = ", "),
           ", but its two levels are ", .want[1], " and ", .want[2],
           ". Name them exactly: c(\"", .want[1], "\" = \"...\", \"",
           .want[2], "\" = \"...\")", call. = FALSE)
    unname(.ol[.want])
  }
  ## moderator level labels: only meaningful for a declared, two-level mod
  if (!is.null(mod_levels)) {
    if (!nzchar(mod))
      stop("`mod_levels` labels the moderator's two levels, but no `mod = ` ",
           "was declared", call. = FALSE)
    mod_levels <- .lvl_resolve(mod_levels, mod, "mod_levels")
  }
  frame <- structure(list(
    data = base, M = M, map = mp, sets = sets, systems = sys_tab,
    modules = mods, selection = cs,
    polarity_table = if (is.null(cs)) NULL else cs$polarity,
    module_evidence = if (is.null(cs)) NULL else
      cs$modules[as.character(cs$modules$system_id) %in%
                   as.character(sys_tab$system_id), , drop = FALSE],
    outcome = outcome, covariates = covariates, chip = chip,
    chip_random = chip_random, chip_effect = chip_effect,
    ## names for the two levels of a two-level outcome, e.g. c("T1", "T4") or
    ## c("female", "male"). Without them the sentence has to fall back to
    ## "time = 0" / "time = 1", which is unambiguous but ugly.
    outcome_levels = if (is.null(outcome_levels)) NULL else
      if (is.list(outcome_levels)) {
        ## named list form (PI, 2026-08-29): per-outcome level labels, so a
        ## multi-outcome frame can label each two-level outcome -
        ## list(pills_past_T1 = c("not taking pills", "taking pills"))
        if (is.null(names(outcome_levels)) ||
            !all(names(outcome_levels) %in% outcome))
          stop("`outcome_levels` given as a list must be named by outcome; ",
               "unknown name(s): ",
               paste(setdiff(names(outcome_levels), outcome), collapse = ", "),
               call. = FALSE)
        stats::setNames(lapply(names(outcome_levels), function(.oc)
          .lvl_resolve(outcome_levels[[.oc]], .oc, "outcome_levels")),
          names(outcome_levels))
      } else {
        if (length(outcome_levels) == 2L && length(outcome) == 1L)
          .lvl_resolve(outcome_levels, outcome, "outcome_levels")
        else {
          .ol <- as.character(outcome_levels)
          if (length(.ol) != 2L)
            stop("`outcome_levels` must be exactly two names for the two ",
                 "levels of a two-level outcome (or a named list of such ",
                 "pairs, one per outcome)", call. = FALSE)
          .ol
        }
      },
    outcome_kind = stats::setNames(lapply(outcome, function(.o)
      .frame_outcome_kind(base, .o, blockv, .frame_which_build(names(data)))),
      outcome),
    shared_probes = shared_pairs, mod_levels = mod_levels,
    ## which INSTALLATION built this frame: a frame object outlives package
    ## reinstalls inside a long R session, and stored fields (outcome_kind,
    ## labels, ...) then lag behind new fixes. dmsa_report() compares this
    ## stamp with the running installation and says so (PI's pills battery
    ## drew a stale linear display exactly this way, 2026-08-29).
    built_stamp = tryCatch(
      utils::packageDescription("dmsa")[["Built"]] %||% NA_character_,
      error = function(e) NA_character_),
    labels = c(.lab_pairs(outcome, outcome_label, "outcome_label"),
               .lab_pairs(covariates, covariate_labels, "covariate_labels"),
               if (nzchar(mod) && !is.null(mod_label))
                 stats::setNames(as.character(mod_label)[1], mod) else NULL,
               if (nzchar(mod2) && !is.null(mod2_label))
                 stats::setNames(as.character(mod2_label)[1], mod2) else NULL),
    block = blk, block_cols = blockv,
    levels = c(system = isTRUE(system), module = isTRUE(module),
               gene = isTRUE(gene), probe = isTRUE(probe)),
    moderation = isTRUE(moderation), mod = mod, mod2 = mod2,
    type = type,
    ## one family per outcome (spec 50). A single-outcome frame stores the bare
    ## string it always stored, so nothing downstream has to special-case the
    ## ordinary run; multi-outcome frames store it named by outcome.
    outcome_type = if (length(outcome_type) == 1L) unname(outcome_type)
                   else outcome_type,
    frame_role = frame_role,
    ## spec 51: what is actually fitted, per path, in words. The main path's
    ## direction is fixed by the estimand; the others follow `frame_role`.
    model_orientation = .orient,
    cpg_map = cpg_map, B = as.integer(B), alpha = alpha, fdr = fdr, seed = seed,
    correction = correction, weighting = weighting, w_floor = w_floor,
    palette = palette,
    plots = plots, tables = tables, summary = summary,
    plot_type = plot_type, table_type = table_type, outdir = outdir,
    gene_models = gene_models, progress = progress, beep = beep,
    corrections = cor, map_conflicts = conflicts, untestable = untestable,
    ## spec 1/15/17: how direction was resolved, and - on the pair path - the
    ## full CpG x gene pair ledger (every pair, used or not, with its reason)
    direction_source = .src, tissue = if (!is.null(pair_rt)) tissue else NULL,
    reference_name = if (!is.null(pair_rt)) {
      if (inherits(.ref, "dmsa_reference")) (.ref$name %||% "reference")
      else "reference" } else NULL,
    pair_ledger = if (!is.null(pair_rt)) pair_rt$ledger else NULL,
    measurements = if (!is.null(pair_rt)) pair_rt$measurements else NULL,
    kappa = kappa_, log2_permutations = as.numeric(logperm),
    block_sizes = bt, pilot = NULL), class = "dmsa_frame")

  ## ---- TEST DRIVE 5: B = 49 pilot on the smallest system ------------------
  smallest <- names(sets)[which.min(vapply(sets, function(s)
    length(s$columns), integer(1)))]
  t0 <- Sys.time()
  pr <- tryCatch(
    .report_gene_level(frame, outcome[1], smallest, B = 49L, quiet = TRUE),
    error = function(e) e)
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(pr, "error"))
    stop("pilot run (B = 49, system ", smallest, ") failed: ",
         conditionMessage(pr), call. = FALSE)
  if (!all(is.finite(pr$p_omnibus[pr$n_probes > 0])))
    stop("pilot run produced non-finite p-values; inspect the frame",
         call. = FALSE)
  units_total <- sum(vapply(sets, function(s) length(unique(s$map$gene)),
                            integer(1)))
  eta_min <- dt / 49 * B * length(outcome) *
    max(1, length(sets)) / max(1, 1) / 60
  frame$pilot <- list(seconds_b49 = dt, eta_minutes = round(eta_min, 1),
                      system = smallest, ok = TRUE)
  frame
}

#' Re-run the frame's test drive
#'
#' @param frame a \code{dmsa_frame}.
#' @return the frame's corrections table, invisibly (printed first).
#' @examples
#' set.seed(1)
#' map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
#'                   probe = c("cg01", "cg02"), column = c("cg01", "cg02"),
#'                   best_direction = c(-1, 1), p_plus = c(.1, .9))
#' dat <- data.frame(anx = rnorm(40), cov1 = rnorm(40), cID = rep(1:20, each = 2),
#'                   cg01 = plogis(rnorm(40)), cg02 = plogis(rnorm(40)))
#' fr <- dmsa_frame(dat, map = map, outcome = "anx", covariates = "cov1",
#'                  B = 49, seed = 1, outdir = tempfile())
#'
#' ## beta values were detected and converted; every such fix is on the record
#' dmsa_test_drive(fr)
#' @export
dmsa_test_drive <- function(frame) {
  stopifnot(inherits(frame, "dmsa_frame"))
  cat("Test drive -", nrow(frame$corrections), "recorded action(s)\n")
  if (nrow(frame$corrections)) print(frame$corrections, row.names = FALSE)
  cat(sprintf("design condition number: %.1f\n", frame$kappa))
  cat(sprintf("block permutations: ~2^%.0f\n", frame$log2_permutations))
  if (!is.null(frame$pilot))
    cat(sprintf("pilot (B = 49, system %s): %.1fs -> full run ETA ~%.0f min\n",
                frame$pilot$system, frame$pilot$seconds_b49,
                frame$pilot$eta_minutes))
  invisible(frame$corrections)
}

#' @export
print.dmsa_frame <- function(x, ...) {
  cat("DMSA frame -", nrow(x$data), "rows,", nrow(x$map),
      "mapped probes,", nrow(x$systems), "system(s)\n")
  lv <- names(x$levels)[x$levels]
  cat("  levels:   ", paste(lv, collapse = " > "), "\n")
  fam <- vapply(x$sets, function(s) length(unique(s$map$gene)), integer(1))
  cat("  families: ", paste(sprintf("%s (%d genes, %d probes)",
      vapply(x$sets, `[[`, "", "system"), fam,
      vapply(x$sets, function(s) length(s$columns), integer(1))),
      collapse = "; "), "\n")
  ## spec 50: with one family per outcome, print the outcome WITH its family
  ## rather than one family for the whole frame.
  .ot <- x$outcome_type
  if (length(.ot) > 1L)
    cat("  outcome(s):", paste(sprintf("%s [%s]", x$outcome,
        .ot[x$outcome]), collapse = ", "),
        if (x$frame_role == "outcome") " [frame_role = outcome: these PREDICT the frame]" else "",
        "\n")
  else
  cat("  outcome(s):", paste(x$outcome, collapse = ", "),
      if (x$frame_role == "outcome") " [frame_role = outcome: these PREDICT the frame]" else "",
      "\n")
  cat("  model:    ", paste(unique(.ot), collapse = "/"), x$type,
      "| outcome enters as a PREDICTOR of methylation",
      if (any(.ot != "gaussian"))
        "\n             (so a binary/categorical outcome is a group contrast on a continuous response - no link function is needed or used)" else "",
      "\n            ",
      if (x$moderation) sprintf("| moderation: frame x %s%s", x$mod,
        if (nzchar(x$mod2)) paste0(" x ", x$mod2) else "") else "",
      "\n")
  ## spec 51: name the model each path fits, rather than leaving the reader to
  ## infer it from `frame_role` - which governs only some of them.
  if (!is.null(x$model_orientation)) {
    .mo <- x$model_orientation
    cat("  fits:      main   ", .mo$main, "\n", sep = "")
    if (!is.na(.mo$shape))
      cat("             shape  ", .mo$shape, "\n", sep = "")
    if (!is.na(.mo$moderation))
      cat("             moder. ", .mo$moderation,
          "  [frame_role = ", x$frame_role, "]\n", sep = "")
  }
  if (!is.null(x$module_evidence) && nrow(x$module_evidence)) {
    cat("  sets:     ", x$selection$name, "\n", sep = "")
    .cas_evidence_banner(x$module_evidence, indent = "  ")
    .cas_polarity_banner(list(polarity = x$polarity_table), indent = "  ")
    ## spec 40: unresolved polarity is said at the frame, per system
    .pa <- tryCatch(.rp_polarity_audit(x), error = function(e) NULL)
    if (!is.null(.pa) && any(.pa$n_polarity_unresolved > 0)) {
      .u <- .pa[.pa$n_polarity_unresolved > 0, , drop = FALSE]
      cat(sprintf(paste0("  polarity unresolved: %s - these genes are weighted 0 ",
                         "in the system score (never assumed +1); the gene ",
                         "level is unaffected\n"),
                  paste(sprintf("%s %d of %d testable gene(s) (%s)", .u$system,
                                .u$n_polarity_unresolved, .u$n_genes_testable,
                                gsub(";", ", ", .u$genes_unresolved)),
                        collapse = "; ")))
    }
  }
  if (!is.null(x$direction_source) && x$direction_source == "cpgdirection")
    cat("  direction: cpgdirection pairs (", x$tissue, "), reference: ",
        x$reference_name, " | ",
        sum(x$pair_ledger$used), " of ", nrow(x$pair_ledger),
        " pair(s) used - dmsa_coverage() itemises\n", sep = "")
  cat("  map:      cpg_map =", x$cpg_map,
      if (nrow(x$map_conflicts)) sprintf("(%d unit(s) DISAGREE between maps - see $map_conflicts)",
                                         nrow(x$map_conflicts)) else "(maps agree)", "\n")
  cat("  engine:   ", (x$weighting %||% "combined"),
      switch(x$weighting %||% "combined",
        combined = "(flat + reliability fused per lens; joint-null maxT)",
        reliability = "(centrality-weighted)", flat = "(equal weights)"), "\n")
  cat("  inference: B =", x$B, "|", x$correction,
      "within level-local families | blocks:",
      if (length(x$block_cols)) paste(x$block_cols, collapse = ", ") else "none",
      sprintf("(~2^%.0f permutations)\n", x$log2_permutations))
  if (!is.null(x$pilot))
    cat(sprintf("  pilot:     B = 49 in %.1fs -> ETA ~%.0f min for the full report\n",
                x$pilot$seconds_b49, x$pilot$eta_minutes))
  if (nrow(x$corrections)) {
    cat("  corrections (test drive):\n")
    for (i in seq_len(nrow(x$corrections)))
      cat(sprintf("   - %s: %s -> %s\n", x$corrections$field[i],
                  x$corrections$issue[i], x$corrections$action[i]))
  }
  invisible(x)
}
