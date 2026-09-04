## ---------------------------------------------------------------------------
## Shared plumbing of the path functions (dmsa_modmed / dmsa_mediate with
## vector roles, dmsa_transmit, dmsa_concordance): the direction map, the
## measurement rule, the tones, the chip transform, the exchangeability
## structure. Everything here is INPUT preparation; the statistical engine
## (mediate.R: .med_measure / .med_moments / .med_reg / .med_fl) is untouched.
##
## Measurement rule v1.2 (PREREG v1.2 section 4c, PI ruling 2026-09-04): every
## module tone is chip quasi-demeaned, then residualised on the technical
## terms of the contract (`measure_covariates`) PLUS the sex column PLUS the
## global mean M-value of the array; sex and the global mean are also fixed
## effects of every structural equation. `sex = NULL` and
## `global_mean = FALSE` switch the two additions off; the output's
## `$measurement_rule` says which rule ran.
## ---------------------------------------------------------------------------

## the global mean M-value: over the methylation matrix when one is given,
## otherwise over every cg-site column of `data` (documented on the output)
.paths_global_mean <- function(data, methylation, beta_input) {
  if (!is.null(methylation)) {
    M <- as.matrix(methylation); mode(M) <- "numeric"
    if (beta_input) M <- dmsa_mvalues(M)
    return(list(value = rowMeans(M, na.rm = TRUE), over = sprintf("%d methylation columns", ncol(M))))
  }
  cg <- names(data)[.map_is_site(names(data))]
  if (!length(cg)) return(NULL)      # nothing to average: the caller carries on without the term
  M <- as.matrix(data[, cg, drop = FALSE]); mode(M) <- "numeric"
  if (beta_input) M <- dmsa_mvalues(M)
  list(value = rowMeans(M, na.rm = TRUE), over = sprintf("%d cg-site columns of data", ncol(M)))
}

## the direction map for the path functions. "auto" = cpgdirection pairs when
## the package is installed (as dmsa_frame()), the bundled snapshot otherwise
## WITH A RECORDED NOTE; "alpha" = the snapshot; a data.frame = a user map.
.paths_map <- function(map, available_columns, tissue = "blood", map_level = "full",
                       reference = "alpha") {
  if (is.data.frame(map)) {
    mp <- .frame_maps(map)[[map_level]]
    return(list(mp = mp, source = "user map", note = NULL))
  }
  if (identical(map, "auto") && is.null(.cpgd("cpg_gene_pairs")) && is.null(getOption("dmsa.pair_table", NULL))) {
    sn <- .paths_snapshot(tissue, map_level)
    if (!is.null(sn)) return(list(mp = sn, source = sprintf("bundled cpgdirection snapshot (coverage_v5_cpgdirection, tissue %s; %d pairs)", tissue, nrow(sn)),
                                  note = "cpgdirection is not installed: the bundled snapshot of its pair directions for the Alpha cascade was used (regenerated before release); install cpgdirection to resolve directions for your own columns"))
    mp <- .frame_maps("alpha")[[map_level]]
    return(list(mp = mp, source = "bundled snapshot (coverage_v4_full)",
                note = "cpgdirection is not installed: falling back to the bundled per-column Alpha snapshot; install cpgdirection for pair-level direction resolution"))
  }
  if (identical(map, "alpha")) {
    mp <- .frame_maps("alpha")[[map_level]]
    return(list(mp = mp, source = "bundled snapshot (coverage_v4_full)", note = NULL))
  }
  if (identical(map, "cpgdirection") && is.null(.cpgd("cpg_gene_pairs")) && is.null(getOption("dmsa.pair_table", NULL)))
    stop("map = \"cpgdirection\" but the cpgdirection package (with cpg_gene_pairs) is not installed; BiocManager::install(\"teindor/cpgdirection\")", call. = FALSE)
  if (!identical(map, "auto") && !identical(map, "cpgdirection"))
    stop("`map` must be \"auto\", \"cpgdirection\", \"alpha\" or a direction-map data.frame", call. = FALSE)
  .ref <- if (inherits(reference, "dmsa_reference")) reference else if (identical(reference, "alpha")) alpha_reference() else
    stop("`reference` must be \"alpha\" or a dmsa_reference object", call. = FALSE)
  rt <- .frame_pairs_runtime(available_columns, .ref, tissue = tissue,
                             pairs = getOption("dmsa.pair_table", NULL))
  mp <- rt[[map_level]]
  list(mp = mp, source = sprintf("cpgdirection pairs (tissue %s; %d usable pairs, %d probes)", tissue, nrow(rt$full), length(unique(rt$full$probe))), note = NULL, ledger = rt$ledger)
}

## the regenerated snapshot of cpgdirection's pair directions for the whole
## Alpha cascade (blood and nasal epithelium), the fallback of the path
## functions when cpgdirection is not installed (PI ruling 3, 2026-09-04)
.paths_snapshot <- function(tissue = "blood", map_level = "full") {
  f <- system.file("extdata", "coverage_v5_cpgdirection.csv.gz", package = "dmsa")
  if (!nzchar(f)) return(NULL)
  sn <- utils::read.csv(gzfile(f), stringsAsFactors = FALSE)
  if (!tissue %in% sn$tissue) tissue <- "blood"
  sn <- sn[sn$tissue == tissue, , drop = FALSE]; sn$system_id <- as.character(sn$system_id)
  if (map_level == "confidence") sn <- sn[sn$direction_tier %in% c("M", "S1"), , drop = FALSE]
  sn
}

## module tones for one system, two builders:
##   "aligned"  - the package rule (dmsa_align x dmsa_scores expected score), as dmsa_modmed() 0.99.12-13
##   "majority" - the v12 tone of the 2026i proofs of concept: the module's
##                majority-direction probes, mean of the scaled M-values,
##                signed by the majority direction; floor min_module_probes
## Both: chip quasi-demeaned, then residualised on CVres (the measurement rule).
.paths_tones <- function(s, ME, mp, cas, sysid, pol_tab, min_module_probes, module_weights, chipv, CVres,
                         module_score = c("aligned", "majority"), anchor = NULL) {
  module_score <- match.arg(module_score)
  if (module_score == "aligned")
    return(.mm_tones(s, ME, mp, cas, sysid, pol_tab, min_module_probes, module_weights, chipv, CVres))
  sm <- mp[mp$system_id == sysid[[s]], , drop = FALSE]
  gm <- unique(cas[cas$system_short == s, c("module_id", "module", "gene"), drop = FALSE])
  sm$module_id <- gm$module_id[match(sm$gene, gm$gene)]
  sm$col <- ifelse(sm$probe %in% colnames(ME), sm$probe, ifelse(sm$column %in% colnames(ME), sm$column, NA_character_))
  sm <- sm[!is.na(sm$col) & !is.na(sm$module_id) & !duplicated(sm$col), , drop = FALSE]
  mods <- unique(sm$module_id)
  if (!is.null(anchor) && !is.null(anchor[[s]])) mods <- c(intersect(anchor[[s]], mods), setdiff(mods, anchor[[s]]))
  tl <- list(); np <- integer(); cov_rows <- list()
  for (md in mods) {
    cm <- sm[sm$module_id == md, , drop = FALSE]; lab <- gm$module[match(md, gm$module_id)]
    d <- cm$best_direction; dmaj <- sign(sum(d, na.rm = TRUE)); if (is.na(dmaj) || dmaj == 0) dmaj <- 1
    cols <- cm$col[is.na(d) | d == dmaj]
    X <- ME[, cols, drop = FALSE]; X <- X[, colSums(is.na(X)) == 0, drop = FALSE]
    used <- ncol(X) >= min_module_probes
    if (used) { tl[[md]] <- dmaj * rowMeans(scale(X)); np <- c(np, nrow(cm)) }
    cov_rows[[md]] <- data.frame(system = s, module = lab, probes = nrow(cm), used = used)
  }
  if (length(tl) < 1L) stop("system '", s, "' has no usable module tone at ", min_module_probes, "+ probes", call. = FALSE)
  TT <- do.call(cbind, tl); colnames(TT) <- names(tl)
  if (!is.null(chipv)) {
    chipv <- droplevels(as.factor(chipv)); nc <- as.numeric(table(chipv)); names(nc) <- levels(chipv)
    gi <- as.character(chipv); X0 <- matrix(1, nrow(TT), 1)
    for (j in seq_len(ncol(TT))) {
      gam <- tryCatch(.ri_reml(TT[, j], X0, chipv), error = function(e) 0)
      if (gam > 0) TT[, j] <- as.numeric(.ri_transform(TT[, j, drop = FALSE], chipv, nc, gi, gam))
    }
  }
  if (!is.null(CVres)) { qrc <- qr(CVres); TT <- apply(TT, 2, function(v) qr.resid(qrc, v)) }
  attr(TT, "mod_w") <- switch(module_weights, equal = NULL, probes = np, sqrt_probes = sqrt(np))
  attr(TT, "n_probes") <- np
  list(tones = TT, coverage = do.call(rbind, cov_rows))
}

## chip quasi-demeaning of one vector ((1 | chip) by REML, as the engine does for tones)
.paths_qd <- function(v, chipv) {
  if (is.null(chipv)) return(as.numeric(v))
  chipv <- droplevels(as.factor(chipv)); nc <- as.numeric(table(chipv)); names(nc) <- levels(chipv)
  gi <- as.character(chipv); X0 <- matrix(1, length(v), 1)
  gam <- tryCatch(.ri_reml(as.numeric(v), X0, chipv), error = function(e) 0)
  if (gam > 0) as.numeric(.ri_transform(cbind(as.numeric(v)), chipv, nc, gi, gam)) else as.numeric(v)
}

## the analysable frame of one contract: complete cases, the M-value matrix,
## the map, the measurement-rule design, the tones of the requested systems,
## the chip-transformed columns and covariates, the exchangeability structure
.paths_setup <- function(data, methylation = NULL, map = "auto", systems = character(), columns = character(),
                         covariates = character(), measure_covariates = NULL, sex = NULL, global_mean = TRUE,
                         chip = NULL, block = NULL, strata = NULL, beta_input = TRUE, tissue = "blood",
                         map_level = "full", min_module_probes = 3L, module_weights = "equal",
                         module_score = "aligned", anchor = NULL, reference = "alpha", extra = character(),
                         need_tones = TRUE, pairs = c("cascade", "all")) {
  pairs <- match.arg(pairs)
  if (inherits(data, "dmsa_import")) { if (is.null(methylation)) methylation <- data$methylation; data <- data$data }
  data <- as.data.frame(data)
  sets <- dmsa_sets("alpha"); cas <- sets$cascade
  allsys <- unique(cas$system_short)
  bad <- setdiff(systems, allsys)
  if (length(bad)) stop("not a system short name: ", paste(bad, collapse = ", "), call. = FALSE)
  ## the global mean is computed on ALL rows before the complete-case cut
  gm <- if (isTRUE(global_mean)) .paths_global_mean(data, methylation, beta_input) else NULL
  if (!is.null(gm)) data[[".global_mean"]] <- gm$value
  if (!is.null(sex) && !sex %in% names(data)) stop("sex column '", sex, "' not in data", call. = FALSE)
  add <- c(if (!is.null(sex)) sex, if (!is.null(gm)) ".global_mean")
  covariates <- unique(c(covariates, add))
  mcv <- unique(c(if (is.null(measure_covariates)) setdiff(covariates, add) else measure_covariates, add))
  need <- unique(c(columns, covariates, mcv, chip, block, strata, extra))
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("column(s) not in data: ", paste(miss, collapse = ", "), call. = FALSE)
  keep <- stats::complete.cases(data[, need, drop = FALSE])
  ME <- NULL; mapinfo <- list(mp = NULL, source = "none", note = NULL)
  if (length(systems)) {
    ME0 <- if (is.null(methylation)) {
      pc <- names(data)[.map_is_site(names(data))]
      if (!length(pc)) stop("no cg-site columns found in data; pass `methylation`", call. = FALSE)
      as.matrix(data[, pc, drop = FALSE])
    } else as.matrix(methylation)
    if (nrow(ME0) != nrow(data)) stop("methylation has ", nrow(ME0), " rows but data has ", nrow(data), call. = FALSE)
    mode(ME0) <- "numeric"
    if (beta_input) ME0 <- dmsa_mvalues(ME0)
    mapinfo <- .paths_map(map, colnames(ME0), tissue = tissue, map_level = map_level, reference = reference)
    mp <- mapinfo$mp; mp$system_id <- as.character(mp$system_id)
    ## pairs = "cascade": only the CpG x gene assignments of the selection
    ## cascade (the curated membership; map v5 of the 2026i runs) - a probe
    ## that pairs with a second gene in the reference is not counted twice.
    ## "all" keeps every usable pair, as dmsa_frame() does (co-effects).
    if (pairs == "cascade" && !is.null(cas$cpg)) {
      k <- paste(mp$probe, mp$gene) %in% paste(cas$cpg, cas$gene)
      mapinfo$source <- sprintf("%s; cascade pairs only (%d of %d)", mapinfo$source, sum(k), nrow(mp))
      mp <- mp[k, , drop = FALSE]
    }
    mapinfo$mp <- mp
    pc <- intersect(unique(c(mp$probe, mp$column)), colnames(ME0))
    ME <- ME0[keep, pc, drop = FALSE]
  }
  data <- data[keep, , drop = FALSE]; n <- nrow(data)
  if (n < 30L) stop("only ", n, " complete rows; this engine is not built for n < 30", call. = FALSE)
  chipv <- if (!is.null(chip)) droplevels(as.factor(data[[chip]])) else NULL
  sysid <- stats::setNames(as.character(sets$systems$system_id), sets$systems$system_short)
  pol_tab <- if (inherits(sets$polarity, "dmsa_polarity")) sets$polarity$polarity else sets$polarity
  CVres <- if (length(mcv)) { CV <- cbind(1, as.matrix(data[, mcv, drop = FALSE])); mode(CV) <- "numeric"
    if (all(is.finite(CV))) CV else NULL } else NULL
  tones <- list(); coverage <- list()
  if (need_tones) for (s in systems) {
    tt <- .paths_tones(s, ME, mapinfo$mp, cas, sysid, pol_tab, min_module_probes, module_weights, chipv, CVres,
                       module_score = module_score, anchor = anchor)
    tones[[s]] <- tt$tones; coverage[[s]] <- tt$coverage
  }
  dat2 <- data
  if (!is.null(chipv)) for (v in unique(c(columns, covariates))) dat2[[v]] <- .paths_qd(dat2[[v]], chipv)
  blockv <- if (!is.null(block)) as.factor(data[[block]]) else NULL
  stratav <- if (!is.null(strata)) as.factor(data[[strata]]) else NULL
  covM <- if (length(covariates)) { M <- as.matrix(dat2[, covariates, drop = FALSE]); mode(M) <- "numeric"; M } else NULL
  rule <- sprintf("tones chip %s, residualised on %s; equations carry %s",
                  if (is.null(chipv)) "untransformed" else "quasi-demeaned",
                  if (length(mcv)) paste(sub("^\\.global_mean$", "global mean", mcv), collapse = " + ") else "nothing",
                  if (length(covariates)) paste(sub("^\\.global_mean$", "global mean", covariates), collapse = " + ") else "no covariates")
  list(data = data, dat2 = dat2, ME = ME, mp = mapinfo$mp, map_source = mapinfo$source, map_note = mapinfo$note,
       cas = cas, sets = sets, sysid = sysid, pol_tab = pol_tab, n = n, keep = keep, chipv = chipv,
       tones = tones, coverage = if (length(coverage)) do.call(rbind, coverage) else NULL,
       covM = covM, covariates = covariates, measure_covariates = mcv, CVres = CVres,
       blocks = blockv, stratav = stratav, rule = rule,
       global_mean = if (is.null(gm)) NULL else gm$over, systems = systems, columns = columns,
       min_module_probes = min_module_probes, module_weights = module_weights, module_score = module_score, anchor = anchor)
}

## ---------------------------------------------------------------------------
## dmsa_contract(): one frame's analysis contract for the two-frame functions
## ---------------------------------------------------------------------------

#' Declare one frame's analysis contract
#'
#' The path functions that join two frames (\code{\link{dmsa_transmit}}) or
#' read one frame in two roles (\code{\link{dmsa_concordance}}) take each
#' frame as a contract: the data, the methylation matrix, the structural
#' covariates (fixed effects of the equations), the technical terms
#' residualised from every score before measurement, the sex column, the
#' chip column and the tissue the direction map is resolved for. Measurement
#' rule v1.2 adds the sex column and the array's global mean M-value to both
#' lists; \code{sex = NULL} or \code{global_mean = FALSE} switch that off.
#'
#' @param data data.frame (or a \code{dmsa_import}): ids, phenotypes,
#'   covariates and - unless \code{methylation} is given - the cg-site columns.
#' @param methylation optional numeric matrix of methylation values (rows as
#'   in \code{data}, columns named by CpG id).
#' @param covariates structural covariates: fixed effects of every equation.
#' @param measure_covariates technical terms residualised from the module,
#'   gene and probe scores before measurement; \code{NULL} = the structural
#'   covariates (as \code{dmsa_modmed}).
#' @param sex the sex column (added to both lists under rule v1.2).
#' @param global_mean add the global mean M-value to both lists (rule v1.2).
#' @param chip the chip / BeadChip column: (1 | chip) by REML quasi-demeaning.
#' @param block exchangeability blocks for the permutation (family id).
#' @param strata permutation strata (E15).
#' @param tissue the tissue the cpgdirection map is resolved for.
#' @param beta_input \code{TRUE} when the methylation values are beta values.
#' @return a \code{dmsa_contract} list.
#' @export
dmsa_contract <- function(data, methylation = NULL, covariates = character(), measure_covariates = NULL,
                          sex = NULL, global_mean = TRUE, chip = NULL, block = NULL, strata = NULL,
                          tissue = "blood", beta_input = TRUE) {
  structure(list(data = data, methylation = methylation, covariates = covariates, measure_covariates = measure_covariates,
                 sex = sex, global_mean = global_mean, chip = chip, block = block, strata = strata, tissue = tissue,
                 beta_input = beta_input), class = "dmsa_contract")
}

.as_contract <- function(x, what) {
  if (inherits(x, "dmsa_contract")) return(x)
  if (is.data.frame(x) || inherits(x, "dmsa_import")) return(dmsa_contract(x))
  stop("`", what, "` must be a dmsa_contract(), a data.frame or a dmsa_import", call. = FALSE)
}
