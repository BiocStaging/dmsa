## ---------------------------------------------------------------------------
## The full-hierarchy path functions (PI rule 2026-09-05): "system to system"
## is system > module > gene > probe against system > module > gene > probe,
## and dyadic concordance is the same ladder. dmsa_transmit() tests the same
## unit on both ends of a parent -> child link; dmsa_concordance() tests the
## same unit in both members of a couple, then what the concordance predicts.
##
## Scores per level, every one under the measurement rule (chip
## quasi-demeaned, residualised on the contract's technical terms + sex +
## global mean): probe = the CpG M-value (observed); gene = the aligned
## dmsa_scores() expected score (composite at composite_rel); module = the
## module tone (composite at composite_rel); system = the ladder
## (congeneric -> tau -> single at composite_rel) of the module tones.
##
## Tests: one shared permutation stream per run, the engine's HC
## errors-in-variables statistic (the pivot of .med_fl, E15) on every unit,
## maxT within the level-local family (all units of the system at that
## level, as dmsa_levels()), and - because the stream is shared - a
## tree-wide maxT over all units of all levels (`p_tree`). The engine
## (mediate.R) is untouched: .hier_fl() re-computes .med_fl()'s statistic
## against an EXTERNAL stream so that the units share it (tested identical).
## ---------------------------------------------------------------------------

## residualise one score under the rule: chip quasi-demean, then CVres
.hier_resid <- function(v, chipv, CVres) {
  v <- .paths_qd(as.numeric(v), chipv)
  if (is.null(CVres)) v - mean(v) else qr.resid(qr(CVres), v)
}

## the four levels for one system in one frame
.hier_scores <- function(st, s, composite_rel = 0.70, levels = c("probe", "gene", "module", "system")) {
  mp <- st$mp; cas <- st$cas; ME <- st$ME
  sm <- mp[mp$system_id == st$sysid[[s]], , drop = FALSE]
  gm <- unique(cas[cas$system_short == s, c("module_id", "module", "gene"), drop = FALSE])
  sm$module_id <- gm$module_id[match(sm$gene, gm$gene)]
  sm$col <- ifelse(sm$probe %in% colnames(ME), sm$probe, ifelse(sm$column %in% colnames(ME), sm$column, NA_character_))
  ## every CpG x gene pair is kept (a CpG can regulate more than one gene);
  ## the probe, module and system levels use each probe once
  sm <- sm[!is.na(sm$col) & !is.na(sm$module_id), , drop = FALSE]
  smu <- sm[!duplicated(sm$col), , drop = FALSE]
  P <- ME[, smu$col, drop = FALSE]; ok <- colSums(is.na(P)) == 0
  P <- P[, ok, drop = FALSE]; smu <- smu[ok, , drop = FALSE]; colnames(P) <- smu$probe; sm <- sm[sm$col %in% smu$col, , drop = FALSE]
  pg <- tapply(sm$gene, sm$probe, function(g) paste(unique(g), collapse = "/"))
  out <- list(); meta <- list(probe_gene = stats::setNames(as.character(pg[smu$probe]), smu$probe), probe_module = stats::setNames(smu$module_id, smu$probe), pairs = sm[, c("probe", "gene", "module_id")])
  if ("probe" %in% levels) out$probe <- apply(P, 2, .hier_resid, chipv = st$chipv, CVres = st$CVres)
  if ("gene" %in% levels) {
    genes <- names(which(table(sm$gene) >= 2L))
    G <- vapply(genes, function(g) {
      k <- sm$gene == g
      al <- tryCatch(suppressMessages(suppressWarnings(dmsa_align(
        data.frame(cpg = sm$probe[k], d = sm$best_direction[k], p_plus = sm$p_plus[k]),
        genes = rep(g, sum(k)), level = "gene", polarity = NULL, missing_polarity = "drop"))), error = function(e) NULL)
      if (is.null(al) || !any(al$usable)) return(rep(NA_real_, nrow(P)))
      sc <- tryCatch(suppressWarnings(dmsa_scores(P[, sm$probe[k], drop = FALSE], al, flavours = "aligned")$aligned), error = function(e) rep(NA_real_, nrow(P)))
      if (all(is.na(sc))) sc else .hier_resid(sc, st$chipv, st$CVres) }, numeric(nrow(P)))
    G <- matrix(G, nrow = nrow(P), dimnames = list(NULL, genes))
    out$gene <- G[, colSums(is.na(G)) == 0, drop = FALSE]
    meta$gene_module <- stats::setNames(gm$module_id[match(colnames(out$gene), gm$gene)], colnames(out$gene))
  }
  TT <- st$tones[[s]]; TT <- TT[, colSums(is.na(TT)) == 0, drop = FALSE]
  if ("module" %in% levels) out$module <- TT
  meta$module_name <- stats::setNames(gm$module[match(colnames(TT), gm$module_id)], colnames(TT))
  ms <- .med_measure(TT, "congeneric", label = s, single_rel = composite_rel, mod_w = attr(TT, "mod_w"))
  if ("system" %in% levels) out$system <- matrix(ms$f, ncol = 1, dimnames = list(NULL, s))
  meta$system_rung <- ms$model; meta$system_omega <- ms$omega; meta$system_V <- ms$V
  meta$n <- c(probe = ncol(P), gene = if (is.null(out$gene)) NA_integer_ else ncol(out$gene), module = ncol(TT), system = 1L)
  list(scores = out, meta = meta)
}

## reliability entering the corrected moments at one level
.hier_rel <- function(level, meta, composite_rel) {
  if (level == "probe") return(1)
  if (level != "system") return(composite_rel)
  if (meta$system_rung == "single") composite_rel else 1 / (1 + meta$system_V)
}

## shared permutation stream: exchange among blocks of equal size (the
## scheme of .med_fl block_mode = "exchange"; blocks NULL = free permutation)
.hier_perms <- function(blocks, n, B, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(blocks)) return(lapply(seq_len(B), function(b) sample.int(n)))
  rws <- split(seq_len(n), blocks); grp <- split(seq_along(rws), lengths(rws))
  lapply(seq_len(B), function(b) { ix <- integer(n)
    for (g in grp) ix[unlist(rws[g], use.names = FALSE)] <- unlist(rws[g[sample.int(length(g))]], use.names = FALSE); ix })
}

## the HC errors-in-variables t of .med_fl() for one term, observed and under
## each permutation of the stream (Freedman-Lane on the reduced model);
## `contrast` = a vector over xn replaces the single term (c'b / sqrt(c'Vc))
.hier_fl <- function(D, S, Sraw, xn, yn, term, n, perms, contrast = NULL) {
  Sxx <- S[xn, xn, drop = FALSE]
  Sxxi <- tryCatch(solve(Sxx), error = function(e) solve(Sxx + diag(1e-6 * mean(diag(Sxx)), nrow(Sxx))))
  So <- Sraw[xn, xn, drop = FALSE]; Ofs <- So - Sxx; Xc <- D[, xn, drop = FALSE]; y <- D[, yn]
  cv <- if (is.null(contrast)) { z <- stats::setNames(numeric(length(xn)), xn); z[term] <- 1; z } else contrast
  stat <- function(yv) {
    Sxy <- as.numeric(crossprod(Xc, yv)) / (n - 1); b <- as.numeric(Sxxi %*% Sxy); e <- as.numeric(yv - Xc %*% b)
    Sc <- Xc * e + matrix(as.numeric(Ofs %*% b), n, length(xn), byrow = TRUE)
    V <- Sxxi %*% crossprod(Sc) %*% Sxxi / (n - 1)^2
    cb <- sum(cv * b); vv <- as.numeric(t(cv) %*% V %*% cv)
    c(cb, cb / sqrt(if (!is.finite(vv) || vv <= 0) 1e-12 else vv)) }
  ob <- stat(y)
  red <- if (is.null(contrast)) setdiff(xn, term) else xn[cv == 0]
  fit <- if (length(red)) as.numeric(D[, red, drop = FALSE] %*% .med_reg(Sraw, red, yn, n)$b) else rep(0, n)
  r <- y - fit
  tn <- vapply(perms, function(ix) { yp <- fit + r[ix]; stat(yp - mean(yp))[2] }, numeric(1))
  list(b = ob[1], t = ob[2], tn = tn)
}
.hier_raw_p <- function(tobs, TN) vapply(seq_along(tobs), function(u) (1 + sum(abs(TN[, u]) >= abs(tobs[u]))) / (nrow(TN) + 1), numeric(1))
.hier_maxT <- function(tobs, TN) { mx <- apply(abs(TN), 1, max, na.rm = TRUE); vapply(abs(tobs), function(t) (1 + sum(mx >= t)) / (length(mx) + 1), numeric(1)) }

## the standardized design of one unit on both ends, with the offsets
.hier_design <- function(y, X, cov, rel, Vy = NULL) {
  n <- length(y); corrected <- is.finite(rel) && rel < 1; V <- if (corrected) (1 - rel) / rel else 0
  ys <- as.numeric(scale(y)); Xs <- scale(X)
  if (corrected) { ys <- ys / sqrt(rel); Xs <- Xs / sqrt(rel) }
  C <- if (is.null(cov)) NULL else scale(cov)          # standardized: the PD guard of .med_moments is scale-aware
  D <- cbind(y = ys, Xs, C); D <- scale(D, scale = FALSE); colnames(D) <- c("y", colnames(X), colnames(cov))
  voff <- if (corrected) stats::setNames(rep(V, ncol(X)), colnames(X)) else numeric(0)
  mm <- .med_moments(D, voffset = voff, products = list(), n = n)
  list(D = D, S = mm$S, Sraw = mm$Sobs, shrink = mm$shrink, n = n, xn = c(colnames(X), colnames(cov)))
}

## ---------------------------------------------------------------------------
## dmsa_transmit(): child_u ~ mother_u + father_u + covariates at every level
## ---------------------------------------------------------------------------

#' Parent-to-child transmission at every level of the hierarchy
#'
#' Tests whether the parents' methylation in a biological system predicts
#' the child's methylation in the same system with the SAME UNIT on both
#' ends at every level: probe, gene, module, system. Per unit the estimand is
#' the paired slope \code{child_u ~ mother_u + father_u + covariates} (the
#' child's structural covariates plus the three arrays' global means), tested
#' with the engine's HC errors-in-variables pivot (E15) on one shared
#' permutation stream (family-block exchange), corrected by maxT within the
#' level-local family (all units of the system at that level, as
#' \code{dmsa_levels()}) and, across the four levels, by a tree-wide maxT
#' (\code{p_tree}). Father beside mother throughout: by the registered rule a
#' father effect of the mother's size or larger argues against intrauterine
#' transmission (genotype, household or the paternal germline remain).
#'
#' @param child,parents \code{\link{dmsa_contract}}s (or data.frames) of the
#'   child frame and the parent frame.
#' @param link the family id column present in both frames.
#' @param parent_role the column of the parent frame that says which parent
#'   a row is; \code{roles} maps its values to mother and father.
#' @param roles named vector \code{c(mother = ..., father = ...)} of
#'   \code{parent_role} values.
#' @param systems system short names (as in \code{dmsa_sets("alpha")}).
#' @param levels the levels to test.
#' @param map \code{"auto"} (cpgdirection pairs when installed, the bundled
#'   snapshot otherwise), \code{"alpha"}, \code{"cpgdirection"} or a map.
#' @param composite_rel the fixed reliability of gene and module composites
#'   and of a system whose ladder ends at the single rung (PI ruling .70).
#' @param B permutations of the shared stream.
#' @param seed seed of the stream.
#' @param min_module_probes,module_weights,module_score,anchor the module
#'   tone builder (\code{module_score = "majority"} is the v12 tone of the
#'   2026i proofs of concept; \code{"aligned"} the package rule).
#' @param map_level \code{"full"} or \code{"confidence"}.
#' @param pairs \code{"cascade"} keeps only the CpG x gene assignments of the
#'   selection cascade (the curated membership; a probe pairing with a second
#'   reference gene is not counted twice); \code{"all"} keeps every usable
#'   pair as \code{dmsa_frame()} does.
#' @return a \code{dmsa_transmit} object: \code{$units} (one row per system x
#'   level x unit: mother and father slopes, raw p, maxT p, tree p),
#'   \code{$named}, \code{$summary}, \code{$measurement}, \code{$design}.
#' @export
dmsa_transmit <- function(child, parents, link = "cID", parent_role = "sex", roles = c(mother = 2, father = 1),
                          systems, levels = c("probe", "gene", "module", "system"), map = "auto",
                          composite_rel = 0.70, B = 1999L, seed = 1L, min_module_probes = 3L,
                          module_weights = c("equal", "sqrt_probes", "probes"), module_score = c("majority", "aligned"),
                          anchor = NULL, map_level = c("full", "confidence"), pairs = c("cascade", "all")) {
  module_weights <- match.arg(module_weights); module_score <- match.arg(module_score); map_level <- match.arg(map_level); pairs <- match.arg(pairs)
  levels <- match.arg(levels, c("probe", "gene", "module", "system"), several.ok = TRUE)
  if (!all(c("mother", "father") %in% names(roles))) stop("roles must name mother and father", call. = FALSE)
  ch <- .as_contract(child, "child"); pa <- .as_contract(parents, "parents")
  stc <- .paths_setup(ch$data, ch$methylation, map = map, systems = systems, covariates = ch$covariates,
                      measure_covariates = ch$measure_covariates, sex = ch$sex, global_mean = ch$global_mean, chip = ch$chip,
                      beta_input = ch$beta_input, tissue = ch$tissue, map_level = map_level, min_module_probes = min_module_probes,
                      module_weights = module_weights, module_score = module_score, anchor = anchor, extra = link, pairs = pairs)
  stp <- .paths_setup(pa$data, pa$methylation, map = map, systems = systems, covariates = pa$covariates,
                      measure_covariates = pa$measure_covariates, sex = pa$sex, global_mean = pa$global_mean, chip = pa$chip,
                      beta_input = pa$beta_input, tissue = pa$tissue, map_level = map_level, min_module_probes = min_module_probes,
                      module_weights = module_weights, module_score = module_score, anchor = anchor, extra = c(link, parent_role), pairs = pairs)
  pk <- paste(stp$data[[link]], stp$data[[parent_role]])
  mi <- match(paste(stc$data[[link]], roles[["mother"]]), pk); fi <- match(paste(stc$data[[link]], roles[["father"]]), pk)
  ok <- !is.na(mi) & !is.na(fi)
  if (sum(ok) < 30L) stop("only ", sum(ok), " complete triads; the engine is not built for n < 30", call. = FALSE)
  K <- which(ok); mi <- mi[ok]; fi <- fi[ok]; n <- length(K); fam <- as.factor(stc$data[[link]][K])
  ## equation covariates: the child's (chip-transformed) + the parents' global means
  cov <- stc$covM[K, , drop = FALSE]
  if (!is.null(stp$data[[".global_mean"]])) {
    gmm <- .paths_qd(stp$data[[".global_mean"]][mi], if (is.null(stc$chipv)) NULL else stc$chipv[K])
    gmf <- .paths_qd(stp$data[[".global_mean"]][fi], if (is.null(stc$chipv)) NULL else stc$chipv[K])
    cov <- cbind(cov, global_mean_mother = gmm, global_mean_father = gmf)
  }
  colnames(cov) <- sub("^\\.global_mean$", "global_mean_child", colnames(cov))
  perms <- .hier_perms(fam, n, B, seed)
  units <- list(); meas <- list(); nfam <- 0L
  for (s in systems) {
    SC <- .hier_scores(stc, s, composite_rel, levels); SP <- .hier_scores(stp, s, composite_rel, levels)
    meas[[s]] <- data.frame(system = s, side = c("child", "parents"),
                            rung = c(SC$meta$system_rung, SP$meta$system_rung), omega = c(SC$meta$system_omega, SP$meta$system_omega),
                            probes = c(SC$meta$n[["probe"]], SP$meta$n[["probe"]]), genes = c(SC$meta$n[["gene"]], SP$meta$n[["gene"]]),
                            modules = c(SC$meta$n[["module"]], SP$meta$n[["module"]]), stringsAsFactors = FALSE)
    TNm_all <- list(); TNf_all <- list(); rows <- list()
    for (L in levels) {
      C <- SC$scores[[L]][K, , drop = FALSE]; Pm <- SP$scores[[L]]; Pf <- SP$scores[[L]]
      un <- intersect(colnames(C), colnames(Pm)); if (!length(un)) next
      rel <- .hier_rel(L, SP$meta, composite_rel)
      res <- lapply(un, function(u) {
        des <- .hier_design(C[, u], cbind(mother = Pm[mi, u], father = Pf[fi, u]), cov, rel)
        m <- .hier_fl(des$D, des$S, des$Sraw, des$xn, "y", "mother", n, perms)
        f <- .hier_fl(des$D, des$S, des$Sraw, des$xn, "y", "father", n, perms)
        list(m = m, f = f, shrink = des$shrink) })
      TNm <- do.call(cbind, lapply(res, function(r) r$m$tn)); TNf <- do.call(cbind, lapply(res, function(r) r$f$tn))
      tm <- vapply(res, function(r) r$m$t, numeric(1)); tf <- vapply(res, function(r) r$f$t, numeric(1))
      rows[[L]] <- data.frame(system = s, level = L, unit = un,
                              gene = if (L == "probe") unname(SP$meta$probe_gene[un]) else if (L == "gene") un else NA_character_,
                              module = if (L == "probe") unname(SP$meta$probe_module[un]) else if (L == "gene") unname(SP$meta$gene_module[un]) else if (L == "module") un else NA_character_,
                              b_mother = vapply(res, function(r) r$m$b, numeric(1)), t_mother = tm, p_mother = .hier_raw_p(tm, TNm), p_mother_maxT = .hier_maxT(tm, TNm),
                              b_father = vapply(res, function(r) r$f$b, numeric(1)), t_father = tf, p_father = .hier_raw_p(tf, TNf), p_father_maxT = .hier_maxT(tf, TNf),
                              shrink = vapply(res, function(r) r$shrink, numeric(1)), stringsAsFactors = FALSE)
      TNm_all[[L]] <- TNm; TNf_all[[L]] <- TNf; nfam <- nfam + 2L
    }
    tb <- do.call(rbind, rows); rownames(tb) <- NULL
    TNm <- do.call(cbind, TNm_all); TNf <- do.call(cbind, TNf_all)
    tb$p_mother_tree <- .hier_maxT(tb$t_mother, TNm); tb$p_father_tree <- .hier_maxT(tb$t_father, TNf)
    units[[s]] <- tb
  }
  units <- do.call(rbind, units); rownames(units) <- NULL
  units$module_name <- NA_character_
  named <- units[units$p_mother_maxT < .05 | units$p_father_maxT < .05, , drop = FALSE]
  summ <- do.call(rbind, lapply(split(units, list(units$system, units$level), drop = TRUE), function(d) data.frame(
    system = d$system[1], level = d$level[1], units = nrow(d),
    named_mother = sum(d$p_mother_maxT < .05), min_p_mother_maxT = min(d$p_mother_maxT),
    named_father = sum(d$p_father_maxT < .05), min_p_father_maxT = min(d$p_father_maxT),
    named_tree = sum(d$p_mother_tree < .05 | d$p_father_tree < .05), stringsAsFactors = FALSE)))
  rownames(summ) <- NULL
  summ <- summ[order(match(summ$system, systems), match(summ$level, c("probe", "gene", "module", "system"))), ]
  structure(list(call = match.call(), units = units, named = named, summary = summ, measurement = do.call(rbind, meas),
                 design = list(n = n, families = nlevels(fam), B = B, seed = seed, pivot = "hc", blocks = "family-exchange",
                               families_tested = nfam, composite_rel = composite_rel, module_score = module_score,
                               map_source = stp$map_source, map_note = stp$map_note,
                               rule_child = stc$rule, rule_parents = stp$rule, covariates = colnames(cov)),
                 systems = systems, levels = levels), class = "dmsa_transmit")
}

#' @export
print.dmsa_transmit <- function(x, ...) {
  d <- x$design
  cat("Parent -> child transmission, full hierarchy (same unit on both ends)\n")
  cat(sprintf("  n = %d triads (%d families)   B = %d   pivot: hc   blocks: family-exchange   maxT within the level-local family; p_tree across levels\n", d$n, d$families, d$B))
  cat(sprintf("  map: %s\n  rule (child): %s\n  rule (parents): %s\n", d$map_source, d$rule_child, d$rule_parents))
  if (!is.null(d$map_note)) cat("  note: ", d$map_note, "\n", sep = "")
  cat(sprintf("  %d level-local families tested (mother and father separately); no cross-family correction beyond p_tree\n", d$families_tested))
  cat("  measurement:\n"); m <- x$measurement
  for (i in seq_len(nrow(m))) cat(sprintf("    %-10s %-8s %s (omega %.2f); probes %d, genes %s, modules %d\n", m$system[i], m$side[i], m$rung[i], m$omega[i], m$probes[i], m$genes[i], m$modules[i]))
  cat("  per level:\n"); s <- x$summary
  for (i in seq_len(nrow(s))) cat(sprintf("    %-10s %-7s units %4d | mother named %d (min maxT p %.3f) | father named %d (min maxT p %.3f) | tree-named %d\n",
                                          s$system[i], s$level[i], s$units[i], s$named_mother[i], s$min_p_mother_maxT[i], s$named_father[i], s$min_p_father_maxT[i], s$named_tree[i]))
  if (nrow(x$named)) {
    cat("  named units (maxT p < .05 within the level):\n"); nm <- x$named
    for (i in seq_len(nrow(nm))) {
      kind <- if (nm$level[i] == "probe") "locus" else if (nm$level[i] == "system") "tone" else nm$level[i]
      who <- if (nm$p_mother_maxT[i] < .05 && nm$p_father_maxT[i] < .05) "both parents" else if (nm$p_mother_maxT[i] < .05) "mother" else "father"
      read <- if (who == "both parents" || (nm$p_father_maxT[i] < .05 && abs(nm$b_father[i]) >= abs(nm$b_mother[i]))) "not intrauterine (genotype / household / paternal line)" else "mother-carried"
      cat(sprintf("    %-10s %-7s %-14s mother %.2f (%.3f) father %.2f (%.3f) tree %.3f/%.3f | %s: %s; %s\n", nm$system[i], nm$level[i], nm$unit[i],
                  nm$b_mother[i], nm$p_mother_maxT[i], nm$b_father[i], nm$p_father_maxT[i], nm$p_mother_tree[i], nm$p_father_tree[i], kind, who, read))
    }
  } else cat("  nothing named at any level\n")
  invisible(x)
}

## ---------------------------------------------------------------------------
## dmsa_concordance(): within-couple r per unit at every level against a
## batch-matched re-pairing, then the outcome surface and the exposure form
## ---------------------------------------------------------------------------

#' Dyadic concordance at every level of the hierarchy, and what it predicts
#'
#' Step 0: for every unit of every level, the within-couple correlation of
#' the two partners' scores against pseudo-couples re-paired WITHIN
#' \code{batch} (the processing plate); the re-pairing is the permutation
#' stream, so maxT within the level-local family and the tree-wide
#' \code{p_tree} come from the same draws. \code{batch} has no default: when
#' the partners of a couple share a processing batch, a null that ignores it
#' names batch as concordance (Alpha 2026i: the oxytocin "tone concordance"
#' was plate). Pass \code{batch = FALSE} to declare that no batch is shared.
#'
#' With \code{outcome}: per unit at every level, the dyadic response
#' surface on stacked actor-partner rows (actor, partner, actor^2, actor x
#' partner, partner^2 + covariates) and the incongruence-line contrast
#' \code{a4 = b3 - b4 + b5} by Freedman-Lane with the HC contrast pivot; the
#' probe level is observed, gene / module / system carry the measurement
#' offsets (MATH 1.6) with the uncorrected surface beside them. With
#' \code{exposure} (a 0/1 column of the partner named by \code{exposed}): the
#' slope change \code{gamma3} of the exposed partner's unit on the other
#' partner's unit (MATH 7.3) per unit at every level, the exposure never
#' permuted.
#'
#' @param data a \code{\link{dmsa_contract}} (or data.frame) of the frame.
#' @param couple the couple id column; \code{partner_by} the column that
#'   distinguishes the two partners, \code{roles = c(a = ..., b = ...)} its
#'   values (named as you like: \code{c(woman = 2, man = 1)}).
#' @param batch the processing-batch column shared within a couple, or
#'   \code{FALSE}.
#' @param outcome an outcome column (of the actor rows).
#' @param exposure a 0/1 exposure column; \code{exposed} names the role it
#'   belongs to.
#' @param NP re-pairings of step 0; \code{B} permutations of the surface and
#'   exposure tests.
#' @inheritParams dmsa_transmit
#' @return a \code{dmsa_concordance} object: \code{$step0}, \code{$surface},
#'   \code{$exposure} (each one row per system x level x unit),
#'   \code{$summary}, \code{$measurement}, \code{$design}.
#' @export
dmsa_concordance <- function(data, couple = "cID", partner_by = "sex", roles = c(woman = 2, man = 1), systems,
                             levels = c("probe", "gene", "module", "system"), batch = NULL, outcome = NULL,
                             exposure = NULL, exposed = names(roles)[1], map = "auto", composite_rel = 0.70,
                             B = 1999L, NP = 2000L, seed = 1L, min_module_probes = 3L,
                             module_weights = c("equal", "sqrt_probes", "probes"), module_score = c("majority", "aligned"),
                             anchor = NULL, map_level = c("full", "confidence"), pairs = c("cascade", "all")) {
  module_weights <- match.arg(module_weights); module_score <- match.arg(module_score); map_level <- match.arg(map_level); pairs <- match.arg(pairs)
  levels <- match.arg(levels, c("probe", "gene", "module", "system"), several.ok = TRUE)
  if (length(roles) != 2L || is.null(names(roles))) stop("roles must be a named vector of the two partner_by values", call. = FALSE)
  if (is.null(batch)) stop("`batch` has no default: name the processing-batch column the partners of a couple share (the re-pairing null is drawn within it), or pass batch = FALSE to declare that no batch is shared", call. = FALSE)
  ct <- .as_contract(data, "data")
  extra <- c(couple, partner_by, if (is.character(batch)) batch)
  st <- .paths_setup(ct$data, ct$methylation, map = map, systems = systems, covariates = ct$covariates,
                     measure_covariates = ct$measure_covariates, sex = ct$sex, global_mean = ct$global_mean, chip = ct$chip,
                     beta_input = ct$beta_input, tissue = ct$tissue, map_level = map_level, min_module_probes = min_module_probes,
                     module_weights = module_weights, module_score = module_score, anchor = anchor, extra = extra, pairs = pairs)
  dd <- st$data; ra <- names(roles)[1]; rb <- names(roles)[2]
  cid <- split(seq_len(nrow(dd)), dd[[couple]]); cid <- cid[lengths(cid) == 2L]
  cid <- cid[vapply(cid, function(i) setequal(dd[[partner_by]][i], roles), logical(1))]
  ia <- vapply(cid, function(i) i[dd[[partner_by]][i] == roles[[1]]], integer(1)); ib <- vapply(cid, function(i) i[dd[[partner_by]][i] == roles[[2]]], integer(1))
  ncp <- length(ia)
  if (ncp < 30L) stop("only ", ncp, " complete couples", call. = FALSE)
  bt <- if (is.character(batch)) { bv <- dd[[batch]]; if (!all(bv[ia] == bv[ib])) warning("partners do not always share `batch`; the re-pairing is drawn within the ", rb, "'s batch", call. = FALSE); as.character(bv[ib]) } else rep("all", ncp)
  set.seed(seed)
  grp <- split(seq_len(ncp), bt)
  RP <- lapply(seq_len(NP), function(k) { ix <- seq_len(ncp); for (g in grp) ix[g] <- g[sample.int(length(g))]; ix })
  SC <- list(); meas <- list(); step0 <- list(); nfam <- c(step0 = 0L, surface = 0L, exposure = 0L)
  for (s in systems) {
    S <- .hier_scores(st, s, composite_rel, levels); SC[[s]] <- S
    meas[[s]] <- data.frame(system = s, rung = S$meta$system_rung, omega = S$meta$system_omega, probes = S$meta$n[["probe"]],
                            genes = S$meta$n[["gene"]], modules = S$meta$n[["module"]], stringsAsFactors = FALSE)
    rows <- list(); RN_all <- list()
    for (L in levels) {
      X <- S$scores[[L]]; if (is.null(X) || !ncol(X)) next
      A <- scale(X[ia, , drop = FALSE]); Bm <- scale(X[ib, , drop = FALSE])
      r0 <- colSums(A * Bm) / (ncp - 1)
      RN <- matrix(vapply(RP, function(ix) colSums(A * Bm[ix, , drop = FALSE]) / (ncp - 1), numeric(ncol(X))), nrow = NP, byrow = TRUE)
      mx <- apply(RN, 1, max)
      rows[[L]] <- data.frame(system = s, level = L, unit = colnames(X),
                              gene = if (L == "probe") unname(S$meta$probe_gene[colnames(X)]) else if (L == "gene") colnames(X) else NA_character_,
                              module = if (L == "probe") unname(S$meta$probe_module[colnames(X)]) else if (L == "gene") unname(S$meta$gene_module[colnames(X)]) else if (L == "module") colnames(X) else NA_character_,
                              r = r0, null_mean = colMeans(RN),
                              p = vapply(seq_along(r0), function(u) (1 + sum(RN[, u] >= r0[u])) / (NP + 1), numeric(1)),
                              p_maxT = vapply(r0, function(r) (1 + sum(mx >= r)) / (NP + 1), numeric(1)), stringsAsFactors = FALSE)
      RN_all[[L]] <- RN; nfam[["step0"]] <- nfam[["step0"]] + 1L
    }
    tb <- do.call(rbind, rows); rownames(tb) <- NULL
    mxt <- apply(do.call(cbind, RN_all), 1, max); tb$p_tree <- vapply(tb$r, function(r) (1 + sum(mxt >= r)) / (NP + 1), numeric(1))
    step0[[s]] <- tb
  }
  step0 <- do.call(rbind, step0); rownames(step0) <- NULL
  ## ---- the outcome surface --------------------------------------------------
  surface <- NULL
  if (!is.null(outcome)) {
    if (!outcome %in% names(dd)) stop("outcome column '", outcome, "' not in data", call. = FALSE)
    yv <- as.numeric(dd[[outcome]]); ka <- !is.na(yv[ia]); kb <- !is.na(yv[ib])
    rows <- c(ia[ka], ib[kb]); part <- c(ib[ka], ia[kb]); fam <- as.factor(c(which(ka), which(kb))); nS <- length(rows)
    if (nS < 30L) stop("only ", nS, " actor rows with the outcome", call. = FALSE)
    chS <- if (is.null(st$chipv)) NULL else st$chipv[rows]
    y <- as.numeric(scale(.paths_qd(yv[rows], chS)))
    covS <- if (is.null(st$covM)) NULL else { M <- as.matrix(st$data[rows, st$covariates, drop = FALSE]); mode(M) <- "numeric"; apply(M, 2, .paths_qd, chipv = chS) }
    if (!is.null(covS)) { covS <- matrix(covS, nrow = nS, dimnames = list(NULL, st$covariates)); colnames(covS) <- sub("^\\.global_mean$", "global_mean", colnames(covS)) }
    permsS <- .hier_perms(fam, nS, B, seed + 1L)
    a4_unit <- function(sc, rel) {
      corrected <- is.finite(rel) && rel < 1; act <- as.numeric(scale(sc[rows])); par <- as.numeric(scale(sc[part]))
      if (corrected) { act <- act / sqrt(rel); par <- par / sqrt(rel) }; Vt <- if (corrected) (1 - rel) / rel else 0
      D <- cbind(act = act, par = par, act2 = act^2, actpar = act * par, par2 = par^2, if (is.null(covS)) NULL else scale(covS), y = y)
      D <- scale(D, scale = FALSE); colnames(D) <- c("act", "par", "act2", "actpar", "par2", colnames(covS), "y")
      mm <- if (corrected) .mm_moments_sq(D, voffset = c(act = Vt, par = Vt), products = list(act2 = c("act", "act"), actpar = c("act", "par"), par2 = c("par", "par")), n = nS)
        else list(S = crossprod(D) / (nS - 1), Sobs = crossprod(D) / (nS - 1), shrink = 1)
      xn <- setdiff(colnames(D), "y"); cv <- stats::setNames(numeric(length(xn)), xn); cv[c("act2", "actpar", "par2")] <- c(1, -1, 1)
      r <- .hier_fl(D, mm$S, mm$Sobs, xn, "y", NULL, nS, permsS, contrast = cv)
      list(a4 = r$b, t = r$t, tn = r$tn, shrink = mm$shrink) }
    out <- list()
    for (s in systems) { S <- SC[[s]]; rows_s <- list(); TN_all <- list()
      for (L in levels) { X <- S$scores[[L]]; if (is.null(X) || !ncol(X)) next
        rel <- .hier_rel(L, S$meta, composite_rel)
        res <- lapply(colnames(X), function(u) a4_unit(X[, u], rel)); TN <- do.call(cbind, lapply(res, `[[`, "tn")); tt <- vapply(res, `[[`, numeric(1), "t")
        tb <- data.frame(system = s, level = L, unit = colnames(X), a4 = vapply(res, `[[`, numeric(1), "a4"), t = tt, p = .hier_raw_p(tt, TN), p_maxT = .hier_maxT(tt, TN),
                         shrink = vapply(res, `[[`, numeric(1), "shrink"), stringsAsFactors = FALSE)
        if (L != "probe") { resU <- lapply(colnames(X), function(u) a4_unit(X[, u], 1)); TNU <- do.call(cbind, lapply(resU, `[[`, "tn")); ttU <- vapply(resU, `[[`, numeric(1), "t")
          tb$a4_uncorrected <- vapply(resU, `[[`, numeric(1), "a4"); tb$p_uncorrected <- .hier_raw_p(ttU, TNU); tb$p_uncorrected_maxT <- .hier_maxT(ttU, TNU)
        } else { tb$a4_uncorrected <- tb$a4; tb$p_uncorrected <- tb$p; tb$p_uncorrected_maxT <- tb$p_maxT }
        rows_s[[L]] <- tb; TN_all[[L]] <- TN; nfam[["surface"]] <- nfam[["surface"]] + 1L }
      tb <- do.call(rbind, rows_s); tb$p_tree <- .hier_maxT(tb$t, do.call(cbind, TN_all)); out[[s]] <- tb }
    surface <- do.call(rbind, out); rownames(surface) <- NULL
    attr(surface, "design") <- list(rows = nS, couples = nlevels(fam), outcome = outcome, covariates = colnames(covS))
  }
  ## ---- the exposure form ----------------------------------------------------
  expo <- NULL
  if (!is.null(exposure)) {
    if (!exposure %in% names(dd)) stop("exposure column '", exposure, "' not in data", call. = FALSE)
    if (!exposed %in% names(roles)) stop("`exposed` must be one of ", paste(names(roles), collapse = ", "), call. = FALSE)
    ie <- if (exposed == ra) ia else ib; io <- if (exposed == ra) ib else ia
    ev <- suppressWarnings(as.numeric(dd[[exposure]][ie])); ev[!ev %in% c(0, 1)] <- NA; kp <- !is.na(ev)
    iep <- ie[kp]; iop <- io[kp]; Pv <- ev[kp]; np <- length(iep)
    if (np < 30L || sum(Pv) < 5L || sum(1 - Pv) < 5L) stop("exposure: ", np, " couples with a 0/1 value (", sum(Pv), " exposed)", call. = FALSE)
    set.seed(seed + 2L); permsP <- lapply(seq_len(B), function(b) sample.int(np))
    cvn <- setdiff(st$covariates, ct$sex)
    covP <- if (length(cvn)) { M <- cbind(as.matrix(st$data[iep, cvn, drop = FALSE]), as.matrix(st$data[iop, cvn, drop = FALSE])); mode(M) <- "numeric"
      colnames(M) <- c(paste0(sub("^\\.global_mean$", "global_mean", cvn), "_", exposed), paste0(sub("^\\.global_mean$", "global_mean", cvn), "_", setdiff(names(roles), exposed))); M } else NULL
    g3_unit <- function(sc, rel) {
      e <- as.numeric(scale(sc[iep])); o <- as.numeric(scale(sc[iop])); corrected <- is.finite(rel) && rel < 1
      if (corrected) { e <- e / sqrt(rel); o <- o / sqrt(rel) }; Vt <- if (corrected) (1 - rel) / rel else 0
      Pc <- Pv - mean(Pv); D <- cbind(E = e, O = o, P = Pc, `O:P` = o * Pc, if (is.null(covP)) NULL else scale(covP)); D <- scale(D, scale = FALSE)
      colnames(D) <- c("E", "O", "P", "O:P", colnames(covP))
      mm <- .med_moments(D, voffset = if (corrected) c(O = Vt) else numeric(0), products = list(`O:P` = c("O", "P")), n = np)
      r <- .hier_fl(D, mm$S, mm$Sobs, c("O", "P", "O:P", colnames(covP)), "E", "O:P", np, permsP); c(r, shrink = mm$shrink) }
    out <- list()
    for (s in systems) { S <- SC[[s]]; rows_s <- list(); TN_all <- list()
      for (L in levels) { X <- S$scores[[L]]; if (is.null(X) || !ncol(X)) next
        rel <- .hier_rel(L, S$meta, composite_rel)
        res <- lapply(colnames(X), function(u) g3_unit(X[, u], rel)); TN <- do.call(cbind, lapply(res, `[[`, "tn")); tt <- vapply(res, `[[`, numeric(1), "t")
        tb <- data.frame(system = s, level = L, unit = colnames(X), gamma3 = vapply(res, `[[`, numeric(1), "b"), t = tt, p = .hier_raw_p(tt, TN), p_maxT = .hier_maxT(tt, TN),
                         module = if (L == "probe") unname(S$meta$probe_module[colnames(X)]) else if (L == "gene") unname(S$meta$gene_module[colnames(X)]) else if (L == "module") colnames(X) else NA_character_,
                         shrink = vapply(res, `[[`, numeric(1), "shrink"), stringsAsFactors = FALSE)
        rows_s[[L]] <- tb; TN_all[[L]] <- TN; nfam[["exposure"]] <- nfam[["exposure"]] + 1L }
      tb <- do.call(rbind, rows_s); tb$p_tree <- .hier_maxT(tb$t, do.call(cbind, TN_all)); out[[s]] <- tb }
    expo <- do.call(rbind, out); rownames(expo) <- NULL
    attr(expo, "design") <- list(couples = np, exposed = sum(Pv), exposure = exposure, exposed_role = exposed, covariates = colnames(covP))
  }
  summ <- do.call(rbind, lapply(split(step0, list(step0$system, step0$level), drop = TRUE), function(d) data.frame(
    system = d$system[1], level = d$level[1], units = nrow(d), max_r = max(d$r), null_at_max = d$null_mean[which.max(d$r)],
    named = sum(d$p_maxT < .05), min_p_maxT = min(d$p_maxT), named_tree = sum(d$p_tree < .05), stringsAsFactors = FALSE)))
  rownames(summ) <- NULL; summ <- summ[order(match(summ$system, systems), match(summ$level, c("probe", "gene", "module", "system"))), ]
  structure(list(call = match.call(), step0 = step0, named = step0[step0$p_maxT < .05, , drop = FALSE], summary = summ,
                 surface = surface, exposure = expo, measurement = do.call(rbind, meas),
                 design = list(couples = ncp, NP = NP, B = B, seed = seed, batch = if (is.character(batch)) batch else "none declared (batch = FALSE)",
                               batches = length(grp), families_tested = nfam, composite_rel = composite_rel, module_score = module_score,
                               map_source = st$map_source, map_note = st$map_note, rule = st$rule, roles = roles),
                 systems = systems, levels = levels), class = "dmsa_concordance")
}

#' @export
print.dmsa_concordance <- function(x, ...) {
  d <- x$design
  cat("Dyadic concordance, full hierarchy (same unit in both partners)\n")
  cat(sprintf("  %d couples   NP = %d re-pairings within %s (%d batches)   maxT within the level-local family; p_tree across levels\n", d$couples, d$NP, d$batch, d$batches))
  cat(sprintf("  map: %s\n  rule: %s\n", d$map_source, d$rule))
  if (!is.null(d$map_note)) cat("  note: ", d$map_note, "\n", sep = "")
  cat("  measurement:\n"); m <- x$measurement
  for (i in seq_len(nrow(m))) cat(sprintf("    %-11s %s (omega %.2f); probes %d, genes %s, modules %d\n", m$system[i], m$rung[i], m$omega[i], m$probes[i], m$genes[i], m$modules[i]))
  cat(sprintf("  step 0 (%d level-local families):\n", d$families_tested[["step0"]])); s <- x$summary
  for (i in seq_len(nrow(s))) cat(sprintf("    %-11s %-7s units %4d | max r %.3f (null %.3f) | named %d (min maxT p %.3f) | tree-named %d\n", s$system[i], s$level[i], s$units[i], s$max_r[i], s$null_at_max[i], s$named[i], s$min_p_maxT[i], s$named_tree[i]))
  if (nrow(x$named)) { cat("  named units:\n"); nm <- x$named
    for (i in seq_len(nrow(nm))) { above <- x$step0[x$step0$system == nm$system[i] & x$step0$level != "probe" & x$step0$p_maxT < .05, ]
      kind <- if (nm$level[i] == "probe") { if (nrow(above)) "locus, with a named unit above it" else "locus (nothing named above it)" } else if (nm$level[i] == "system") "diffuse tone" else paste0(nm$level[i], "-level")
      cat(sprintf("    %-11s %-7s %-14s r %.2f (null %.2f) maxT %.3f tree %.3f | %s\n", nm$system[i], nm$level[i], nm$unit[i], nm$r[i], nm$null_mean[i], nm$p_maxT[i], nm$p_tree[i], kind)) }
  } else cat("  step 0: nothing named at any level\n")
  if (!is.null(x$surface)) { a <- attr(x$surface, "design"); sf <- x$surface
    cat(sprintf("  surface -> %s (%d actor rows, %d couples, B = %d; %d families): a4 on the incongruence line\n", a$outcome, a$rows, a$couples, d$B, d$families_tested[["surface"]]))
    for (s in x$systems) for (L in x$levels) { k <- sf$system == s & sf$level == L; if (!any(k)) next
      cat(sprintf("    %-11s %-7s units %4d | min maxT p %.3f, named %d%s\n", s, L, sum(k), min(sf$p_maxT[k]), sum(sf$p_maxT[k] < .05),
                  if (L != "probe") sprintf(" | uncorrected: min maxT p %.3f, named %d", min(sf$p_uncorrected_maxT[k]), sum(sf$p_uncorrected_maxT[k] < .05)) else "")) }
    nm <- sf[sf$p_maxT < .05, ]; for (i in seq_len(nrow(nm))) cat(sprintf("    named: %s %s %s a4 %.2f maxT %.3f tree %.3f\n", nm$system[i], nm$level[i], nm$unit[i], nm$a4[i], nm$p_maxT[i], nm$p_tree[i])) }
  if (!is.null(x$exposure)) { a <- attr(x$exposure, "design"); ex <- x$exposure
    cat(sprintf("  exposure %s of the %s (%d couples, %d exposed; %d families): gamma3 = slope change on the partner's unit\n", a$exposure, a$exposed_role, a$couples, a$exposed, d$families_tested[["exposure"]]))
    for (s in x$systems) for (L in x$levels) { k <- ex$system == s & ex$level == L; if (!any(k)) next; i <- which(k)[which.min(ex$p[k])]
      cat(sprintf("    %-11s %-7s units %4d | top %s gamma3 %.2f raw p %.4f maxT %.3f | named %d\n", s, L, sum(k), ex$unit[i], ex$gamma3[i], ex$p[i], ex$p_maxT[i], sum(ex$p_maxT[k] < .05))) } }
  invisible(x)
}
