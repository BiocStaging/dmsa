## ---------------------------------------------------------------------------
## dmsa_modmed(): general moderated mediation over any mix of OBSERVED and
## LATENT roles (X, M, W, Y each a methylation system or a data column),
## with the moderator on the a-path, the b-path or both (dual stage), an
## optional latent outcome, an optional curvilinear safeguard, and a
## full-pipeline family-cluster bootstrap with the measurement rung held
## fixed. Built ONLY on the shipped engine's primitives (.med_measure,
## .med_moments, .med_reg, .med_fl, .med_collin, .ri_reml/.ri_transform):
## nothing in mediate.R changes, and with every role a system, w_on = "a"
## and quadratic = FALSE the numbers are those of dmsa_mediate() model 7
## (tested). Math: MATH_v1_general_modmed sections 1-5; gates:
## GATE_REPORT_v3, GATE_REPORT_E14, GATE_REPORT_hc.
## ---------------------------------------------------------------------------

## corrected moments with SQUARES allowed (MATH 1.6). Identical to
## .med_moments() when no product is a square; used only when quadratic = TRUE.
.mm_moments_sq <- function(D, voffset, products, n, alpha = 0) {
  k <- ncol(D); nm <- colnames(D)
  Sobs <- crossprod(D) / (n - 1)
  O <- matrix(0, k, k, dimnames = list(nm, nm))
  V <- stats::setNames(numeric(k), nm); V[names(voffset)] <- voffset
  sig <- function(a, b) Sobs[a, b] - if (a == b) V[[a]] else 0
  for (a in names(voffset)) O[a, a] <- voffset[[a]]
  for (pn in names(products)) {
    ab <- products[[pn]]; a <- ab[1]; b <- ab[2]
    O[pn, pn] <- if (a == b) 4 * sig(a, a) * V[[a]] + 2 * V[[a]]^2 else
      sig(a, a) * V[[b]] + V[[a]] * sig(b, b) + V[[a]] * V[[b]]
    for (qn in setdiff(names(products), pn)) {
      cd <- products[[qn]]
      if (a == b && cd[1] == cd[2]) next               # two squares: unbiased
      if (a == b) {                                     # square vs product
        if (a %in% cd) { oth <- setdiff(cd, a)
          O[pn, qn] <- O[qn, pn] <- 2 * sig(a, oth) * V[[a]] }
      } else if (cd[1] == cd[2]) {
        next                                            # handled from the square's side
      } else {
        shared <- intersect(ab, cd)
        if (length(shared) == 1L) { oth <- c(setdiff(ab, shared), setdiff(cd, shared))
          O[pn, qn] <- O[qn, pn] <- sig(oth[1], oth[2]) * V[[shared]] }
      }
    }
  }
  tf <- 1 - alpha / (n - 1); shrink <- 1
  repeat {
    Sstar <- Sobs - (tf * shrink) * O
    ev <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)$values
    ok_pd <- min(ev) > 1e-6 * max(ev)
    ok_cond <- ok_pd && { evR <- eigen(stats::cov2cor(Sstar), symmetric = TRUE, only.values = TRUE)$values
      min(evR) > 0 && max(evR) / min(evR) <= 1e3 }
    if (ok_cond || shrink < 1e-3) break
    shrink <- shrink * 0.9
  }
  list(S = Sstar, Sobs = Sobs, O = O, shrink = shrink * tf, lambda_triggered = shrink < 1)
}

## module tones for one system: the same steps as dmsa_mediate() (probe
## universe and directions from the analysed map, module membership from
## the cascade, aligned scores, chip quasi-demeaning, composition
## residualisation) - kept as a helper here so that mediate.R is untouched.
.mm_tones <- function(s, ME, mp, cas, sysid, pol_tab, min_module_probes,
                      module_weights, chipv, CVres) {
  sm <- mp[mp$system_id == sysid[[s]], , drop = FALSE]
  ps <- if (is.null(pol_tab)) NULL else {
    pp <- if ("system_id" %in% names(pol_tab))
      pol_tab[as.character(pol_tab$system_id) == sysid[[s]], , drop = FALSE] else pol_tab
    if (nrow(pp)) pp else NULL
  }
  gm <- unique(cas[cas$system_short == s, c("module_id", "module", "gene"), drop = FALSE])
  sm$module_id <- gm$module_id[match(sm$gene, gm$gene)]
  sm$module <- gm$module[match(sm$gene, gm$gene)]
  sm$col <- ifelse(sm$probe %in% colnames(ME), sm$probe,
                   ifelse(sm$column %in% colnames(ME), sm$column, NA_character_))
  sm <- sm[!is.na(sm$col) & !is.na(sm$module_id) & !duplicated(sm$col), , drop = FALSE]
  tl <- list(); np <- integer(); cov_rows <- list()
  for (md in unique(gm$module_id)) {
    cm <- sm[sm$module_id == md, , drop = FALSE]
    lab <- gm$module[match(md, gm$module_id)]
    if (nrow(cm) < min_module_probes) {
      cov_rows[[md]] <- data.frame(system = s, module = lab, probes = nrow(cm), used = FALSE); next }
    al <- if (is.null(ps)) NULL else tryCatch(suppressMessages(suppressWarnings(
      dmsa_align(data.frame(cpg = cm$probe, d = cm$best_direction, p_plus = cm$p_plus),
                 genes = cm$gene, level = "system", polarity = ps, missing_polarity = "drop"))),
      error = function(e) NULL)
    sc <- if (is.null(al) || sum(al$usable) < min_module_probes) NULL else
      tryCatch(suppressWarnings(dmsa_scores(ME[, cm$col[match(al$probe, cm$probe)], drop = FALSE],
                                            al, flavours = "aligned")$aligned), error = function(e) NULL)
    if (!is.null(sc)) { tl[[lab]] <- sc; np <- c(np, nrow(cm)) }
    cov_rows[[md]] <- data.frame(system = s, module = lab, probes = nrow(cm), used = !is.null(sc))
  }
  if (length(tl) < 1L)
    stop("system '", s, "' has no usable module tone at ", min_module_probes, "+ probes", call. = FALSE)
  TT <- do.call(cbind, tl)
  mw <- switch(module_weights, equal = NULL, probes = np, sqrt_probes = sqrt(np))
  if (!is.null(chipv)) {
    nc <- as.numeric(table(chipv)); names(nc) <- levels(chipv)
    gi <- as.character(chipv); X0 <- matrix(1, nrow(TT), 1)
    for (j in seq_len(ncol(TT))) {
      gam <- tryCatch(.ri_reml(TT[, j], X0, chipv), error = function(e) 0)
      if (gam > 0) TT[, j] <- as.numeric(.ri_transform(TT[, j, drop = FALSE], chipv, nc, gi, gam))
    }
  }
  if (!is.null(CVres)) { qrc <- qr(CVres); TT <- apply(TT, 2, function(v) qr.resid(qrc, v)) }
  attr(TT, "mod_w") <- mw
  list(tones = TT, coverage = do.call(rbind, cov_rows))
}

## one measured role: list(f, V, role, model, omega, ...) in the engine's
## unit-latent-variance metric
.mm_measure_role <- function(nm, spec, tones, dat2, measurement, single_rel, composite_rel, n, rung = NULL) {
  if (spec$kind == "system") {
    TT <- tones[[nm]]
    if (spec$role == "composite") {
      ms <- .med_measure(TT, "single", label = nm, single_rel = composite_rel[[nm]], mod_w = attr(TT, "mod_w"))
      ms$role <- sprintf("O(rel %.2f)", composite_rel[[nm]])
    } else {
      ms <- .med_measure(TT, if (is.null(rung)) measurement else rung, label = nm, single_rel = single_rel, mod_w = attr(TT, "mod_w"))
      ms$role <- sprintf("L(%s, omega %.2f)", ms$model, ms$omega)
    }
    ms$kind <- "system"; ms
  } else {
    v <- as.numeric(dat2[[nm]]); f <- v - mean(v)
    if (spec$role == "composite") {
      rel <- composite_rel[[nm]]; s2 <- rel * stats::var(f)
      list(f = f / sqrt(s2), V = (1 - rel) / rel, role = sprintf("O(rel %.2f)", rel), model = "column", omega = rel, p = 1L, determinacy = NA_real_, notes = character(), kind = "column")
    } else list(f = f, V = 0, role = "O", model = "column", omega = NA_real_, p = 1L, determinacy = NA_real_, notes = character(), kind = "column")
  }
}

## structural stage: the general O/L core (MATH sections 1-4)
.mm_core <- function(meas, x, m, w, y, w_on, direct_moderation, quadratic, covM, covn,
                     B, blocks, block_mode, conditional_at, alpha_ssc, pivot = "hc", strata = NULL) {
  n <- length(meas[[y]]$f)
  nm <- c(x, m, if (!is.null(w)) w, y)
  F <- vapply(nm, function(v) meas[[v]]$f - mean(meas[[v]]$f), numeric(n)); colnames(F) <- nm
  voff <- vapply(nm, function(v) meas[[v]]$V, numeric(1)); names(voff) <- nm
  voff <- voff[voff > 0]
  prods <- list(); D <- F
  addcol <- function(D, v, pn) { D <- cbind(D, v); colnames(D)[ncol(D)] <- pn; D }
  if (!is.null(w)) {
    if (w_on %in% c("a", "both") || direct_moderation) {
      pn <- paste0(x, ":", w); D <- addcol(D, F[, x] * F[, w], pn); prods[[pn]] <- c(x, w) }
    if (w_on %in% c("b", "both")) {
      pn <- paste0(m, ":", w); D <- addcol(D, F[, m] * F[, w], pn); prods[[pn]] <- c(m, w) }
  }
  if (quadratic) for (v in unique(c(x, w, m))) {
    pn <- paste0(v, "^2"); D <- addcol(D, F[, v]^2, pn); prods[[pn]] <- c(v, v) }
  if (!is.null(covM)) { C <- scale(as.matrix(covM), scale = FALSE); D <- cbind(D, C)
    colnames(D)[seq(ncol(D) - ncol(C) + 1, ncol(D))] <- covn }
  D <- scale(D, scale = FALSE)
  mm <- if (quadratic) .mm_moments_sq(D, voffset = voff, products = prods, n = n, alpha = alpha_ssc) else
    .med_moments(D, voffset = voff, products = prods, n = n, alpha = alpha_ssc)
  Vof <- function(v) if (v %in% names(voff)) voff[[v]] else 0
  xw <- if (!is.null(w)) paste0(x, ":", w) else character(); mw_ <- if (!is.null(w)) paste0(m, ":", w) else character()
  a_terms <- c(x, if (!is.null(w)) w, if (w_on %in% c("a", "both")) xw)
  y_terms <- c(m, x, if (!is.null(w)) w, if (w_on %in% c("b", "both")) mw_, if (direct_moderation) xw)
  sq_m <- if (quadratic) paste0(unique(c(x, w)), "^2") else character()
  sq_y <- if (quadratic) paste0(unique(c(m, x, w)), "^2") else character()
  mx <- c(a_terms, sq_m, covn); yx <- c(y_terms, sq_y, covn)
  ## the same call order as dmsa_mediate(): a-terms, b[m], the direct terms
  paths <- list()
  for (tm in a_terms) paths[[paste0("a[", tm, "]")]] <-
    .med_fl(D, mm$S, mx, m, tm, n, B, blocks, Vy = Vof(m), block_mode = block_mode, Sobs = mm$Sobs, pivot = pivot, strata = strata)
  paths[[paste0("b[", m, "]")]] <-
    .med_fl(D, mm$S, yx, y, m, n, B, blocks, Vy = Vof(y), block_mode = block_mode, Sobs = mm$Sobs, pivot = pivot, strata = strata)
  for (tm in c(x, if (!is.null(w)) w)) paths[[paste0("c'[", tm, "]")]] <-
    .med_fl(D, mm$S, yx, y, tm, n, B, blocks, Vy = Vof(y), block_mode = block_mode, Sobs = mm$Sobs, pivot = pivot, strata = strata)
  if (w_on %in% c("b", "both")) paths[[paste0("b[", mw_, "]")]] <-
    .med_fl(D, mm$S, yx, y, mw_, n, B, blocks, Vy = Vof(y), block_mode = block_mode, Sobs = mm$Sobs, pivot = pivot, strata = strata)
  if (direct_moderation) paths[[paste0("c'[", xw, "]")]] <-
    .med_fl(D, mm$S, yx, y, xw, n, B, blocks, Vy = Vof(y), block_mode = block_mode, Sobs = mm$Sobs, pivot = pivot, strata = strata)
  if (quadratic) for (tm in sq_m) paths[[paste0("a[", tm, "]")]] <-
    .med_fl(D, mm$S, mx, m, tm, n, B, blocks, Vy = Vof(m), block_mode = block_mode, Sobs = mm$Sobs, pivot = pivot, strata = strata)
  g <- function(k) paths[[k]]$b; gp <- function(k) paths[[k]]$p
  a1 <- g(paste0("a[", x, "]")); b1 <- g(paste0("b[", m, "]"))
  a3 <- if (w_on %in% c("a", "both")) g(paste0("a[", xw, "]")) else 0
  b3 <- if (w_on %in% c("b", "both")) g(paste0("b[", mw_, "]")) else 0
  sdw <- if (!is.null(w)) sqrt(max(mm$S[w, w], 0)) else NA_real_
  theta <- function(Wsd) (a1 + a3 * Wsd * sdw) * (b1 + b3 * Wsd * sdw)
  cond <- if (!is.null(w)) data.frame(w_sd = conditional_at, indirect = theta(conditional_at)) else
    data.frame(w_sd = 0, indirect = a1 * b1)
  p_a1 <- gp(paste0("a[", x, "]")); p_b1 <- gp(paste0("b[", m, "]"))
  idx <- list(IE = list(est = a1 * b1, p_js = max(p_a1, p_b1), kind = "a1*b1 (indirect effect at W = 0)", parts = c(paste0("a[", x, "]"), paste0("b[", m, "]"))))
  if (!is.null(w)) {
    pa <- if (w_on %in% c("a", "both")) gp(paste0("a[", xw, "]")) else NA_real_
    pb <- if (w_on %in% c("b", "both")) gp(paste0("b[", mw_, "]")) else NA_real_
    if (w_on == "a") idx$IMM <- list(est = a3 * b1, p_js = max(pa, p_b1), kind = "index of moderated mediation a3*b1 (first stage)", parts = c(paste0("a[", xw, "]"), paste0("b[", m, "]")))
    if (w_on == "b") idx$IMM <- list(est = a1 * b3, p_js = max(pb, p_a1), kind = "index of moderated mediation a1*b3 (second stage)", parts = c(paste0("a[", x, "]"), paste0("b[", mw_, "]")))
    if (w_on == "both") {
      pj <- c(I_a = max(pa, p_b1), I_b = max(pb, p_a1), I_ab = max(pa, pb))
      ph <- stats::p.adjust(pj, "holm")
      idx$IMM <- list(est = theta(1) - theta(-1), p_js = min(ph),
                      products = c(I_a = a3 * b1, I_b = a1 * b3, I_ab = a3 * b3), p_products = pj, p_holm = ph,
                      kind = "dual stage: delta theta = theta(+1 SD) - theta(-1 SD); test = Holm union over the joint-significance p of I_a, I_b, I_ab")
    }
  }
  for (i in seq_along(idx)) if (is.null(idx[[i]]$p_holm)) idx[[i]]$p_holm <- idx[[i]]$p_js
  list(paths = paths, indices = idx, conditional = cond, moments = mm, n = n,
       terms = list(M_model = mx, Y_model = yx),
       collinearity = list(M_model = .med_collin(mm$S, mx), Y_model = .med_collin(mm$S, yx)),
       estimates = c(a1 = a1, a3 = a3, b1 = b1, b3 = b3, IE = a1 * b1, theta_lo = theta(-1), theta_hi = theta(1), dtheta = theta(1) - theta(-1)))
}

#' General moderated mediation over observed and latent roles
#'
#' The general form of \code{\link{dmsa_mediate}}: each of the predictor(s)
#' \code{x}, the mediator(s) \code{m}, the moderator \code{w} and the outcome
#' \code{y} is EITHER a methylation system (measured as a latent tone by the
#' ladder, or as a fixed-reliability composite) OR a column of \code{data}
#' (observed, or a fixed-reliability questionnaire composite). The moderator
#' can act on the a-path (Hayes model 7), the b-path (model 14) or both
#' (model 58, dual stage), with an optional direct-effect moderation (models
#' 8 / 15 / 59). \code{x} and \code{m} may be vectors: several predictors
#' and several PARALLEL mediators (one M-model per mediator, one Y-model),
#' every X_i -> M_j -> Y product tested by joint significance with Holm over
#' the products, and - with two or more latent mediators - the connected
#' measurement of the mediators (MATH 1.7, \code{theta_mm}). Estimation is
#' the corrected-moments (SAM / Croon) engine of \code{dmsa_mediate}: latent
#' products and, with \code{quadratic = TRUE}, squares receive their
#' fourth-moment offsets; inference is Freedman-Lane per path with the
#' errors-in-variables sandwich pivot and family-block exchange; the
#' dual-stage moderation is tested by the Holm union over the
#' joint-significance p of the three products. With every role a system,
#' \code{w_on = "a"} and \code{quadratic = FALSE} the numbers are exactly
#' those of \code{dmsa_mediate(model = "7")}.
#'
#' Measurement rule v1.2 (the package default since 0.99.14): every module
#' tone is chip quasi-demeaned and residualised on \code{measure_covariates}
#' (the technical terms; \code{NULL} = \code{covariates}) PLUS the
#' \code{sex} column PLUS the array's global mean M-value, and sex and the
#' global mean also enter every equation as covariates. \code{sex = NULL} /
#' \code{global_mean = FALSE} switch the additions off; \code{$rule} says
#' which rule ran.
#'
#' @inheritParams dmsa_mediate
#' @param x,m,w,y Character: a system short name (see \code{dmsa_systems()})
#'   or a column of \code{data}. \code{x} and \code{m} may be vectors;
#'   \code{w = NULL} is plain mediation.
#' @param map \code{"auto"} (default): cpgdirection pair directions when the
#'   package is installed, the bundled Alpha snapshot otherwise (recorded on
#'   the output); \code{"alpha"} forces the snapshot; a data.frame is a user
#'   map.
#' @param roles Optional named character vector over \code{c("x","m","w","y")}
#'   with values \code{"latent"} (a system measured by the ladder; the
#'   default for systems), \code{"composite"} (a system as a single-indicator
#'   composite, or a column, at the fixed reliability \code{composite_rel}),
#'   or \code{"observed"} (a column treated as error-free; the default for
#'   columns). A role applies to every variable of that letter.
#' @param composite_rel Reliability for \code{"composite"} roles: a single
#'   value or a vector named by role (\code{x}, \code{m}, ...). The
#'   registered choice for a methylation composite is .70.
#' @param w_on Where the moderator acts: \code{"a"} (first stage),
#'   \code{"b"} (second stage) or \code{"both"} (dual stage).
#' @param direct_moderation Add \code{x:w} to the Y-model (models 8/15/59).
#' @param quadratic Add \code{x^2}, \code{w^2} to the M-model and
#'   \code{m^2}, \code{x^2}, \code{w^2} to the Y-model, with their offsets.
#'   Gate v3: an omitted curvilinear main effect of the size of the paths
#'   produces a spurious interaction 55\% of the time; with the squares the
#'   level is nominal. An interaction claim should survive this.
#' @param sex the sex column (rule v1.2); \code{global_mean} adds the global
#'   mean M-value (over the \code{methylation} matrix, or over every cg-site
#'   column of \code{data}).
#' @param tissue the tissue the cpgdirection map is resolved for.
#' @param module_score \code{"majority"} (default since 0.99.17: the module's
#'   majority-direction probes, mean of the scaled M-values, signed by that
#'   direction - the v12 tone of the 2026i registered runs) or
#'   \code{"aligned"} (the direction-weighted expected score, the default of
#'   0.99.12-16). The registered runs, and every number of the methods paper,
#'   are \code{"majority"}; \code{anchor} an optional named vector of
#'   the module placed first per system; \code{pairs} \code{"cascade"}
#'   (only the selection cascade's CpG x gene assignments) or \code{"all"}.
#' @param theta_mm \code{"auto"} fits the connected measurement of every
#'   pair of latent mediators (lavaan) and applies theta_MM by the .10 rule;
#'   \code{"none"} skips it.
#' @param sensitivities a character vector among \code{"homoscedastic"}
#'   (pivot), \code{"unrestricted"} (no blocks), \code{"composite_rel"}
#'   (every system as a .70 composite), \code{"quadratic"} (the safeguard),
#'   \code{"global_mean_off"} (rule v1.1): each is re-run at \code{B_sens}
#'   permutations and returned as one table (\code{$sensitivities}).
#' @param boot_rung \code{"fixed"} holds every system's measurement rung at
#'   its full-sample value inside the bootstrap and re-estimates loadings
#'   only; \code{"ladder"} reruns the ladder per draw.
#' @return An object of class \code{dmsa_modmed}: the roles and their
#'   measurement, the paths (a-, b- and direct c-prime terms), the
#'   indices (\code{IE}, \code{IMM}, per X -> M pair when vectors), the
#'   conditional indirect effect over \code{conditional_at}, the
#'   collinearity diagnostics, the moment guard, the connected-measurement
#'   record, the sensitivity table, and (with \code{boot > 0}) percentile
#'   and basic bootstrap intervals.
#' @export
dmsa_modmed <- function(data, methylation = NULL, map = "auto",
                        x, m, w = NULL, y, roles = NULL, composite_rel = 0.70,
                        w_on = c("both", "a", "b"), direct_moderation = FALSE,
                        quadratic = FALSE, covariates = character(),
                        measurement = c("congeneric", "tau", "single"),
                        min_module_probes = 3L, single_rel = 0.70,
                        module_weights = c("equal", "sqrt_probes", "probes"),
                        measure_covariates = NULL, sex = NULL, global_mean = TRUE,
                        chip = NULL, block = NULL,
                        perm_blocks = c("none", "chip"),
                        map_level = c("full", "confidence"), tissue = "blood",
                        module_score = c("majority", "aligned"), anchor = NULL, pairs = c("cascade", "all"),
                        B = 1999L, boot = 0L, boot_rung = c("fixed", "ladder"),
                        alpha_ssc = 0, conditional_at = c(-1, 0, 1),
                        beta_input = TRUE, seed = NULL,
                        pivot = c("hc", "homoscedastic"), strata = NULL,
                        theta_mm = c("auto", "none"), sensitivities = NULL, B_sens = 1999L) {
  measurement <- match.arg(measurement); perm_blocks <- match.arg(perm_blocks); pivot <- match.arg(pivot)
  map_level <- match.arg(map_level); module_weights <- match.arg(module_weights); module_score <- match.arg(module_score); pairs <- match.arg(pairs)
  w_on <- match.arg(w_on); boot_rung <- match.arg(boot_rung); theta_mm <- match.arg(theta_mm)
  .args <- mget(setdiff(names(formals(dmsa_modmed)), "sensitivities"))
  if (is.null(w)) { w_on <- "a"; direct_moderation <- FALSE }   # plain mediation
  if (!is.null(seed)) {
    .old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(.old_seed)) assign(".Random.seed", .old_seed, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  cl <- match.call()
  if (inherits(data, "dmsa_import")) { if (is.null(methylation)) methylation <- data$methylation; data <- data$data }
  data <- as.data.frame(data)
  if (inherits(map, "dmsa_sets")) map <- "alpha"
  if (!(is.data.frame(map) || (is.character(map) && length(map) == 1L && map %in% c("auto", "alpha", "cpgdirection"))))
    stop("`map` must be \"auto\", \"alpha\", \"cpgdirection\" or a direction-map data.frame", call. = FALSE)
  cas <- dmsa_sets("alpha")$cascade

  ## -- roles ----------------------------------------------------------------
  x <- as.character(x); m <- as.character(m); y <- as.character(y); if (!is.null(w)) w <- as.character(w)
  if (length(y) != 1L || (!is.null(w) && length(w) != 1L)) stop("y and w must be single names", call. = FALSE)
  letter <- c(stats::setNames(rep("x", length(x)), x), stats::setNames(rep("m", length(m)), m), if (!is.null(w)) stats::setNames("w", w), stats::setNames("y", y))
  vals <- names(letter)
  if (anyDuplicated(vals)) stop("x, m, w and y must be distinct", call. = FALSE)
  allsys <- unique(cas$system_short)
  spec <- lapply(vals, function(v) {
    r <- letter[[v]]; kind <- if (v %in% allsys) "system" else if (v %in% names(data)) "column" else
      stop("'", v, "' is neither a system short name nor a column of data", call. = FALSE)
    role <- if (!is.null(roles) && r %in% names(roles)) roles[[r]] else if (kind == "system") "latent" else "observed"
    if (!role %in% c("latent", "composite", "observed")) stop("roles must be 'latent', 'composite' or 'observed'", call. = FALSE)
    if (kind == "column" && role == "latent") stop("a data column cannot be 'latent'; use 'composite' with composite_rel", call. = FALSE)
    if (kind == "system" && role == "observed") stop("a system cannot be 'observed' (error-free); use 'latent' or 'composite'", call. = FALSE)
    list(name = v, kind = kind, role = role)
  })
  names(spec) <- vals
  crel <- if (length(composite_rel) == 1L && is.null(names(composite_rel))) stats::setNames(rep(composite_rel, length(vals)), vals) else {
    cr <- stats::setNames(rep(NA_real_, length(vals)), vals)
    for (r in names(composite_rel)) cr[letter == r] <- composite_rel[[r]]
    cr }
  for (v in vals) if (spec[[v]]$role == "composite" && !is.finite(crel[[v]]))
    stop("composite_rel needed for role '", v, "'", call. = FALSE)
  systems <- vals[vapply(spec, function(s) s$kind == "system", logical(1))]
  columns <- vals[vapply(spec, function(s) s$kind == "column", logical(1))]
  if (!is.null(block) && perm_blocks == "chip") stop("`block` and perm_blocks = 'chip' cannot be combined", call. = FALSE)

  ## -- the frame: complete cases, map, rule v1.2, tones, chip transform ----
  st <- .paths_setup(data, methylation, map = map, systems = systems, columns = columns, covariates = covariates,
                     measure_covariates = measure_covariates, sex = sex, global_mean = global_mean, chip = chip, block = block,
                     strata = strata, beta_input = beta_input, tissue = tissue, map_level = map_level,
                     min_module_probes = min_module_probes, module_weights = module_weights, module_score = module_score,
                     anchor = anchor, pairs = pairs)
  n <- st$n; tones <- st$tones; dat2 <- st$dat2; chipv <- st$chipv; blockv <- st$blocks; covM <- st$covM; covariates <- st$covariates
  blocks <- if (!is.null(blockv)) blockv else if (perm_blocks == "chip") {
    if (is.null(chipv)) stop("perm_blocks = 'chip' needs `chip`", call. = FALSE); chipv } else NULL
  bmode <- if (!is.null(blockv)) "exchange" else "within"
  stratav <- st$stratav
  covn <- sub("^\\.global_mean$", "global_mean", covariates); if (!is.null(covM)) colnames(covM) <- covn

  ## -- measurement and structural stage ---------------------------------
  meas <- lapply(vals, function(v) .mm_measure_role(v, spec[[v]], tones, dat2, measurement, single_rel, crel, n))
  names(meas) <- vals
  thm <- NULL
  latm <- m[vapply(m, function(v) spec[[v]]$kind == "system" && spec[[v]]$role == "latent", logical(1))]
  if (theta_mm == "auto" && length(latm) >= 2L) {
    thm <- list(); pr <- utils::combn(latm, 2, simplify = FALSE)
    for (p2 in pr) thm[[paste(p2, collapse = ",")]] <- .mm_theta_connected(tones[[p2[1]]], tones[[p2[2]]], meas[[p2[1]]], meas[[p2[2]]], p2)
  }
  core <- .mm_core_v(meas, x, m, w, y, w_on, direct_moderation, quadratic, covM, covn,
                     B, blocks, bmode, conditional_at, alpha_ssc, pivot = pivot, strata = stratav, theta_mm = thm)

  ## -- full-pipeline bootstrap, rung fixed ---------------------------------
  ci <- NULL
  if (boot > 0L) {
    rungs <- lapply(vals, function(v) if (spec[[v]]$kind == "system" && spec[[v]]$role == "latent" && boot_rung == "fixed") meas[[v]]$model else NULL)
    names(rungs) <- vals
    est_of <- function(cc) c(vapply(cc$paths, function(p) p$b, numeric(1)), cc$estimates)
    obs <- est_of(core)
    BM <- matrix(NA_real_, boot, length(obs), dimnames = list(NULL, names(obs)))
    brows <- if (!is.null(blockv)) split(seq_len(n), blockv) else NULL
    stepped <- 0L
    for (bb in seq_len(boot)) {
      ix <- if (is.null(brows)) sample.int(n, n, replace = TRUE) else
        unlist(brows[sample.int(length(brows), length(brows), replace = TRUE)], use.names = FALSE)
      cc <- tryCatch({
        tb <- lapply(tones, function(TT) { mw <- attr(TT, "mod_w"); T2 <- TT[ix, , drop = FALSE]; attr(T2, "mod_w") <- mw; T2 })
        mb <- lapply(vals, function(v) .mm_measure_role(v, spec[[v]], tb, dat2[ix, , drop = FALSE], measurement, single_rel, crel, length(ix), rung = rungs[[v]]))
        names(mb) <- vals
        stepped <- stepped + any(vapply(vals, function(v) !is.null(rungs[[v]]) && mb[[v]]$model != rungs[[v]], logical(1)))
        .mm_core_v(mb, x, m, w, y, w_on, direct_moderation, quadratic, if (is.null(covM)) NULL else covM[ix, , drop = FALSE], covn,
                   0L, NULL, "within", conditional_at, alpha_ssc, theta_mm = thm)
      }, error = function(e) NULL)
      if (!is.null(cc)) { eb <- est_of(cc); BM[bb, names(eb)] <- eb }
    }
    ok <- stats::complete.cases(BM)
    lo <- apply(BM[ok, , drop = FALSE], 2, stats::quantile, 0.025); hi <- apply(BM[ok, , drop = FALSE], 2, stats::quantile, 0.975)
    ci <- list(percentile = cbind(lower = lo, upper = hi), basic = cbind(lower = 2 * obs - hi, upper = 2 * obs - lo),
               estimate = obs, n_ok = sum(ok), draws = boot, rung = boot_rung, stepped = stepped)
  }

  ## -- the sensitivity battery: the same call, one change each, at B_sens ---
  sens <- NULL
  if (length(sensitivities)) {
    sensitivities <- match.arg(sensitivities, c("homoscedastic", "unrestricted", "composite_rel", "quadratic", "global_mean_off"), several.ok = TRUE)
    base <- .args; base$boot <- 0L; base$B <- B_sens
    rows <- list()
    for (sv in sensitivities) {
      a <- base
      if (sv == "homoscedastic") a$pivot <- "homoscedastic"
      if (sv == "unrestricted") { a$block <- NULL; a$perm_blocks <- "none" }
      if (sv == "composite_rel") { a$roles <- c(x = "composite", m = "composite", w = "composite", y = "composite")[unique(letter)]
        for (v in vals) if (spec[[v]]$kind == "column" && spec[[v]]$role == "observed") a$roles[letter[[v]]] <- "observed"; a$composite_rel <- composite_rel }
      if (sv == "quadratic") a$quadratic <- TRUE
      if (sv == "global_mean_off") a$global_mean <- FALSE
      f <- tryCatch(do.call(dmsa_modmed, a), error = function(e) NULL)
      if (is.null(f)) { rows[[sv]] <- data.frame(model = sv, term = NA_character_, estimate = NA_real_, p = NA_real_, stringsAsFactors = FALSE); next }
      rows[[sv]] <- rbind(data.frame(model = sv, term = names(f$paths), estimate = vapply(f$paths, function(p) p$b, numeric(1)), p = vapply(f$paths, function(p) p$p, numeric(1)), stringsAsFactors = FALSE),
                          data.frame(model = sv, term = names(f$indices), estimate = vapply(f$indices, function(i) i$est, numeric(1)), p = vapply(f$indices, function(i) i$p_js, numeric(1)), stringsAsFactors = FALSE))
    }
    sens <- do.call(rbind, rows); rownames(sens) <- NULL; attr(sens, "B") <- B_sens
  }

  structure(list(
    call = cl, n = n, B = B, x = x, m = m, w = w, y = y, w_on = w_on,
    direct_moderation = direct_moderation, quadratic = quadratic,
    roles = stats::setNames(vapply(vals, function(v) meas[[v]]$role, character(1)), vals),
    measurement = stats::setNames(lapply(meas, function(ms) ms[intersect(c("kind", "role", "model", "p", "V", "omega", "determinacy", "lambda", "theta", "notes"), names(ms))]), vals),
    coverage = st$coverage, rule = st$rule, map_source = st$map_source, map_note = st$map_note, global_mean = st$global_mean,
    paths = core$paths, indices = core$indices, conditional = core$conditional,
    estimates = core$estimates, terms = core$terms, ci = ci,
    collinearity = core$collinearity, theta_mm = thm, theta_applied = core$theta_applied,
    sensitivities = sens,
    perm_blocks = if (!is.null(blockv)) "family-exchange" else perm_blocks,
    pivot = pivot, strata = strata,
    ## the corrected SD of every modeled column (sqrt of the diagonal of the
    ## corrected moment matrix; a latent role is 1 by construction). Kept so
    ## that a path can be reported on the standardized scale as well as on
    ## the raw metric of the observed variables - the registered runs and the
    ## methods paper print the standardized one (0.99.17).
    scales = sqrt(pmax(diag(core$moments$S), 0)),
    moments = list(shrink = core$moments$shrink, lambda_triggered = core$moments$lambda_triggered)),
    class = "dmsa_modmed")
}

#' @export
print.dmsa_modmed <- function(x, ...) {
  stage <- switch(x$w_on, a = "first stage (Hayes 7)", b = "second stage (Hayes 14)", both = "dual stage (Hayes 58)")
  if (is.null(x$w)) stage <- "plain mediation"
  cat("General O/L moderated mediation - ", stage, if (x$direct_moderation) " + direct-effect moderation" else "", if (x$quadratic) " + quadratic safeguard" else "", "\n", sep = "")
  cat(sprintf("  n = %d   B = %d permutations   blocks: %s   pivot: %s%s\n", x$n, x$B, x$perm_blocks, x$pivot, if (!is.null(x$strata)) paste0("   strata: ", x$strata) else ""))
  if (!is.null(x$rule)) cat("  rule: ", x$rule, "\n  map:  ", x$map_source, if (!is.null(x$map_note)) paste0(" (", x$map_note, ")") else "", "\n", sep = "")
  cat("  roles:\n")
  for (v in names(x$roles)) {
    ms <- x$measurement[[v]]
    cat(sprintf("    %-14s %-22s%s\n", v, x$roles[[v]], if (ms$kind == "system") sprintf("  %d modules%s", ms$p, if (length(ms$notes)) "  [stepped]" else "") else ""))
  }
  if (x$moments$lambda_triggered)
    cat(sprintf("  NOTE: offset shrunk to %.2f of full correction to keep the moment matrix PD\n", x$moments$shrink))
  for (mdl in names(x$collinearity)) {
    v <- x$collinearity[[mdl]]$vif; if (is.null(v)) next
    flag <- if (max(v) > 10 || x$collinearity[[mdl]]$kappa > 30) "  <-- COLLINEARITY WARNING" else ""
    cat(sprintf("  collinearity %-8s max VIF %.1f (%s)  condition no. %.1f%s\n", mdl, max(v), names(v)[which.max(v)], x$collinearity[[mdl]]$kappa, flag))
  }
  cat("  structural paths (Freedman-Lane p):\n")
  for (nm in names(x$paths)) { p <- x$paths[[nm]]; cat(sprintf("    %-20s b = %+.3f (SE %.3f)  p = %.4f\n", nm, p$b, p$se, p$p)) }
  cat("  indices:\n")
  for (nm in names(x$indices)) { i <- x$indices[[nm]]
    cat(sprintf("    %-6s est = %+.4f  p = %.4f   %s\n", nm, i$est, i$p_js, i$kind))
    if (!is.null(i$products)) for (k in names(i$products))
      cat(sprintf("      %-5s %+.4f  p_js = %.4f  p_holm = %.4f\n", k, i$products[[k]], i$p_products[[k]], i$p_holm[[k]])) }
  if (!is.null(x$theta_mm)) for (k in names(x$theta_mm)) { th <- x$theta_mm[[k]]
    cat(sprintf("  connected measurement %s: theta %.3f (ratio %.2f, phi %.2f) - %s\n", k, th$theta, th$ratio, th$phi, th$note)) }
  if (!is.null(x$w)) { cat("  conditional indirect effect:\n"); cd <- x$conditional; for (i in seq_len(nrow(cd)))
    cat(sprintf("    %sW = %+.1f SD   theta = %+.4f\n", if ("x" %in% names(cd)) sprintf("%s -> %s   ", cd$x[i], cd$m[i]) else "", cd$w_sd[i], cd$indirect[i])) }
  if (!is.null(x$sensitivities)) { cat(sprintf("  sensitivities (B = %d):\n", attr(x$sensitivities, "B"))); sn <- x$sensitivities
    for (md in unique(sn$model)) { k <- sn$model == md & !is.na(sn$term) & grepl("^(IE|IMM)", sn$term)
      cat(sprintf("    %-16s %s\n", md, if (any(k)) paste(sprintf("%s %+.3f (p %.3f)", sn$term[k], sn$estimate[k], sn$p[k]), collapse = "; ") else "failed")) } }
  if (!is.null(x$ci))
    cat(sprintf("  bootstrap: %d/%d full-pipeline resamples usable, rung %s%s (percentile and basic CIs stored)\n",
                x$ci$n_ok, x$ci$draws, x$ci$rung, if (x$ci$stepped) sprintf(", %d draws stepped down a rung", x$ci$stepped) else ""))
  invisible(x)
}
