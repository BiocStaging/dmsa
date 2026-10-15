# ============================================================================
# dmsa_triangulate(): THREE COMPLEMENTARY LENSES, ONE FAMILY, ONE PERMUTATION
#
# WHY THIS EXISTS. AVP's p-value appeared to wander from .0095 to .24 across
# analyses. It did not. AVP's EVIDENCE never changed once: eight CpGs, all
# carrying d = -1, all slopes the same sign, five of eight below .04, and a
# gene-level marginal p of .010 in every single analysis that was run. What
# changed was the multiplicity bill, and it changed because three different
# things were being varied at once:
#
#   1. WHICH STATISTIC   pooled aligned burden, composite score, or quadratic
#   2. WHICH FAMILY      42 genes in one system, or 110-150 across nineteen
#   3. WHICH CORRECTION  maxT on the raw statistic, or minP on calibrated p
#
# Those are three different questions, not three answers to one question. Fix
# them in advance and there is exactly one number per lens. That is what this
# function does.
#
# THE THREE LENSES are chosen to fail in different ways, so agreement between
# them means something:
#
#   coherence  the pooled aligned burden. sum_j m_j b_j / sqrt(sum_j m_j^2 w_j).
#              Rewards CROSS-PROBE SIGN CONCORDANCE, which is the whole premise
#              of a directional method. Blind to a signal whose probes disagree.
#              This is the lens that found AVP.
#
#   composite  the unit's aligned score as a single predictor in one ordinary
#              model. Rewards a UNIFIED SHIFT of the unit that survives as one
#              regressor, and it is the only lens that carries to interactions,
#              categorical outcomes and mediation. Discards concordance, because
#              averaging correlated probes gains far less than sqrt(K). This is
#              the lens that found PLAGL1.
#
#   diffuse    globaltest's quadratic score, weighted by DIRECTION CONFIDENCE:
#                 Q = sum_j |m_j| (r_j' x)^2 / (sum_j |m_j| * x'x)
#              Directionless in the aggregation - it detects signal spread over
#              many probes whichever way each one points - but it spends its
#              weight on the probes whose direction call is trustworthy, which
#              plain globaltest does not. Rewards DIFFUSE signal. This is the
#              family of statistic that found the serotonin and avoidance
#              signals no other engine saw.
#
# All three share one permutation matrix, so their p-values are comparable and
# an omnibus over them is exact. The omnibus is ACAT, which is valid under
# arbitrary dependence and therefore costs no multiplicity for hedging.
# ============================================================================

#' Three complementary lenses on the same family
#'
#' @param M Numeric matrix, samples x probes, methylation on the analysis scale
#'   (probes are the responses).
#' @param data data.frame of predictors, rows aligned to \code{M}.
#' @param rhs Right-hand-side terms, e.g. \code{c("anx","avo","sex_c","chip_f")}
#'   or with an interaction \code{c("anx*ace", ...)}.
#' @param term Design column to test, as \code{model.matrix} names it.
#' @param units Character vector, one entry per column of \code{M}, naming the
#'   unit that column belongs to (typically the gene).
#' @param alignment \code{dmsa_align()} result covering \code{M}'s columns.
#' @param block Permutation blocks: a column name of \code{data} (or several,
#'   combined by interaction), or a vector of labels, one per row. Rows in the
#'   same block travel together under permutation.
#' @param id Optional participant identifier - a column name or per-row
#'   vector. Purely a safety declaration: if any participant contributes more
#'   than one row, \code{block} must keep those rows together, and this stops
#'   the analysis when it would not.
#' @param B Permutations.
#' @param winsor MADs for winsorisation; \code{NULL} to disable.
#' @param method Pooling method for \code{dmsa_test()}.
#' @param correction \code{"maxT"} (default; keeps the information advantage of
#'   a unit with many concordant probes) or \code{"minP"}.
#' @param seed Optional integer.
#' @param weighting Character. Probe weighting engine within a unit: \code{"combined"} (default) fuses the flat and reliability statistics on one shared permutation stream, \code{"flat"} weights every usable aligned probe equally, \code{"reliability"} weights each probe by its item-rest correlation with the rest of its unit. Weights are computed from methylation alone, so the permutation null is unaffected.
#' @param w_floor Numeric. Lower bound applied to reliability weights before normalisation, so a single poorly behaved probe cannot be driven to zero influence. Ignored when \code{weighting = "flat"}.
#' @param ri_group Optional grouping factor for a one-way random intercept,
#'   one entry per row of \code{M} (typically the chip). When supplied with
#'   more than one level, the aligned response and every design column are
#'   quasi-demeaned by a REML-estimated shrinkage factor before any
#'   permutation - the GLS solution for \code{(1 | group)}, at a fraction of
#'   the cost of a mixed fit per probe. \code{NULL} (default) leaves the
#'   design untouched. The transform is applied identically under permutation,
#'   so it moves power, not type-I error.
#' @return data.frame, one row per unit: the raw and family-adjusted p under
#'   each lens, the ACAT omnibus, the direction, and how many probes agreed.
#'   \code{p_union_exact} is the second-level Westfall-Young minP calibration
#'   of the any-lens rule: an EXACT family-wise p for the claim "some lens
#'   carries this unit", computed by ranking the unit's smallest per-lens
#'   family-adjusted p against the null distribution of the family-wide best
#'   adjusted p (one value per permutation draw, min over units and lenses;
#'   returned as \code{attr(., "union_null_min")}). It adapts to the measured
#'   lens dependence - perfectly correlated lenses cost nothing extra,
#'   independent ones the full threefold - and \code{mean(union_null_min <
#'   alpha)} is the realized family-wise error of naming at the per-lens
#'   adjusted-p threshold \code{alpha}. References: Westfall & Young (1993);
#'   Ge, Dudoit & Speed (2003). \code{p_omnibus_adj} family-corrects the
#'   ACAT omnibus the same way (permutation minP over each draw's per-unit
#'   omnibus, \code{attr(., "omnibus_null_min")}): an exact family-wise p
#'   that defends against BOTH multiplicities - the cross-lens combination
#'   is absorbed by ACAT (Liu & Xie 2020), the cross-unit search by the
#'   calibration - and is most powerful when all three lenses agree.
#'   \code{q_unit} is the Benjamini-Hochberg q-value of \code{p_unit} within
#'   the family, and \code{q_unit_dep} its dependence-calibrated version: the
#'   share of permutation draws in which BH, run at the level needed to name
#'   this unit, would have named something in a family without signal
#'   (\code{attr(., "unit_null_qmin")} holds each draw's smallest BH-adjusted
#'   p). Naming at \code{q_unit_dep <= fdr} is the package's naming rule
#'   (PI ruling 2026-10-15); it fires in a fraction \code{fdr} of null
#'   families under the family's measured lens and unit dependence.
#'   References: Benjamini & Hochberg (1995); Yekutieli & Benjamini (1999).
#'   \code{attr(., "probe_fits")} is the per-probe fit the lenses were built
#'   from (\code{column}, \code{b}, \code{se} of each standardised,
#'   covariate-adjusted, random-intercept-transformed probe on the tested
#'   term; its own attributes carry the residual df and the RI gamma), so a
#'   figure can draw exactly the numbers the test used.
#' @examples
#' set.seed(1)
#' n <- 100
#' ## half the probes raise expression when methylated, half lower it: on the
#' ## methylation scale they cancel, on the expression scale they agree
#' d <- rep(c(1, -1), each = 3)
#' f <- rnorm(n)
#' M <- sapply(d, function(dj) dj * .6 * f + sqrt(1 - .6^2) * rnorm(n))
#' colnames(M) <- paste0("cg", 1:6)
#' dat <- data.frame(y = .5 * f + rnorm(n), cv = rnorm(n))
#' al <- dmsa_align(data.frame(cpg = colnames(M), d = d,
#'                             p_plus = ifelse(d > 0, .9, .1)),
#'                  genes = rep("OXTR", 6), level = "gene")
#'
#' res <- dmsa_triangulate(M, dat, rhs = c("y", "cv"), term = "y",
#'                         units = rep("OXTR", 6), alignment = al, B = 99, seed = 1)
#' res[, c("unit", "n_probes", "concordance", "direction", "p_omnibus")]
#' @export
dmsa_triangulate <- function(M, data, rhs, term, units, alignment, block = NULL,
                             id = NULL, B = 1999, winsor = 3,
                             method = c("expected", "fixed"),
                             correction = c("maxT", "minP"),
                             weighting = c("combined", "reliability", "flat"),
                             w_floor = 1.5, seed = NULL, ri_group = NULL) {
  method <- match.arg(method); correction <- match.arg(correction)
  weighting <- match.arg(weighting)
  if (!is.null(seed)) {
    ## restore the caller's RNG state on exit: a permutation seed is for
    ## reproducing THIS result, not for silently reseeding the user's session.
    .old_seed <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(.old_seed))
      assign(".Random.seed", .old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv()))
        ## E10: a true restore of "no prior RNG state" - otherwise a fresh
        ## session leaves deterministically seeded mid-stream of `seed`
        rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  M <- as.matrix(M); data <- as.data.frame(data)
  block    <- .dmsa_rows(block,    data, "block")
  ri_group <- .dmsa_rows(ri_group, data, "ri_group")
  id       <- .dmsa_rows(id,       data, "id")
  .dmsa_check_block(block)
  ## the participant declaration: if anyone contributes more than one row,
  ## those rows are correlated, and a permutation that separates them is
  ## anticonservative. Measured on two-wave data with ICC ~ .6 the true-null
  ## rejection rate roughly triples. So repeats REQUIRE a block that keeps
  ## each participant's rows together.
  if (!is.null(id)) {
    per <- table(id)
    if (max(per) > 1L) {
      if (is.null(block))
        stop(sum(per > 1L), " participant(s) contribute more than one row ",
             "(max ", max(per), "), but block = NULL would permute their ",
             "correlated rows independently. Set block to the participant ",
             "id (block = id), or to a grouping at least as coarse.",
             call. = FALSE)
      spans <- tapply(as.character(block), id, function(z) length(unique(z)))
      if (any(spans > 1L))
        stop(sum(spans > 1L), " participant id(s) span more than one ",
             "permutation block, so a permutation can still separate a ",
             "participant's rows. Each id must sit inside one block.",
             call. = FALSE)
    }
  }
  al <- as.data.frame(alignment)
  if (length(units) != ncol(M) || nrow(al) != ncol(M))
    stop("units, alignment and M must describe the same probes in order",
         call. = FALSE)
  ## E9: same positional-join guard as dmsa_scores()
  if (!is.null(al$probe)) {
    .ord <- .dmsa_align_order(colnames(M), al$probe)
    if (is.numeric(.ord)) al <- al[.ord, , drop = FALSE]
  }
  if (!is.null(winsor)) M <- apply(M, 2, function(y) {
    md <- stats::median(y, na.rm = TRUE); s <- stats::mad(y, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(y)
    pmin(pmax(y, md - winsor * s), md + winsor * s) })
  Y <- scale(M); ok <- apply(Y, 2, function(z) all(is.finite(z)))
  Y <- Y[, ok, drop = FALSE]; units <- units[ok]; al <- al[ok, , drop = FALSE]
  mlt <- 2 * as.numeric(al$p_s_plus) - 1
  mlt[!is.finite(mlt) | !al$usable] <- 0
  ## reliability (centrality) weights - outcome-free, so the permutation null is
  ## exact and weighting changes power only. Needed by the "reliability" and
  ## "combined" engines; the flat engine uses all-ones.
  relw <- dmsa_relweights(Y, units, mlt, w_floor = w_floor,
    weighting = if (weighting == "flat") "flat" else "reliability")
  X <- stats::model.matrix(stats::reformulate(rhs), data)
  for (j in setdiff(colnames(X), "(Intercept)")) {
    s <- stats::sd(X[, j]); if (is.finite(s) && s > 0)
      X[, j] <- (X[, j] - mean(X[, j])) / s }

  ## CHIP AS A RANDOM INTERCEPT (Alpha contract: `(1 | chip_T1)`).
  ## The GLS solution for a one-way random intercept is OLS on quasi-demeaned
  ## data, so one REML fit gives a linear transform applied to the response and
  ## every design column BEFORE any permutation. Validity does not depend on
  ## getting gamma exactly right: the permutation null is computed under the
  ## same transform, so gamma affects power, not type-I error. gamma is
  ## estimated once per tested family, on that family's aligned score.
  ri_gamma <- NA_real_
  if (!is.null(ri_group)) {
    grp <- as.factor(ri_group)
    if (nlevels(grp) > 1L && nlevels(grp) < nrow(X)) {
      ncg <- as.numeric(table(grp)); names(ncg) <- levels(grp)
      gig <- as.character(grp)
      s0 <- as.numeric(scale(rowMeans(sweep(Y, 2, mlt, "*"))))
      if (all(is.finite(s0))) {
        ri_gamma <- .ri_reml(s0, X, grp)
        if (is.finite(ri_gamma) && ri_gamma > 0) {
          X <- .ri_transform(X, grp, ncg, gig, ri_gamma)
          Y <- .ri_transform(Y, grp, ncg, gig, ri_gamma)
        }
      }
    }
  }
  fi <- which(colnames(X) == term)
  if (!length(fi)) stop("term '", term, "' is not a design column", call. = FALSE)
  n <- nrow(X); XtXi <- solve(crossprod(X)); H <- XtXi %*% t(X)
  dfr <- n - ncol(X); vff <- XtXi[fi, fi]
  Zo <- X[, -fi, drop = FALSE]; PZ <- Zo %*% solve(crossprod(Zo), t(Zo))
  Fit <- PZ %*% Y; Res <- Y - PZ %*% Y
  xr <- as.numeric(X[, fi] - PZ %*% X[, fi]); xrss <- sum(xr^2)
  blk <- if (is.null(block)) seq_len(n) else block
  rws <- split(seq_len(n), blk); strata <- split(seq_along(rws), lengths(rws))
  PM <- matrix(0L, B, n)
  for (b in seq_len(B)) { ix <- integer(n)
    for (st in strata) ix[unlist(rws[st], use.names = FALSE)] <-
      unlist(rws[st[sample.int(length(st))]], use.names = FALSE)
    PM[b, ] <- ix }
  bse <- function(Ym) { bh <- H %*% Ym; r <- Ym - X %*% bh
    list(b = bh[fi, ], se = sqrt(colSums(r^2) / dfr * vff)) }

  ug <- unique(units); gi <- split(seq_along(units), units)[ug]
  ## Two engines share ONE permutation stream. The per-probe regression
  ## quantities (b, se; cross-products) are engine-INDEPENDENT and computed
  ## once per permutation; only the aligned WEIGHTS differ. relw_e is the
  ## per-probe reliability weight of engine e (all-ones for flat); mltw_e is
  ## the aligned multiplier under that engine. The "combined" engine runs
  ## both and fuses them per lens on the same null (joint-null maxT).
  engines <- switch(weighting,
    flat = list(flat = rep(1, length(mlt))),
    reliability = list(reliability = relw),
    combined = list(flat = rep(1, length(mlt)), reliability = relw))
  ## lens builders, parameterised by the engine's weights ------------------
  lensA_w <- function(L, relw_e) vapply(gi, function(j) {
    a <- al[j, , drop = FALSE]
    if (!any(a$usable & a$s != 0, na.rm = TRUE)) return(NA_real_)
    z <- tryCatch(dmsa_test(L$b[j], L$se[j], a, method = method,
                            w = relw_e[j])$z, error = function(e) NA_real_)
    if (length(z)) abs(z) else NA_real_ }, numeric(1))
  lensC_w <- function(Sc, mltw_e) vapply(gi, function(j) {
    w <- abs(mltw_e[j]); if (!sum(w)) return(NA_real_)
    sum(w * Sc[j]^2) / (sum(w) * xrss) }, numeric(1))
  Sc_of <- function(mltw_e) lapply(gi, function(j) {
    w <- mltw_e[j]; if (all(w == 0)) return(NULL)
    as.numeric(scale(as.numeric(Y[, j, drop = FALSE] %*% w) / sum(abs(w)))) })
  lensB_w <- function(yv, Sc_list_e) vapply(seq_along(gi), function(k) {
    s <- Sc_list_e[[k]]; if (is.null(s)) return(NA_real_)
    Xs <- cbind(s, Zo); A <- tryCatch(solve(crossprod(Xs)),
                                      error = function(e) NULL)
    if (is.null(A)) return(NA_real_)
    bh <- A %*% crossprod(Xs, yv); r <- yv - Xs %*% bh
    abs(bh[1] / sqrt(sum(r^2) / (n - ncol(Xs)) * A[1, 1])) }, numeric(1))
  yv <- as.numeric(X[, fi]); yfit <- PZ %*% yv; yres <- yv - yfit
  Sc_lists <- lapply(engines, function(re) Sc_of(mlt * re))
  ## observed + null lens statistics, per engine, on the shared stream -----
  L0 <- bse(Y); S0 <- as.numeric(crossprod(Res, xr))
  E <- lapply(names(engines), function(e) {
    re <- engines[[e]]; mw <- mlt * re
    list(A_obs = lensA_w(L0, re), C_obs = lensC_w(S0, mw),
         B_obs = lensB_w(yv, Sc_lists[[e]]),
         A_n = matrix(NA_real_, B, length(ug)),
         C_n = matrix(NA_real_, B, length(ug)),
         B_n = matrix(NA_real_, B, length(ug))) })
  names(E) <- names(engines)
  for (b in seq_len(B)) {
    ix <- PM[b, ]
    Lp <- bse(Fit + Res[ix, , drop = FALSE])
    Scp <- as.numeric(crossprod(Res[ix, , drop = FALSE], xr))
    ypv <- yfit + yres[ix]
    for (e in names(engines)) { re <- engines[[e]]; mw <- mlt * re
      E[[e]]$A_n[b, ] <- lensA_w(Lp, re)
      E[[e]]$C_n[b, ] <- lensC_w(Scp, mw)
      E[[e]]$B_n[b, ] <- lensB_w(ypv, Sc_lists[[e]]) }
  }
  ## fuse engines per lens. Single engine -> its own RAW statistic (see the
  ## width-advantage note inside fuse()); combined -> the MEAN of the two
  ## engines' null-standardised statistics - a valid statistic with its own
  ## permutation null, so maxT over units is an exact joint-null FWER.
  ## (The manuscript's equation 5 second line describes an earlier rank-p
  ## ACAT fusion; it was replaced because rank granularity 1/(B+1) floors
  ## the family-wise p at ~family/(B+1) - see the note below. MS text to be
  ## updated by the PI.)
  ## ---- continuous common scale ------------------------------------------
  ## Each lens/engine statistic is standardised by ITS OWN permutation null
  ## (robust centre/scale). This is continuous, so - unlike a rank-p scale,
  ## whose granularity is 1/(B+1) - it cannot produce ties at the extreme.
  ## A rank-p fusion floors the family-wise p at about (family size)/(B+1),
  ## which silently zeroes power for large families at modest B (found on a
  ## K = 42 grid, 15 Aug 2026). Standardising removes the floor entirely.
  zstd <- function(obs, nul) {
    ctr <- apply(nul, 2, stats::median, na.rm = TRUE)
    scl <- apply(nul, 2, stats::mad, na.rm = TRUE)
    bad <- !is.finite(scl) | scl <= 0
    if (any(bad)) { sdv <- apply(nul[, bad, drop = FALSE], 2, stats::sd, na.rm = TRUE)
      sdv[!is.finite(sdv) | sdv <= 0] <- 1; scl[bad] <- sdv }
    list(obs = (obs - ctr) / scl,
         nul = sweep(sweep(nul, 2, ctr, "-"), 2, scl, "/")) }
  fuse <- function(lens) {
    ob <- paste0(lens, "_obs"); nu <- paste0(lens, "_n")
    ## ONE ENGINE: use the raw statistic. Standardising each unit by its own
    ## permutation null equalises the family's spreads, and that is precisely
    ## what maxT's information advantage is made of - a unit whose probes all
    ## agree generates the WIDEST null in its family and is therefore charged
    ## the least. Divide that width out and maxT collapses to minP-like
    ## behaviour. It cost the pre-registered anchor its headline: AVP's
    ## family-adjusted coherence went from .0120 (raw statistic, as reported)
    ## to .1135 (standardised), which is the documented COMBINED value.
    ## Standardisation is needed only to put two engines on one scale, so it
    ## now happens only when there are two.
    if (length(engines) == 1L) {
      e <- names(engines)[1L]
      return(list(obs = E[[e]][[ob]], nul = E[[e]][[nu]]))
    }
    zs <- lapply(names(engines), function(e) zstd(E[[e]][[ob]], E[[e]][[nu]]))
    ## engine fusion: mean of the standardised statistics (equal weight, and
    ## on a scale where the two engines are commensurable)
    list(obs = Reduce(`+`, lapply(zs, `[[`, "obs")) / length(zs),
         nul = Reduce(`+`, lapply(zs, `[[`, "nul")) / length(zs)) }
  FA <- fuse("A"); FB <- fuse("B"); FC <- fuse("C")
  A_obs <- FA$obs; A_n <- FA$nul; B_obs <- FB$obs; B_n <- FB$nul
  C_obs <- FC$obs; C_n <- FC$nul
  L <- L0                                  # for direction/concordance below
  padj <- function(obs, nul) {
    p <- (1 + colSums(sweep(nul, 2, obs, ">="), na.rm = TRUE)) /
      (colSums(is.finite(nul)) + 1)
    ## `bst` is this lens's SECOND-LEVEL ingredient (PI-approved 2026-08-29):
    ## for each permutation draw, the family-best adjusted p that draw would
    ## have received - i.e. the same min-over-units transform the observed
    ## data gets. (For both corrections the family minimum of the step-down
    ## adjusted p's equals the best unit's single-step value, so this B-vector
    ## is exact for the family-min statistic.) The min of `bst` across the
    ## three lenses is one draw of the null "best adjusted p anywhere in the
    ## family, under any lens" - the calibration distribution of the
    ## second-level Westfall-Young minP (Westfall & Young 1993; Ge, Dudoit &
    ## Speed 2003).
    bst <- rep(NA_real_, nrow(nul))
    if (correction == "maxT") {
      mx <- apply(nul, 1, max, na.rm = TRUE)
      mx[!is.finite(mx)] <- NA_real_
      pa <- vapply(obs, function(o) (1 + sum(mx >= o, na.rm = TRUE)) /
                     (sum(is.finite(mx)) + 1), numeric(1))
      smx <- sort(mx[is.finite(mx)]); nf <- length(smx)
      if (nf) bst <- (1 + nf -
                        findInterval(mx, smx, left.open = TRUE)) / (nf + 1)
    } else {
      ## E10 (PI-approved): observed p above uses the (1+count)/(B+1)
      ## convention; the null p's compared against it must use the SAME
      ## convention, with the denominator counting only finite draws (frank
      ## used to divide by the full B even when a column carried NAs, which
      ## understated null p's and made p_adj conservative for the family)
      PNl <- apply(nul, 2, function(v) {
        nf <- sum(is.finite(v))
        (data.table::frank(-v, ties.method = "max", na.last = "keep") + 1) /
          (nf + 1)
      })
      mn <- apply(PNl, 1, min, na.rm = TRUE)
      mn[!is.finite(mn)] <- NA_real_
      pa <- vapply(p, function(q) (1 + sum(mn <= q, na.rm = TRUE)) /
                     (sum(is.finite(mn)) + 1), numeric(1))
      smn <- sort(mn[is.finite(mn)]); nf <- length(smn)
      if (nf) bst <- (1 + findInterval(mn, smn)) / (nf + 1)
    }
    o <- order(p); pa[o] <- cummax(pa[o])
    list(p = p, p_adj = pa, bst = bst)
  }
  PA <- padj(A_obs, A_n); PB <- padj(B_obs, B_n); PC <- padj(C_obs, C_n)
  ## ---- second-level Westfall-Young minP over the three lenses ------------
  ## (PI-approved 2026-08-29.) The any-lens naming rule takes the MIN over
  ## the three per-lens family-adjusted p's; this calibrates that minimum
  ## exactly, from the same permutation stream: Mb[b] = the smallest adjusted
  ## p that null draw b earns anywhere in the family under any lens.
  ## p_union_exact per unit = share of null draws whose family-wide best is
  ## at or below the unit's own min adjusted p - an EXACT family-wise p for
  ## the union claim, adaptive to the measured lens dependence (perfectly
  ## correlated lenses cost nothing; independent ones the full 3x). The
  ## per-lens adjusted p's and every previously reported number are
  ## untouched; this is an additional, additive quantity.
  Mb <- pmin(PA$bst, PB$bst, PC$bst, na.rm = TRUE)
  Mb[!is.finite(Mb)] <- NA_real_
  .obsmin <- pmin(PA$p_adj, PB$p_adj, PC$p_adj, na.rm = TRUE)
  .nfM <- sum(is.finite(Mb))
  p_union <- vapply(.obsmin, function(q)
    if (!is.finite(q) || !.nfM) NA_real_ else
      (1 + sum(Mb <= q, na.rm = TRUE)) / (.nfM + 1), numeric(1))
  .ou <- order(.obsmin, na.last = NA)
  if (length(.ou)) p_union[.ou] <- cummax(p_union[.ou])
  ## ---- unit-level test on ONE joint null across lenses (and engines) -----
  ## max over the three standardised lens statistics per unit, then maxT over
  ## units: a single family-wise rule valid across lenses AND units. The
  ## per-lens adjusted p-values above correct across units only, so taking
  ## their minimum is NOT family-wise controlled (it inflates roughly
  ## threefold); p_unit_adj is the rule to use for naming a unit.
  ## E3 fix (2026-08-29, PI-approved): the cross-lens max needs the three
  ## lenses on ONE scale. The combined engine's lens statistics are already
  ## null-standardised; the single-engine statistics are raw and
  ## incomparable (|z|- and |t|-like tails ~3+ against a variance-ratio near
  ## 1), so the diffuse lens could never drive the any-lens test - exactly
  ## the signal class it exists to catch. Standardising here touches ONLY
  ## the unit-level statistic: the per-lens tests above keep their raw maxT
  ## and its width advantage (the AVP .0120 anchor is untouched).
  if (length(engines) == 1L) {
    .ZA <- zstd(A_obs, A_n); .ZB <- zstd(B_obs, B_n); .ZC <- zstd(C_obs, C_n)
    Uobs <- pmax(.ZA$obs, .ZB$obs, .ZC$obs, na.rm = TRUE)
    Unul <- pmax(.ZA$nul, .ZB$nul, .ZC$nul, na.rm = TRUE)
  } else {
    Uobs <- pmax(A_obs, B_obs, C_obs, na.rm = TRUE)
    Unul <- pmax(A_n, B_n, C_n, na.rm = TRUE)
  }
  PU_ <- padj(Uobs, Unul)
  acat1 <- function(v) { v <- v[is.finite(v)]; if (!length(v)) return(NA_real_)
    q <- pmin(pmax(v, 1e-15), 1 - 1e-15)
    0.5 - atan(mean(tan((0.5 - q) * pi))) / pi }
  omni <- vapply(seq_along(ug), function(k)
    acat1(c(PA$p[k], PB$p[k], PC$p[k])), numeric(1))
  ## ---- family-corrected omnibus (PI-approved 2026-08-29) -----------------
  ## The ACAT omnibus already absorbs the cross-LENS multiplicity (Liu & Xie
  ## 2020: valid combination under arbitrary dependence); this corrects it
  ## across the FAMILY too, by permutation minP on the same stream: each
  ## null draw gets its own per-unit omnibus (ACAT of that draw's
  ## leave-one-out rank p's - the same (1+count)/(B+1) convention the
  ## observed p uses, so neither side is favoured), the family minimum per
  ## draw is the calibration distribution, and p_omnibus_adj ranks the
  ## observed omnibus against it. Exact family-wise error (measured .050 in
  ## the undermining sims, dmsa_patch/undermine_omnibus.R); powerful
  ## precisely where all lenses agree, which is what the omnibus is FOR.
  ## Additive: nothing previously reported changes value.
  .On <- {
    .pn <- function(nl) {                       # B x U -> leave-one-out p's
      apply(nl, 2, function(v) {
        nf <- sum(is.finite(v))
        (1 + nf - rank(v, ties.method = "min", na.last = "keep")) / (nf + 1)
      })
    }
    .pA <- .pn(A_n); .pB <- .pn(B_n); .pC <- .pn(C_n)
    q <- function(m) pmin(pmax(m, 1e-15), 1 - 1e-15)
    ## mean over the AVAILABLE lenses, exactly as acat1() does for the
    ## observed omnibus, so a unit with a missing lens is treated the same
    ## way on both sides
    .t1 <- tan((0.5 - q(.pA)) * pi); .t2 <- tan((0.5 - q(.pB)) * pi)
    .t3 <- tan((0.5 - q(.pC)) * pi)
    .k <- (!is.na(.t1)) + (!is.na(.t2)) + (!is.na(.t3))
    .z <- function(m) { m[is.na(m)] <- 0; m }
    Tn <- (.z(.t1) + .z(.t2) + .z(.t3)) / ifelse(.k > 0, .k, NA)
    0.5 - atan(Tn) / pi                         # B x U null omnibus p's
  }
  Mo_ <- apply(.On, 1, min, na.rm = TRUE)
  Mo_[!is.finite(Mo_)] <- NA_real_
  .nfO <- sum(is.finite(Mo_))
  omni_adj <- vapply(omni, function(x)
    if (!is.finite(x) || !.nfO) NA_real_ else
      (1 + sum(Mo_ <= x, na.rm = TRUE)) / (.nfO + 1), numeric(1))
  .oo <- order(omni, na.last = NA)
  if (length(.oo)) omni_adj[.oo] <- cummax(omni_adj[.oo])
  dir <- vapply(gi, function(j) {
    w <- mlt[j]; if (all(w == 0)) return(NA_real_)
    sign(sum(w * L$b[j] / L$se[j]^2)) }, numeric(1))
  conc <- vapply(gi, function(j) {
    w <- mlt[j]; k <- w != 0; if (!any(k)) return(NA_real_)
    t <- sign(L$b[j][k] / L$se[j][k]) * sign(w[k])
    max(mean(t > 0), mean(t < 0)) }, numeric(1))
  out <- data.frame(unit = ug,
    n_probes = vapply(gi, function(j) sum(mlt[j] != 0), integer(1)),
    concordance = conc, direction = dir,
    p_coherence = PA$p, p_coherence_adj = PA$p_adj,
    p_composite = PB$p, p_composite_adj = PB$p_adj,
    p_diffuse = PC$p, p_diffuse_adj = PC$p_adj,
    p_unit = PU_$p, p_unit_adj = PU_$p_adj,
    p_union_exact = p_union,
    p_omnibus = omni, p_omnibus_adj = omni_adj, stringsAsFactors = FALSE)
  attr(out, "union_null_min") <- Mb
  out$n_lenses_05 <- rowSums(cbind(out$p_coherence, out$p_composite,
                                   out$p_diffuse) < .05, na.rm = TRUE)
  ## units with NO usable direction-called probe are UNTESTABLE: report NA,
  ## not the p = 1 placeholder the counting formula would produce (their obs
  ## statistic is NA, so every column of the null is NA and the formula
  ## degenerates). Found in the 14 Aug all-systems run - 62% of mapped Alpha
  ## genes are in this state under the confidence map.
  unt <- out$n_probes == 0L
  if (any(unt)) {
    for (pc in c("p_coherence", "p_coherence_adj", "p_composite",
                 "p_composite_adj", "p_diffuse", "p_diffuse_adj",
                 "p_unit", "p_unit_adj", "p_union_exact", "p_omnibus",
                 "p_omnibus_adj"))
      out[[pc]][unt] <- NA_real_
    out$n_lenses_05[unt] <- NA_integer_
  }
  ## ---- the naming statistic (PI ruling 2026-10-15) -----------------------
  ## Benjamini-Hochberg on the unit p within this level-local family, with
  ## the level calibrated to the family's OWN dependence. q_unit is ordinary
  ## BH. q_unit_dep is the share of permutation draws in which BH, run at
  ## exactly the level needed to name this unit, would have named SOMETHING
  ## in a family with no signal at all - the true null false-alarm rate of
  ## that decision under the measured dependence between lenses (already
  ## priced into p_unit) and between units. Naming at q_unit_dep <= fdr
  ## therefore fires in a fraction fdr of null families whatever the
  ## correlation structure; FDR <= fdr when some units carry signal was
  ## verified by simulation on the real 12-gene oxytocin geometry
  ## (ERRCTL_DECISION_v2.md). Untestable units stay NA.
  QD <- .q_unit_dep(out$p_unit, Unul)
  out$q_unit <- QD$q
  out$q_unit_dep <- QD$q_dep
  rownames(out) <- NULL
  ## the reorder is a data.frame subset, which drops custom attributes -
  ## re-attach them after it so union_null_min reaches the caller
  out <- out[order(out$p_omnibus), ]
  attr(out, "family_size") <- length(ug)
  attr(out, "correction") <- correction
  attr(out, "weighting") <- weighting
  attr(out, "union_null_min") <- Mb
  attr(out, "omnibus_null_min") <- Mo_
  attr(out, "unit_null_qmin") <- QD$qmin
  ## spec 52: the per-probe fit the lenses were built from - b and se of
  ## each (standardised, covariate-adjusted, random-intercept-transformed
  ## when one was declared) probe on the tested term. The report's locus
  ## panel draws THESE, so the figure and the test share one estimator.
  ## Attached only; no statistic above depends on it.
  attr(out, "probe_fits") <- structure(
    data.frame(column = colnames(Y), b = as.numeric(L0$b),
               se = as.numeric(L0$se), stringsAsFactors = FALSE,
               row.names = NULL),
    dfr = dfr, ri_gamma = ri_gamma, term = term)
  out
}

## Dependence-calibrated Benjamini-Hochberg for one level-local family.
## p_obs: observed unit p's (NA = untestable). Unul: B x U permutation null of
## the unit statistic (larger = stronger), columns in the order of p_obs.
## Returns the BH q-values, the calibrated q-values and, per null draw, the
## smallest BH-adjusted p that draw would have produced (the least strict BH
## level at which it names anything).
.q_unit_dep <- function(p_obs, Unul) {
  U <- length(p_obs)
  q <- q_dep <- rep(NA_real_, U)
  ok <- is.finite(p_obs)
  if (!any(ok) || is.null(Unul) || !length(Unul))
    return(list(q = q, q_dep = q_dep, qmin = numeric(0)))
  Unul <- as.matrix(Unul)
  q[ok] <- stats::p.adjust(p_obs[ok], method = "BH")
  ## each draw's own unit p's, by the same (1 + count) / (B + 1) convention
  ## the observed p_unit uses (the draw counts itself, as the observed does)
  Pn <- apply(Unul[, ok, drop = FALSE], 2, function(v) {
    nf <- sum(is.finite(v))
    (1 + nf - rank(v, ties.method = "min", na.last = "keep")) / (nf + 1)
  })
  Pn <- matrix(Pn, ncol = sum(ok))
  ## the smallest BH-adjusted p of a draw is min_j p_(j) * k / j over its k
  ## finite p's (sorted ascending) - vectorised, no per-draw p.adjust()
  S <- if (ncol(Pn) > 1L) t(apply(Pn, 1, sort, na.last = TRUE)) else Pn
  S <- matrix(S, nrow = nrow(Pn))
  kb <- rowSums(is.finite(S))
  V <- S * outer(kb, seq_len(ncol(S)), function(k, j) k / j)
  V[!is.finite(V)] <- Inf
  qmin <- pmin(apply(V, 1, min), 1)
  qmin[kb == 0L | !is.finite(qmin)] <- NA_real_
  nf <- sum(is.finite(qmin))
  if (nf) {
    sq <- sort(qmin[is.finite(qmin)])
    q_dep[ok] <- (1 + findInterval(q[ok], sq)) / (nf + 1)
    o <- order(q[ok]); qd <- q_dep[ok]; qd[o] <- cummax(qd[o]); q_dep[ok] <- qd
  }
  list(q = q, q_dep = pmin(q_dep, 1), qmin = qmin)
}
