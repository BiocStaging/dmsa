# ============================================================================
# dmsa_levels(): the powerful form of the gated cascade
#
# WHY NOT THE SCORE. Collapsing a unit's probes into one composite and running
# one regression is the right thing for a MODEL - interactions, categorical
# outcomes, mediation, an effect size to report. It is the wrong thing for a
# TEST, and measurably so. AVP's eight CpGs all carry d = -1 and all move the
# same way; that agreement IS the signal. A composite averages them, and with
# a within-gene correlation near .78 the average gains only about 0.43 of the
# sqrt(8) it would gain from independent probes. Empirically, on oxytocin ->
# attachment anxiety with the family of 42 genes:
#
#     coefficient-pooled aligned statistic   AVP maxT p = .0095
#     composite score, one model per gene    AVP maxT p = .098
#
# A ten-fold power difference, both permutation-exact, because the pooled
# statistic REWARDS CROSS-PROBE CONCORDANCE and the composite does not. For a
# coherence method that is the whole game.
#
# So: TEST with the pooled statistic, MODEL with the score. Both calibrated by
# the same permutation, so they are two readings of one analysis.
#
# THE LEVELS AND THE FAMILIES ARE EXACTLY AS SPECIFIED.
#   system  pre-registered family. One system, no correction.
#   module  only if its system passed; corrected over that system's modules.
#   gene    only inside a passed module; over that module's genes.
#   probe   only inside a passed gene; over that gene's probes.
# Within a family the correction is permutation minP, which exploits the
# correlation the family has by construction instead of assuming it away.
# Each unit's gate is the calibrated minimum of its own pooled statistic and
# the minP over the descendants it is allowed to see.
# ============================================================================

#' Gated hierarchical DMSA on the pooled aligned statistic
#'
#' @param M Numeric matrix, samples x probes, methylation on the analysis scale.
#'   The probes are the RESPONSES, as in \code{dmsa_fit()}.
#' @param data data.frame of predictors, same row order as \code{M}.
#' @param rhs Character vector of right-hand-side terms, e.g.
#'   \code{c("anx","avo","sex_c","chip_f")} or with an interaction
#'   \code{c("anx*aceZ", ...)}.
#' @param term The design column to pool, as \code{model.matrix} names it -
#'   \code{"anx"} or \code{"anx:aceZ"}.
#' @param map data.frame, one row per column of \code{M}, level columns
#'   outermost first (e.g. system, module, gene, probe).
#' @param alignment \code{dmsa_align()} result for \code{M}'s columns in order.
#' @param roots Pre-registered units at the outermost level.
#' @param block Permutation block labels, one per row.
#' @param alpha,B,winsor,seed As in \code{\link{dmsa_gate}}.
#' @param method Pooling method for \code{dmsa_test()}.
#' @param gate \code{"both"}, \code{"sparse"} or \code{"dense"}.
#' @param sparse_reach \code{"children"} or the name of a deeper level a unit's
#'   sparse arm may look straight at.
#' @param family_correction Character. Family-wise correction applied inside each level-local family: \code{"maxT"} (default) Westfall-Young step-down on the strength scale, or \code{"minP"} step-down on marginal p-values. The choice is consequential and should be fixed in advance; see the package vignette and \code{REPRODUCE.md}.
#' @return Object of class \code{dmsa_gate}.
#' @examples
#' set.seed(1)
#' n <- 60; probes <- paste0("cg", 1:6)
#' map <- data.frame(gene = rep(c("OXTR", "AVP"), each = 3), probe = probes)
#' al <- dmsa_align(data.frame(cpg = probes, d = rep(c(1, -1), 3),
#'                             p_plus = rep(c(.9, .1), 3)),
#'                  genes = map$gene, level = "gene")
#' d <- data.frame(anx = rnorm(n), age = rnorm(n),
#'                 cID = rep(seq_len(n / 2), each = 2))
#' ## the probes are the RESPONSES: their aligned coefficients are pooled,
#' ## so cross-probe concordance is what carries the evidence
#' M <- matrix(rnorm(n * 6), n, 6) + outer(0.3 * d$anx, 2 * al$p_s_plus - 1)
#' dmsa_levels(M, d, rhs = c("anx", "age"), term = "anx", map = map,
#'             alignment = al, roots = c("OXTR", "AVP"), block = d$cID,
#'             B = 99, seed = 1)
#' @export
dmsa_levels <- function(M, data, rhs, term, map, alignment, roots, block = NULL,
                      alpha = 0.05, B = 1999, winsor = 3,
                      method = c("expected", "fixed"),
                      gate = c("both", "sparse", "dense"),
                      family_correction = c("maxT", "minP"),
                      sparse_reach = "children", seed = NULL) {
  method <- match.arg(method); gate <- match.arg(gate)
  family_correction <- match.arg(family_correction)
  if (!is.null(seed)) {
    ## restore the caller's RNG state on exit: a permutation seed is for
    ## reproducing THIS result, not for silently reseeding the user's session.
    ## withr sets the seed and, on exit, restores the caller's state - or
    ## removes it when there was none (E10); the draws are those of set.seed().
    withr::local_seed(seed)
  }
  M <- as.matrix(M); data <- as.data.frame(data)
  map <- as.data.frame(map, stringsAsFactors = FALSE); levs <- names(map)
  if (nrow(map) != ncol(M)) stop("map must have one row per column of M",
                                 call. = FALSE)
  al <- as.data.frame(alignment)
  if (nrow(al) != ncol(M)) stop("alignment must cover M's columns in order",
                                call. = FALSE)
  if (!is.null(winsor)) M <- apply(M, 2, function(y) {
    md <- stats::median(y, na.rm = TRUE); s <- stats::mad(y, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(y)
    pmin(pmax(y, md - winsor * s), md + winsor * s) })
  Y <- scale(M)
  ok <- apply(Y, 2, function(z) all(is.finite(z)))
  Y <- Y[, ok, drop = FALSE]; map <- map[ok, , drop = FALSE]; al <- al[ok, ]
  X <- stats::model.matrix(stats::reformulate(rhs), data)
  for (j in setdiff(colnames(X), "(Intercept)")) {
    sdj <- stats::sd(X[, j]); if (is.finite(sdj) && sdj > 0)
      X[, j] <- (X[, j] - mean(X[, j])) / sdj }
  fi <- which(colnames(X) == term)
  if (!length(fi)) stop("term '", term, "' is not a design column. Columns: ",
                        paste(utils::head(colnames(X), 20), collapse = ", "),
                        call. = FALSE)
  n <- nrow(X); XtXi <- solve(crossprod(X)); H <- XtXi %*% t(X)
  dfr <- n - ncol(X); vff <- XtXi[fi, fi]
  bse <- function(Ym) { bh <- H %*% Ym; r <- Ym - X %*% bh
    list(b = bh[fi, ], se = sqrt(colSums(r^2) / dfr * vff)) }
  Zo <- X[, -fi, drop = FALSE]; PZ <- Zo %*% solve(crossprod(Zo), t(Zo))
  Fit <- PZ %*% Y; Res <- Y - Fit
  block <- .dmsa_rows(block, data, "block"); .dmsa_check_block(block)
  blk <- if (is.null(block)) seq_len(n) else block
  rws <- split(seq_len(n), blk); strata <- split(seq_along(rws), lengths(rws))
  PM <- matrix(0L, B, n)
  for (b in seq_len(B)) { ix <- integer(n)
    for (st in strata) ix[unlist(rws[st], use.names = FALSE)] <-
      unlist(rws[st[sample.int(length(st))]], use.names = FALSE)
    PM[b, ] <- ix }

  ## one panel fit per permutation, then every unit is a cheap pooling
  L <- bse(Y)
  LB <- matrix(NA_real_, B, ncol(Y)); LS <- LB
  for (b in seq_len(B)) { Lp <- bse(Fit + Res[PM[b, ], , drop = FALSE])
    LB[b, ] <- Lp$b; LS[b, ] <- Lp$se }
  pool <- function(cols, bb, ss) {
    a <- al[cols, , drop = FALSE]
    if (!any(a$usable & a$s != 0, na.rm = TRUE)) return(NA_real_)
    z <- tryCatch(dmsa_test(bb[cols], ss[cols], a, method = method)$z,
                  error = function(e) NA_real_)
    if (is.null(z) || !length(z)) NA_real_ else z }

  units <- lapply(levs, function(lv) unique(map[[lv]])); names(units) <- levs
  P <- list(); PN <- list(); Z <- list(); ZN <- list()
  for (lv in levs) {
    us <- units[[lv]]
    P[[lv]] <- stats::setNames(rep(NA_real_, length(us)), us)
    Z[[lv]] <- P[[lv]]
    PN[[lv]] <- matrix(NA_real_, B, length(us), dimnames = list(NULL, us))
    ZN[[lv]] <- PN[[lv]]
    for (u in us) {
      cols <- which(map[[lv]] == u)
      zo <- pool(cols, L$b, L$se); if (!is.finite(zo)) next
      zn <- vapply(seq_len(B), function(b) pool(cols, LB[b, ], LS[b, ]),
                   numeric(1))
      an <- abs(zn)
      Z[[lv]][u] <- abs(zo); ZN[[lv]][, u] <- an
      P[[lv]][u] <- (1 + sum(an >= abs(zo), na.rm = TRUE)) /
        (sum(is.finite(an)) + 1)
      PN[[lv]][, u] <- data.table::frank(-an, ties.method = "max",
                                        na.last = "keep") / B
    }
  }
  ## maxT works on the RAW pooled |z|; minP standardises each unit to its own
  ## null first. The difference is not cosmetic. A gene with eight concordant
  ## probes has a WIDER null than a gene with one probe, and maxT lets it keep
  ## that advantage - which is exactly right for a coherence method, because
  ## eight agreeing probes carry more information than one. minP deliberately
  ## discards it and then charges the full family size. On oxytocin -> anxiety
  ## the same AVP finding is .0095 under maxT and .24 under minP over the 42
  ## genes. Both control FWER; they weight the alpha budget differently, and
  ## the choice has to be pre-registered.
  if (family_correction == "maxT") { PP <- Z; PPN <- ZN
  } else { PP <- lapply(P, function(v) -v); PPN <- lapply(PN, function(m) -m) }
  ## gates, deepest level upward
  GP <- P; GPN <- PN; GS <- PP; GSN <- PPN
  for (li in rev(seq_len(length(levs) - 1L))) {
    lv <- levs[li]
    ch <- if (identical(sparse_reach, "children")) levs[li + 1L] else
      if (sparse_reach %in% levs[(li + 1L):length(levs)]) sparse_reach else
        levs[li + 1L]
    for (u in units[[lv]]) {
      kids <- unique(map[[ch]][map[[lv]] == u])
      kids <- kids[kids %in% colnames(GPN[[ch]])]
      kids <- kids[apply(GPN[[ch]][, kids, drop = FALSE], 2,
                         function(x) any(is.finite(x)))]
      if (!length(kids)) next
      om <- max(GS[[ch]][kids], na.rm = TRUE)
      nm <- apply(GSN[[ch]][, kids, drop = FALSE], 1, max, na.rm = TRUE)
      p_sp <- (1 + sum(nm >= om, na.rm = TRUE)) / (sum(is.finite(nm)) + 1)
      rk <- data.table::frank(-nm, ties.method = "max", na.last = "keep") /
        length(nm)
      oa <- switch(gate, dense = P[[lv]][[u]], sparse = p_sp,
                   both = pmin(P[[lv]][[u]], p_sp, na.rm = TRUE))
      na_ <- switch(gate, dense = PN[[lv]][, u], sparse = rk,
                    both = pmin(PN[[lv]][, u], rk, na.rm = TRUE))
      GP[[lv]][u] <- (1 + sum(na_ <= oa, na.rm = TRUE)) /
        (sum(is.finite(na_)) + 1)
      GPN[[lv]][, u] <- data.table::frank(na_, ties.method = "max",
                                          na.last = "keep") / length(na_)
      ## the gate on the strength scale, so a parent can be maxT-corrected too
      GS[[lv]][u] <- -GP[[lv]][u]; GSN[[lv]][, u] <- -GPN[[lv]][, u]
    }
  }
  ## descend
  out <- list(); passed <- roots
  for (li in seq_along(levs)) {
    lv <- levs[li]
    if (li == 1L) { us <- roots; fams <- rep("(pre-registered)", length(us))
    } else {
      keep <- map[[levs[li - 1L]]] %in% passed
      us <- unique(map[[lv]][keep])
      fams <- vapply(us, function(u)
        unique(map[[levs[li - 1L]]][map[[lv]] == u])[1], character(1)) }
    us <- us[us %in% names(GP[[lv]])]
    if (!length(us)) { out[[lv]] <- NULL; break }
    tb <- data.frame(level = lv, unit = us, family = fams[seq_along(us)],
      n_probes = vapply(us, function(u)
        sum(map[[lv]] == u & al$usable & al$s != 0, na.rm = TRUE), integer(1)),
      p_own = P[[lv]][us], p_gate = GP[[lv]][us],
      family_size = NA_integer_, p_adj = NA_real_, stringsAsFactors = FALSE)
    if (li > 1L) tb$family <- fams
    for (fm in unique(tb$family)) {
      k <- which(tb$family == fm & is.finite(tb$p_gate))
      tb$family_size[tb$family == fm] <- length(k)
      if (!length(k)) next
      if (length(k) == 1L) { tb$p_adj[k] <- tb$p_gate[k]; next }
      ## family correction: maxT on the strength scale (raw |z| at the deepest
      ## reach, or -gate p once a gate has been formed), else minP on p
      SS <- if (family_correction == "maxT" && li == length(levs)) Z[[lv]] else
        -GP[[lv]]
      SSN <- if (family_correction == "maxT" && li == length(levs)) ZN[[lv]] else
        -GPN[[lv]]
      nullmax <- apply(SSN[, tb$unit[k], drop = FALSE], 1, max, na.rm = TRUE)
      pa <- vapply(SS[tb$unit[k]], function(z)
        (1 + sum(nullmax >= z, na.rm = TRUE)) / (sum(is.finite(nullmax)) + 1),
        numeric(1))
      o <- order(tb$p_gate[k]); pa[o] <- cummax(pa[o])
      names(pa) <- NULL
      tb$p_adj[k] <- pa
    }
    tb$passed <- is.finite(tb$p_adj) & tb$p_adj < alpha
    rownames(tb) <- NULL; out[[lv]] <- tb
    passed <- tb$unit[tb$passed]
    if (!length(passed)) break
  }
  structure(list(levels = levs, tables = out, alpha = alpha, term = term,
                 B = B, p_own = P, p_gate = GP, z_own = Z,
                 family_correction = family_correction,
                 statistic = "pooled aligned"),
            class = "dmsa_gate")
}
