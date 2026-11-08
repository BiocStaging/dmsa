# ============================================================================
# dmsa_gate(): SYSTEM > MODULE > GENE > PROBE
#   one score and one model per unit
#   correction ONLY within the level and family being tested
#
# THE RULE.
#   system  the pre-registered systems are the first family. One registered
#           system means nothing to correct: alpha is spent on it.
#   module  tested only if its system passed, corrected over THAT system's
#           modules alone. Three modules in HPA is a family of three.
#   gene    tested only inside a passed module, over that module's genes alone.
#   probe   tested only inside a passed gene, over that gene's probes alone.
#           AVP's eight probes are a family of eight, not of 180 or 552.
#
# This is the family-based graphical approach to hierarchically structured
# hypothesis families (Bretz, Maurer, Brannath & Posch 2009; Qiu, Li & Hsu
# 2018): a pre-specified graph of families, alpha allocated inside a family and
# PROPAGATED to the children of a family that passes. It STRONGLY CONTROLS the
# familywise error rate over the whole tree. Correcting a gene against every
# gene on the array is not more rigorous - it is a different and much weaker
# procedure, answering a question nobody asked.
#
# TWO THINGS MAKE IT POWERFUL RATHER THAN JUST TIDY.
#
# 1. EVERY GATE IS AN OMNIBUS. A gatekeeper that only tests whether the whole
#    unit moves together will never open on a signal carried by one child - and
#    then every drop of power below it is unreachable. So each unit is tested
#    twice and combined: a DENSE arm (one model on the unit's own aligned
#    score) and a SPARSE arm (the smallest child p, calibrated against the
#    permutation distribution of that minimum). The two are combined with ACAT,
#    which is valid under arbitrary dependence and therefore costs no
#    multiplicity - unlike min(1, 2*min(p)), which was a Bonferroni over arms.
#
# 2. WITHIN-FAMILY CORRECTION IS PERMUTATION minP, NOT HOLM. Modules of one
#    system, and genes of one module, are correlated by construction. minP
#    against the joint permutation null exploits that correlation; Holm assumes
#    it away and charges the full family size. Holm remains the fallback when
#    no permutation null is available.
#
# Every unit's null is built from the SAME permutation matrix, which is what
# makes the minP families comparable.
# ============================================================================

#' Gated hierarchical DMSA: one score, one model, level-local correction
#'
#' @param data data.frame with outcome, moderators and covariates.
#' @param M Numeric matrix, samples x probes, rows aligned to \code{data}.
#' @param map data.frame, one row per column of \code{M}, with the level columns
#'   ordered outermost first, e.g. \code{system}, \code{module}, \code{gene},
#'   \code{probe}. The innermost column is the probe level.
#' @param alignment \code{dmsa_align()} result for \code{M}'s columns in order,
#'   or a numeric multiplier vector.
#' @param formula Upper-level model formula, using the literal name \code{S}
#'   where the unit's score belongs.
#' @param term Coefficient to test, as \code{model.matrix} names it - e.g.
#'   \code{"S"} or \code{"S:ACE"}.
#' @param roots Units at the outermost level that were PRE-REGISTERED. Their
#'   count is the first family's size; one root means no correction at the top.
#' @param block Permutation blocks: a column name of \code{data} (or several,
#'   combined by interaction), or per-row labels (family/cluster id).
#' @param alpha Level for every family, propagated to children that pass.
#' @param B Permutations.
#' @param winsor Passed to \code{dmsa_scores()}.
#' @param seed Optional integer.
#' @param gate Which arm(s) to run: \code{"both"} (default), \code{"sparse"}
#'   (survivors descend one unit at a time), or \code{"dense"} (a surviving
#'   family descends whole).
#' @param sparse_reach How far a sparse survivor reaches down the hierarchy:
#'   \code{"children"} (default) descends one level; a level name (e.g.
#'   \code{"probe"}) descends straight to that level.
#' @return Object of class \code{dmsa_gate}.
#' @examples
#' set.seed(1)
#' n <- 60; probes <- paste0("cg", 1:6)
#' map <- data.frame(gene = rep(c("OXTR", "AVP"), each = 3), probe = probes)
#' al <- dmsa_align(data.frame(cpg = probes, d = rep(c(1, -1), 3),
#'                             p_plus = rep(c(.9, .1), 3)),
#'                  genes = map$gene, level = "gene")
#' M <- matrix(rnorm(n * 6), n, 6, dimnames = list(NULL, probes))
#' d <- data.frame(age = rnorm(n), cID = rep(seq_len(n / 2), each = 2))
#' d$y <- as.numeric(scale(M %*% (2 * al$p_s_plus - 1))) + rnorm(n)
#' ## the two genes are the pre-registered family; a probe is corrected
#' ## inside its own gene only, never against all six
#' dmsa_gate(d, M, map, al, y ~ S + age, term = "S", roots = c("OXTR", "AVP"),
#'           block = d$cID, B = 99, seed = 1)
#' @export
dmsa_gate <- function(data, M, map, alignment, formula, term, roots,
                      block = NULL, alpha = 0.05, B = 1999, winsor = 3,
                      gate = c("both", "sparse", "dense"),
                      sparse_reach = "children", seed = NULL) {
  gate <- match.arg(gate)
  if (!is.null(seed)) {
    ## restore the caller's RNG state on exit: a permutation seed is for
    ## reproducing THIS result, not for silently reseeding the user's session.
    ## withr sets the seed and, on exit, restores the caller's state - or
    ## removes it when there was none (E10); the draws are those of set.seed().
    withr::local_seed(seed)
  }
  data <- as.data.frame(data); M <- as.matrix(M)
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  levs <- names(map)
  if (length(levs) < 2L)
    stop("map needs at least two level columns, outermost first", call. = FALSE)
  if (nrow(map) != ncol(M))
    stop("map must have one row per column of M", call. = FALSE)
  m <- if (is.numeric(alignment)) alignment else {
    al <- as.data.frame(alignment)
    mm <- 2 * as.numeric(al$p_s_plus) - 1
    us <- if (!is.null(al$usable)) as.logical(al$usable) else rep(TRUE, length(mm))
    mm[!is.finite(mm) | !us] <- 0; mm }
  if (length(m) != ncol(M)) stop("alignment length != ncol(M)", call. = FALSE)

  ## ---- one shared permutation matrix ------------------------------------
  n <- nrow(data)
  block <- .dmsa_rows(block, data, "block"); .dmsa_check_block(block)
  blk <- if (is.null(block)) seq_len(n) else block
  rws <- split(seq_len(n), blk); strata <- split(seq_along(rws), lengths(rws))
  PM <- matrix(0L, B, n)
  for (b in seq_len(B)) { ix <- integer(n)
    for (st in strata) ix[unlist(rws[st], use.names = FALSE)] <-
      unlist(rws[st[sample.int(length(st))]], use.names = FALSE)
    PM[b, ] <- ix }

  ## ---- observed and null t for one unit ---------------------------------
  unit_t <- function(cols) {
    if (!length(cols) || all(m[cols] == 0)) return(NULL)
    S <- tryCatch(dmsa_scores(M[, cols, drop = FALSE], m[cols], winsor = winsor,
                              flavours = "aligned")$aligned,
                  error = function(e) NULL)
    if (is.null(S) || !all(is.finite(S))) return(NULL)
    d2 <- data; d2$S <- S
    mf <- tryCatch(stats::model.frame(formula, d2, na.action = stats::na.omit),
                   error = function(e) NULL)
    if (is.null(mf) || nrow(mf) < 20) return(NULL)
    keep <- match(rownames(mf), rownames(d2))
    X <- stats::model.matrix(formula, mf); y <- stats::model.response(mf)
    fi <- which(colnames(X) == term); if (!length(fi)) return(NULL)
    XtXi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
    if (is.null(XtXi)) return(NULL)
    bh <- XtXi %*% crossprod(X, y); r <- y - X %*% bh
    dfr <- nrow(X) - ncol(X)
    tob <- as.numeric(bh[fi] / sqrt(sum(r^2) / dfr * XtXi[fi, fi]))
    Zo <- X[, -fi, drop = FALSE]
    PZ <- Zo %*% solve(crossprod(Zo), t(Zo))
    fit <- PZ %*% y; res <- y - fit
    tn <- numeric(B)
    for (b in seq_len(B)) {
      ix <- PM[b, keep]; ix <- match(ix, keep); ix[is.na(ix)] <- seq_along(ix)[is.na(ix)]
      yp <- fit + res[ix]
      bp <- XtXi %*% crossprod(X, yp); rp <- yp - X %*% bp
      tn[b] <- bp[fi] / sqrt(sum(rp^2) / dfr * XtXi[fi, fi])
    }
    list(t = tob, tn = tn, n = nrow(X), b = as.numeric(bh[fi]),
         r = tob / sqrt(tob^2 + dfr))
  }
  ## observed p and the B leave-one-out null p's, so minP families are exact
  to_p <- function(z) {
    a <- abs(z$t); an <- abs(z$tn)
    p_obs <- (1 + sum(an >= a, na.rm = TRUE)) / (B + 1)
    p_null <- data.table::frank(-an, ties.method = "max", na.last = "keep") / B
    list(p = p_obs, pn = p_null)
  }
  acat_p <- function(v) { v <- v[is.finite(v)]
    if (!length(v)) return(NA_real_)
    q <- pmin(pmax(v, 1e-15), 1 - 1e-15)
    0.5 - atan(mean(tan((0.5 - q) * pi))) / pi }
  acat_row <- function(Mx) apply(Mx, 1, acat_p)

  ## ---- fit every unit at every level, once ------------------------------
  cat_units <- lapply(levs, function(lv) unique(map[[lv]]))
  names(cat_units) <- levs
  P <- list(); PN <- list()
  for (lv in levs) {
    us <- cat_units[[lv]]
    P[[lv]] <- stats::setNames(rep(NA_real_, length(us)), us)
    PN[[lv]] <- matrix(NA_real_, B, length(us), dimnames = list(NULL, us))
    for (u in us) {
      z <- unit_t(which(map[[lv]] == u)); if (is.null(z)) next
      pp <- to_p(z); P[[lv]][u] <- pp$p; PN[[lv]][, u] <- pp$pn
    }
  }
  ## ---- the omnibus gate: dense arm ACAT sparse arm ----------------------
  ## sparse arm of a unit = minP over its CHILDREN, calibrated against the
  ## permutation distribution of that same minimum.
  GP <- P; GPN <- PN
  for (li in rev(seq_len(length(levs) - 1L))) {
    lv <- levs[li]
    ## WHICH DESCENDANTS THE SPARSE ARM REACHES. "children" walks the tree one
    ## step at a time, which makes the evidence pay a multiplicity toll at EVERY
    ## level - fatal for a signal that lives in one gene four levels down. Naming
    ## a level instead lets a unit's gate look straight at the depth where the
    ## signal is expected, paying that toll once. Localisation afterwards is
    ## still strictly level-local, so the reported families are unchanged.
    ch <- if (identical(sparse_reach, "children")) levs[li + 1L] else
      if (sparse_reach %in% levs[(li + 1L):length(levs)]) sparse_reach else
        levs[li + 1L]
    for (u in cat_units[[lv]]) {
      kids <- unique(map[[ch]][map[[lv]] == u])
      kids <- kids[kids %in% colnames(GPN[[ch]])]
      kids <- kids[apply(GPN[[ch]][, kids, drop = FALSE], 2,
                         function(x) any(is.finite(x)))]
      if (!length(kids)) next
      obs_min <- min(GP[[ch]][kids], na.rm = TRUE)
      null_min <- apply(GPN[[ch]][, kids, drop = FALSE], 1, min, na.rm = TRUE)
      ## the sparse arm, on the p scale, with its own permutation calibration
      rk <- data.table::frank(null_min, ties.method = "max", na.last = "keep") /
        length(null_min)
      ## COMBINE BY CALIBRATED MINIMUM, not by ACAT. A gate's only job is not
      ## to miss anything below it, and ACAT AVERAGES the Cauchy transforms, so
      ## a dead dense arm destroys a live sparse one: ACAT(.40, .055) = .106,
      ## which shut the oxytocin gate on a module sitting at .012. minP over a
      ## handful of correlated arms, calibrated against the joint permutation
      ## null of that same minimum, costs almost nothing and cannot be dragged
      ## up by an uninformative arm. ACAT stays the right choice for a final
      ## SUMMARY, where averaging evidence is what you want.
      p_sp <- (1 + sum(null_min <= obs_min, na.rm = TRUE)) /
        (sum(is.finite(null_min)) + 1)
      obs_arms <- switch(gate, dense = P[[lv]][[u]], sparse = p_sp,
                         both = pmin(P[[lv]][[u]], p_sp, na.rm = TRUE))
      null_arms <- switch(gate, dense = PN[[lv]][, u], sparse = rk,
                          both = pmin(PN[[lv]][, u], rk, na.rm = TRUE))
      GP[[lv]][u] <- (1 + sum(null_arms <= obs_arms, na.rm = TRUE)) /
        (sum(is.finite(null_arms)) + 1)
      GPN[[lv]][, u] <- data.table::frank(null_arms, ties.method = "max",
                                          na.last = "keep") / length(null_arms)
    }
  }
  ## ---- descend, correcting only within level and family -----------------
  out <- list(); passed <- roots
  for (li in seq_along(levs)) {
    lv <- levs[li]
    if (li == 1L) { units <- roots; fams <- rep("(pre-registered)", length(units))
    } else {
      keep <- map[[levs[li - 1L]]] %in% passed
      units <- unique(map[[lv]][keep])
      fams <- vapply(units, function(u)
        unique(map[[levs[li - 1L]]][map[[lv]] == u])[1], character(1))
    }
    units <- units[units %in% names(GP[[lv]])]
    if (!length(units)) { out[[lv]] <- NULL; break }
    tb <- data.frame(level = lv, unit = units, family = fams[match(units, units)],
      n_probes = vapply(units, function(u) sum(map[[lv]] == u & m != 0), integer(1)),
      p_gate = GP[[lv]][units], p_dense = P[[lv]][units],
      family_size = NA_integer_, p_adj = NA_real_, stringsAsFactors = FALSE)
    if (li > 1L) tb$family <- fams
    for (fm in unique(tb$family)) {
      k <- which(tb$family == fm & is.finite(tb$p_gate))
      tb$family_size[tb$family == fm] <- length(k)
      if (!length(k)) next
      if (length(k) == 1L) { tb$p_adj[k] <- tb$p_gate[k]; next }
      ## permutation minP within this family only
      nullmin <- apply(GPN[[lv]][, tb$unit[k], drop = FALSE], 1, min, na.rm = TRUE)
      tb$p_adj[k] <- vapply(tb$p_gate[k], function(p)
        (1 + sum(nullmin <= p, na.rm = TRUE)) / (sum(is.finite(nullmin)) + 1),
        numeric(1))
      ## step-down: an adjusted p can never be smaller than a stronger unit's
      o <- order(tb$p_gate[k]); tb$p_adj[k][o] <- cummax(tb$p_adj[k][o])
    }
    tb$passed <- is.finite(tb$p_adj) & tb$p_adj < alpha
    rownames(tb) <- NULL
    out[[lv]] <- tb
    passed <- tb$unit[tb$passed]
    if (!length(passed)) break
  }
  structure(list(levels = levs, tables = out, alpha = alpha, term = term,
                 B = B, p_all = GP, p_dense = P), class = "dmsa_gate")
}

#' @export
print.dmsa_gate <- function(x, ...) {
  cat("Gated DMSA - one score and one model per unit; correction within level ",
      "and family only\n", sep = "")
  cat("  term '", x$term, "'   alpha = ", x$alpha, "   B = ", x$B, "\n",
      "  each gate = calibrated min(dense score model, minP over children)\n\n", sep = "")
  for (lv in x$levels) {
    tb <- x$tables[[lv]]
    if (is.null(tb)) { cat("  ", lv, ": not reached\n", sep = ""); next }
    cat("  --- ", toupper(lv), " ---\n", sep = "")
    o <- order(tb$p_adj, tb$p_gate)
    own <- if (!is.null(tb$p_own)) tb$p_own else tb$p_dense
    for (i in o) cat(sprintf(
      "   %-30s fam %-24s(%2d)  own %-7s gate %-7s adj %-7s %s\n",
      substr(tb$unit[i], 1, 30), substr(tb$family[i], 1, 24), tb$family_size[i],
      ifelse(is.finite(own[i]), sprintf("%.4f", own[i]), "-"),
      ifelse(is.finite(tb$p_gate[i]),  sprintf("%.4f", tb$p_gate[i]),  "-"),
      ifelse(is.finite(tb$p_adj[i]),   sprintf("%.4f", tb$p_adj[i]),   "-"),
      ifelse(isTRUE(tb$passed[i]), "PASS", "")))
    cat("\n")
  }
  invisible(x)
}
