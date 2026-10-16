# ============================================================================
# LONGITUDINAL mDMSA - two-wave change on ONE measurement scale.
#
# With two waves, time x S x E can be represented as a change-score model only
# when S_before and S_after are defined on the SAME measurement scale. The
# implementation therefore estimates probe centring/scaling and score scaling
# once from the stacked two-wave methylation matrix and applies those locked
# transforms to both waves before Delta-S = S_after - S_before is formed.
#
# Baseline adjustment is NOT a generic "regression-to-the-mean guard". It
# changes the estimand to baseline-adjusted/residualised change and is therefore
# opt-in. Likewise, change in estimated cell composition can be a nuisance or
# part of the biological response; supplying Delta-cell covariates estimates a
# composition-adjusted rather than total bulk-tissue change.
#
# Fixed within-pair roles (mother/father, treated/control matched pairs) need a
# paired analysis. Whole-block residual permutation does not randomise a role
# whose composition is identical in every block. paired = TRUE therefore
# collapses each pair to a focal-minus-reference change and uses a studentised
# Rademacher wild bootstrap across independent pairs, with one shared bootstrap
# stream across units for maxT correction.
# ============================================================================

## Winsorise and standardise the two waves on ONE joint scale. This is the
## critical longitudinal invariant: a common shift between waves must remain a
## shift rather than disappearing because each wave was re-standardised.
.dmsa_change_stack <- function(M_before, M_after, winsor = 3) {
  M_before <- as.matrix(M_before); M_after <- as.matrix(M_after)
  n <- nrow(M_before)
  M <- rbind(M_before, M_after)
  mode(M) <- "numeric"
  if (!is.null(winsor)) M <- apply(M, 2, function(y) {
    md <- stats::median(y, na.rm = TRUE); s <- stats::mad(y, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(y)
    pmin(pmax(y, md - winsor * s), md + winsor * s)
  })
  Z <- scale(M)
  if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1L)
  ok <- apply(Z, 2, function(z) all(is.finite(z)))
  list(Z = Z, ok = ok, n = n)
}

## Project a stacked, jointly-standardised methylation matrix onto one aligned
## score, then standardise THAT score once across both waves. The two returned
## scores therefore share an origin and unit.
.dmsa_change_project <- function(st, multiplier) {
  m <- as.numeric(multiplier)
  m[!is.finite(m)] <- 0
  m[!st$ok] <- 0
  if (!any(m != 0))
    stop("no usable aligned probe after joint scaling", call. = FALSE)
  s <- as.numeric(st$Z %*% m) / sum(abs(m))
  s <- as.numeric(scale(s))
  n <- st$n
  list(S0 = s[seq_len(n)], S1 = s[n + seq_len(n)])
}

## Does every block contain the same multiset of exposure levels? If yes and
## exposure varies within every block, a whole-block permutation cannot create
## the null distribution of that fixed within-block role contrast.
.dmsa_fixed_within_block_role <- function(x, block) {
  if (is.null(block)) return(FALSE)
  rows <- split(seq_along(block), as.character(block))
  if (length(rows) < 2L || any(lengths(rows) < 2L)) return(FALSE)
  lev <- sort(unique(as.character(x)))
  if (length(lev) != 2L) return(FALSE)
  tab <- vapply(rows, function(ii) {
    tt <- table(factor(as.character(x[ii]), levels = lev))
    paste(as.integer(tt), collapse = ":")
  }, character(1))
  all(vapply(rows, function(ii)
    length(unique(as.character(x[ii]))) == 2L, logical(1))) &&
    length(unique(tab)) == 1L
}

## Resolve and validate a two-person paired design. The focal level is explicit
## by design: guessing which of two labels means "mother" or "treated" is not a
## statistical convenience we should allow.
.dmsa_change_pairs <- function(data, exposure, block, pair_focal) {
  if (is.null(block))
    stop("paired = TRUE requires `block` to identify the pair", call. = FALSE)
  if (is.null(pair_focal) || length(pair_focal) != 1L)
    stop("paired = TRUE requires an explicit single `pair_focal` level",
         call. = FALSE)
  rows <- split(seq_len(nrow(data)), as.character(block))
  if (any(lengths(rows) != 2L))
    stop("paired = TRUE requires exactly two rows in every block; block sizes are: ",
         paste(sort(unique(lengths(rows))), collapse = ", "), call. = FALSE)
  x <- data[[exposure]]
  focal <- vapply(rows, function(ii) {
    hit <- which(as.character(x[ii]) == as.character(pair_focal))
    if (length(hit) != 1L) return(NA_integer_)
    ii[hit]
  }, integer(1))
  ref <- vapply(rows, function(ii) {
    hit <- which(as.character(x[ii]) != as.character(pair_focal))
    if (length(hit) != 1L) return(NA_integer_)
    ii[hit]
  }, integer(1))
  if (anyNA(focal) || anyNA(ref))
    stop("each pair must contain exactly one `pair_focal` row and one reference row",
         call. = FALSE)
  list(focal = focal, reference = ref, pair = names(rows))
}

## HC3 studentisation for the intercept of a fixed design. Pair-level
## covariates are centred before this function is called, so the intercept is
## the adjusted mean focal-minus-reference change at the observed mean nuisance
## profile.
.dmsa_change_hc3 <- function(y, X) {
  q <- qr(X)
  if (q$rank < ncol(X))
    stop("paired change design is rank deficient", call. = FALSE)
  b <- qr.coef(q, y)
  r <- y - as.numeric(X %*% b)
  Q <- qr.Q(q, complete = FALSE)
  h <- rowSums(Q^2)
  den <- pmax(1 - h, 1e-8)
  XtXi <- solve(crossprod(X))
  u <- (r / den)^2
  meat <- crossprod(X, X * u)
  V <- XtXi %*% meat %*% XtXi
  se <- sqrt(diag(V))
  list(b = b, se = se, t = b / se, resid = r)
}

## One pair-level unit test. The bootstrap is under H0: intercept = 0, while
## nuisance slopes remain. Rademacher multipliers are supplied by the caller so
## every gene shares the SAME null stream for maxT.
.dmsa_change_paired_fit <- function(y, covariates, W) {
  y <- as.numeric(y)
  C <- as.data.frame(covariates)
  if (ncol(C)) {
    for (v in names(C)) {
      C[[v]] <- as.numeric(C[[v]])
      if (!all(is.finite(C[[v]])))
        stop("paired covariate '", v, "' contains non-finite values",
             call. = FALSE)
      C[[v]] <- C[[v]] - mean(C[[v]])
    }
    X <- cbind(`(Intercept)` = 1, as.matrix(C))
    X0 <- as.matrix(C)
  } else {
    X <- matrix(1, nrow = length(y), ncol = 1L,
                dimnames = list(NULL, "(Intercept)"))
    X0 <- NULL
  }
  obs <- .dmsa_change_hc3(y, X)
  if (is.null(X0)) {
    fit0 <- rep(0, length(y)); e0 <- y
  } else {
    q0 <- qr(X0)
    if (q0$rank < ncol(X0))
      stop("paired reduced design is rank deficient", call. = FALSE)
    b0 <- qr.coef(q0, y)
    fit0 <- as.numeric(X0 %*% b0)
    h0 <- rowSums(qr.Q(q0, complete = FALSE)^2)
    ## HC2 leverage correction of restricted residuals before the wild draw.
    ## The refitted statistic itself is HC3-studentised.
    e0 <- (y - fit0) / sqrt(pmax(1 - h0, 1e-8))
  }
  if (!nrow(W))
    return(list(b = obs$b[1], se = obs$se[1], t = obs$t[1],
                null_t = numeric()))
  tn <- vapply(seq_len(nrow(W)), function(bb) {
    yb <- fit0 + e0 * W[bb, ]
    .dmsa_change_hc3(yb, X)$t[1]
  }, numeric(1))
  list(b = obs$b[1], se = obs$se[1], t = obs$t[1], null_t = tn)
}

#' Two-wave change-score DMSA on a common longitudinal scale
#'
#' Builds aligned scores for both waves using one stacked-wave transform, then
#' tests either outcome ~ Delta-S * exposure or Delta-S ~ exposure. For fixed
#' two-person roles such as mother/father, use \code{paired = TRUE}; the method
#' then tests focal-minus-reference change across independent pairs with an HC3
#' studentised Rademacher wild bootstrap and maxT across units.
#'
#' @param M_before,M_after numeric matrices, same probes in the same order
#'   (rows = the same people in the same order).
#' @param data data.frame with outcome, exposure, moderators and covariates.
#' @param outcome outcome column (role = "predictor").
#' @param exposure exposure/moderator column. With \code{paired = TRUE}, this
#'   is the two-level within-pair role variable.
#' @param mod2 optional second moderator. Not currently available in paired mode.
#' @param units unit label per probe column (genes, modules, or one system).
#' @param alignment alignment as in \code{\link{dmsa_triangulate}}.
#' @param covariates optional change-scale nuisance covariates. Adjusting for
#'   Delta-cell composition estimates composition-adjusted, not total bulk,
#'   change and should be a scientific choice rather than an automatic default.
#' @param baseline include the before-wave aligned score. Default FALSE because
#'   this changes the estimand to baseline-adjusted/residualised change.
#' @param role "predictor": outcome ~ dS * E (* mod2) + covariates;
#'   "outcome": dS ~ E (* mod2) + covariates.
#' @param block exchangeability block. With \code{paired = TRUE}, identifies
#'   the pair and must contain exactly two rows per block.
#' @param paired logical. Use pair-level focal-minus-reference inference instead
#'   of whole-block Freedman-Lane permutation.
#' @param pair_focal explicit focal exposure level in paired mode (e.g. 1 for
#'   mother when exposure is a 0/1 mother indicator).
#' @param winsor MADs for winsorisation on the STACKED two-wave matrix; NULL
#'   disables. The same cut points are therefore applied to both waves.
#' @param B,correction,seed permutation/bootstrap controls; correction applies
#'   across the units of the tested family.
#' @param weighting Probe weighting engine. NULL chooses "combined" when
#'   methylation change is a predictor and "flat" when methylation change is
#'   the outcome. Outcome-side reliability weighting is currently refused until
#'   dedicated null-calibration simulations support it.
#' @param w_floor lower bound for reliability weighting; ignored for flat.
#' @return data.frame of class "dmsa_change": unit, n, n_probes, b, t, p,
#'   p_adj; attributes describe the estimand and inference mode.
#' @examples
#' set.seed(5)
#' n <- 80; g <- rep(c("g1", "g2"), each = 4); P <- length(g)
#' d <- rep(c(1, -1), length.out = P); E <- rnorm(n); y <- rnorm(n)
#' M0 <- matrix(rnorm(n * P), n); M1 <- 0.6 * M0 + matrix(rnorm(n * P), n) * .8
#' M1[, g == "g1"] <- M1[, g == "g1"] + outer(0.5 * y * E, d[g == "g1"])
#' al <- dmsa_align(data.frame(cpg = paste0("cg", 1:P), d = d,
#'                             p_plus = ifelse(d > 0, 0.9, 0.1)), genes = g)
#' dat <- data.frame(y = y, E = E, cv = rnorm(n))
#' dmsa_change(M0, M1, dat, outcome = "y", exposure = "E", units = g,
#'             alignment = al, covariates = "cv", B = 99, seed = 1)
#' @export
dmsa_change <- function(M_before, M_after, data, outcome = NULL, exposure,
                        mod2 = NULL, units, alignment, covariates = NULL,
                        baseline = FALSE, role = c("predictor", "outcome"),
                        block = NULL, paired = FALSE, pair_focal = NULL,
                        winsor = 3, B = 1999,
                        correction = c("maxT", "minP"),
                        weighting = NULL, w_floor = 1.5, seed = 1) {
  role <- match.arg(role); correction <- match.arg(correction)
  if (is.null(weighting))
    weighting <- if (role == "outcome") "flat" else "combined"
  weighting <- match.arg(weighting, c("combined", "reliability", "flat"))
  if (correction == "minP" && identical(weighting, "combined"))
    stop("correction = \"minP\" is not defined for the combined engine in ",
         "dmsa_change(). Use correction = \"maxT\", or a single engine.",
         call. = FALSE)
  if (role == "outcome" && weighting != "flat")
    stop("role = 'outcome' currently requires weighting = 'flat'. Reliability ",
         "weights are estimated from the methylation outcome itself and are ",
         "not enabled for confirmatory change inference until null calibration ",
         "has been demonstrated.", call. = FALSE)
  if (paired && role != "outcome")
    stop("paired = TRUE is currently defined for role = 'outcome' only",
         call. = FALSE)
  if (paired && !is.null(mod2))
    stop("mod2 is not currently implemented in paired change mode",
         call. = FALSE)
  if (paired && weighting != "flat")
    stop("paired change mode currently requires weighting = 'flat'",
         call. = FALSE)
  if (paired && correction != "maxT")
    stop("paired change mode currently supports correction = 'maxT' only",
         call. = FALSE)

  M_before <- as.matrix(M_before); M_after <- as.matrix(M_after)
  stopifnot(ncol(M_before) == ncol(M_after),
            nrow(M_before) == nrow(M_after),
            length(units) == ncol(M_before))
  data <- as.data.frame(data)
  stopifnot(nrow(data) == nrow(M_before))
  block <- .dmsa_rows(block, data, "block"); .dmsa_check_block(block)
  if (!exposure %in% names(data))
    stop("exposure column '", exposure, "' not in data", call. = FALSE)
  if (role == "predictor" && (is.null(outcome) || !outcome %in% names(data)))
    stop("role = 'predictor' needs an `outcome` column", call. = FALSE)
  if (!is.null(mod2) && !mod2 %in% names(data))
    stop("mod2 column '", mod2, "' not in data", call. = FALSE)
  if (!is.null(covariates)) {
    miss <- setdiff(covariates, names(data))
    if (length(miss))
      stop("covariate column(s) absent: ", paste(miss, collapse = ", "),
           call. = FALSE)
  } else {
    message("no change-scale nuisance covariates supplied. This estimates total ",
            "bulk-tissue change; add Delta-cell/technical covariates only when ",
            "the desired estimand is composition/technical-adjusted change.")
  }
  if (baseline)
    message("baseline = TRUE estimates BASELINE-ADJUSTED (residualised) change; ",
            "this is a different estimand from unadjusted change.")

  al <- as.data.frame(alignment)
  if (nrow(al) != ncol(M_before))
    stop("alignment must have one row per methylation matrix column",
         call. = FALSE)
  if (!is.null(al$probe) && anyDuplicated(paste(al$probe, units, sep = "\r")))
    stop("the same CpG appears more than once inside the same unit. A CpG may ",
         "map to several genes, but duplicating one CpG x unit pair would ",
         "pseudo-replicate that measurement.", call. = FALSE)

  ## A fixed within-block role requires the paired estimator; whole-block
  ## residual permutation only moves complete role-pattern vectors around.
  if (!paired && role == "outcome" &&
      .dmsa_fixed_within_block_role(data[[exposure]], block))
    stop("exposure '", exposure, "' is a fixed two-level within-block role: ",
         "every block has the same role composition. Whole-block residual ",
         "permutation cannot calibrate that contrast. Re-run with paired = TRUE ",
         "and an explicit pair_focal level.", call. = FALSE)

  pair_index <- if (paired)
    .dmsa_change_pairs(data, exposure, block, pair_focal) else NULL
  if (paired && length(pair_index$pair) < 8L)
    warning("paired change inference has only ", length(pair_index$pair),
            " independent pairs; HC3/wild-bootstrap tails may be unstable")

  d0 <- data
  if (!paired) {
    d0$.E <- as.numeric(scale(d0[[exposure]]))
    if (!is.null(mod2)) d0$.M2 <- as.numeric(scale(d0[[mod2]]))
  }

  ## One shared Rademacher stream across units -> joint maxT null.
  W <- matrix(numeric(0), 0, if (paired) length(pair_index$pair) else 0)
  if (paired && B > 0) {
    old_seed <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old_seed))
      assign(".Random.seed", old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv()))
        rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(seed)
    W <- matrix(sample(c(-1, 1), B * length(pair_index$pair), replace = TRUE),
                nrow = B, ncol = length(pair_index$pair))
  }

  engs <- if (weighting == "combined") c("flat", "reliability") else weighting
  cauchy <- function(p) tan((0.5 - pmin(pmax(p, 1e-12), 1 - 1e-12)) * pi)
  ug <- unique(units)
  res <- list(); fobs <- list(); fnul <- list()

  for (u in ug) {
    j <- which(units == u)
    aj <- al[j, , drop = FALSE]
    usable <- if (!is.null(aj$usable))
      sum(aj$usable & aj$s != 0, na.rm = TRUE) else
      sum(is.finite(aj$s) & aj$s != 0)
    if (usable < 2) next
    mj <- 2 * as.numeric(aj$p_s_plus) - 1
    us <- if (!is.null(aj$usable))
      as.logical(aj$usable) else rep(TRUE, length(mj))
    mj[!is.finite(mj) | !us] <- 0

    ## ONE transform for both waves. Reliability weights, when allowed, are
    ## estimated from this same stacked matrix.
    st <- .dmsa_change_stack(M_before[, j, drop = FALSE],
                             M_after[, j, drop = FALSE], winsor = winsor)
    rw_rel <- rep(1, length(mj))
    if ("reliability" %in% engs && any(st$ok)) {
      rw_rel[st$ok] <- dmsa_relweights(
        st$Z[, st$ok, drop = FALSE], rep(1L, sum(st$ok)), mj[st$ok],
        w_floor = w_floor, weighting = "reliability")
    }

    if (paired) {
      sc <- .dmsa_change_project(st, mj)
      dS <- sc$S1 - sc$S0
      fi <- pair_index$focal; ri <- pair_index$reference
      ypair <- dS[fi] - dS[ri]
      C <- data.frame(row.names = pair_index$pair)
      if (baseline)
        C$baseline_difference <- sc$S0[fi] - sc$S0[ri]
      if (length(covariates)) {
        for (cv in covariates) {
          vv <- data[[cv]]
          if (!is.numeric(vv))
            stop("paired change covariate '", cv,
                 "' must be numeric so a focal-minus-reference difference is defined",
                 call. = FALSE)
          C[[paste0("d_", cv)]] <- vv[fi] - vv[ri]
        }
      }
      ft <- .dmsa_change_paired_fit(ypair, C, W)
      so <- abs(ft$t); sn <- abs(ft$null_t)
      p0 <- if (length(sn))
        (1 + sum(sn >= so)) / (length(sn) + 1) else NA_real_
      res[[u]] <- data.frame(unit = u, n = length(ypair), n_probes = usable,
        b = ft$b, t = ft$t, p = p0, stringsAsFactors = FALSE)
      fobs[[u]] <- so; fnul[[u]] <- sn
      next
    }

    covs <- c(if (baseline) ".S0", covariates)
    if (role == "predictor") {
      lhs <- outcome
      inter <- if (is.null(mod2)) ".dS * .E" else ".dS * .E * .M2"
      term <- if (is.null(mod2)) ".dS:.E" else ".dS:.E:.M2"
    } else {
      lhs <- ".dS"
      inter <- if (is.null(mod2)) ".E" else ".E * .M2"
      term <- if (is.null(mod2)) ".E" else ".E:.M2"
    }
    ff <- stats::as.formula(
      paste(lhs, "~", paste(c(inter, covs), collapse = " + ")))

    ## Same seed -> same block-permutation stream across units/engines, provided
    ## they use the same complete-case rows.
    fits <- list()
    for (e in engs) {
      meff <- mj * (if (e == "flat") rep(1, length(mj)) else rw_rel)
      sc <- .dmsa_change_project(st, meff)
      d2 <- d0
      d2$.dS <- as.numeric(scale(sc$S1 - sc$S0))
      d2$.S0 <- as.numeric(scale(sc$S0))
      ft <- tryCatch(
        dmsa_model(ff, d2, term, block = block, B = B, seed = seed,
                   nulls = TRUE),
        error = function(e) NULL)
      if (is.null(ft)) { fits <- NULL; break }
      fits[[e]] <- ft
    }
    if (is.null(fits) || !length(fits)) next
    ref <- fits[[1]]
    if (length(fits) == 1L) {
      so <- abs(ref$t); sn <- abs(ref$null_t)
    } else {
      rp <- lapply(fits, function(f) {
        pool <- c(abs(f$t), abs(f$null_t)); m <- length(pool)
        r <- (m - data.table::frank(pool, ties.method = "min") + 1) / m
        list(po = r[1], pn = r[-1])
      })
      so <- Reduce(`+`, lapply(rp, function(r) cauchy(r$po)))
      sn <- Reduce(`+`, lapply(rp, function(r) cauchy(r$pn)))
    }
    res[[u]] <- data.frame(unit = u, n = ref$n, n_probes = usable,
      b = ref$b, t = ref$t,
      p = (1 + sum(sn >= so)) / (length(sn) + 1),
      stringsAsFactors = FALSE)
    fobs[[u]] <- so; fnul[[u]] <- sn
  }

  if (!length(res))
    stop("no unit had >= 2 usable direction-called probes", call. = FALSE)
  out <- do.call(rbind, res); rownames(out) <- NULL
  SN <- do.call(cbind, fnul); so_v <- unlist(fobs)

  ## B=0 is descriptive only.
  if (!length(fnul[[1]])) {
    out$p_adj <- NA_real_
  } else if (!is.null(SN) && ncol(SN) > 1) {
    if (correction == "minP" && length(engs) == 1L && !paired) {
      PN <- apply(SN, 2, function(v)
        (data.table::frank(-v, ties.method = "max") + 1) /
          (sum(is.finite(v)) + 1))
      mn <- apply(PN, 1, min, na.rm = TRUE)
      out$p_adj <- vapply(out$p, function(q)
        (1 + sum(mn <= q, na.rm = TRUE)) /
          (sum(is.finite(mn)) + 1), numeric(1))
    } else {
      mx <- apply(SN, 1, max, na.rm = TRUE)
      mx[!is.finite(mx)] <- NA_real_
      out$p_adj <- vapply(so_v, function(o)
        (1 + sum(mx >= o, na.rm = TRUE)) /
          (sum(is.finite(mx)) + 1), numeric(1))
    }
    o <- order(out$p, na.last = NA)
    if (length(o)) out$p_adj[o] <- cummax(out$p_adj[o])
  } else {
    out$p_adj <- out$p
  }

  attr(out, "term") <- if (paired) {
    paste0("paired focal-minus-reference dS (focal = ", pair_focal, ")")
  } else if (role == "predictor") {
    if (is.null(mod2)) "dS x E (two-wave change on common scale)"
    else "dS x E x mod2"
  } else paste("exposure ->", "dS")
  attr(out, "role") <- role
  attr(out, "correction") <- correction
  attr(out, "paired") <- paired
  attr(out, "baseline_adjusted") <- baseline
  attr(out, "inference") <- if (paired)
    "pair-level HC3 studentised Rademacher wild bootstrap" else
    "Freedman-Lane residual permutation (approximate unless null exchangeability is exact)"
  attr(out, "scale") <-
    "probe and score transforms estimated once on stacked before+after waves"
  class(out) <- c("dmsa_change", "data.frame")
  out[order(out$p), ]
}

#' @export
print.dmsa_change <- function(x, ...) {
  cat("DMSA change-score test -", attr(x, "term"), "\n")
  cat("  scale:", attr(x, "scale"), "\n")
  cat("  inference:", attr(x, "inference"), "\n")
  cat("  baseline adjustment:",
      if (isTRUE(attr(x, "baseline_adjusted"))) "YES" else "NO", "\n")
  cat("  family of", nrow(x), "unit(s),", attr(x, "correction"),
      "across the family\n")
  for (i in seq_len(nrow(x)))
    cat(sprintf(
      "  %-28s n %4d  probes %3d  b %+.3f  t %+.2f  p %.4f  adj %.4f%s\n",
      x$unit[i], x$n[i], x$n_probes[i], x$b[i], x$t[i], x$p[i],
      x$p_adj[i],
      if (is.finite(x$p_adj[i]) && x$p_adj[i] < .05) "  *" else ""))
  invisible(x)
}
