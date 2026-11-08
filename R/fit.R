# dmsa fitting engine
# --------------------------------------------------------------------------
# One entry point that owns the whole chain: complete-case selection under the
# design, M-transform, per-probe fitting (mixed where the design says so),
# alignment-weighted pooling, and block permutation matched to the estimator
# actually used for the null.

#' M-value transform with clamping
#' @param beta numeric matrix or vector of beta values in (0,1)
#' @param eps clamp bound
#' @return A numeric matrix or vector of M-values, the same shape as
#'   \code{beta}, on the logit scale.
#' @examples
#' beta <- matrix(runif(24, 0.05, 0.95), nrow = 6,
#'                dimnames = list(NULL, paste0("cg", 1:4)))
#' round(head(dmsa_mvalues(beta), 3), 2)
#' ## clamping is what keeps a beta of exactly 0 or 1 from becoming +/-Inf
#' dmsa_mvalues(c(0, 0.5, 1))
#' @export
dmsa_mvalues <- function(beta, eps = 1e-4) {
  b <- pmin(pmax(beta, eps), 1 - eps)
  log2(b / (1 - b))
}

## block-permutation index: blocks of equal size are swapped whole, rows keep
## their within-block order. Cross-sectionally this swaps couples; with
## repeated measures it swaps whole families across all their person-time rows.
dmsa_block_index <- function(block, B) {
  blk <- as.character(block)
  if (length(unique(blk)) == 1L && length(blk) > 1L)
    stop("every row falls in ONE permutation block - no permutation is ",
         "possible and every p-value would be 1.", call. = FALSE)
  rows <- split(seq_along(blk), blk)
  sz <- lengths(rows)
  by_size <- split(names(rows), sz)
  lapply(seq_len(B), function(b) {
    idx <- seq_along(blk)
    for (s in names(by_size)) {
      f <- by_size[[s]]
      if (length(f) < 2) next
      pm <- sample(f)
      for (k in seq_along(f)) idx[rows[[f[k]]]] <- rows[[pm[k]]]
    }
    idx
  })
}

#' Fit, pool and permute a DMSA set test under a declared design
#'
#' @param data data.frame with one row per observation.
#' @param probes Character vector of methylation columns, in the SAME order as
#'   the rows of \code{alignment}.
#' @param alignment A \code{dmsa_align()} result.
#' @param design A \code{dmsa_design()} (see \code{alpha_design()}).
#' @param method \code{"expected"} (certainty-weighted) or \code{"fixed"}.
#' @param B Permutations. \code{0} skips the permutation (descriptive only).
#' @param engine \code{"auto"} reports from \code{lmer} when the design
#'   declares random effects and \pkg{lme4} is installed, otherwise
#'   \code{lm}. The permutation null is ALWAYS built with the matched
#'   \code{lm}, and the observed \code{lm} statistic is the one compared to it
#'   - mixing estimators between observed and null is the classic way to
#'   manufacture significance. The agreement between the two is returned.
#' @param beta_input If TRUE (default) \code{probes} hold beta values and are
#'   converted to M; set FALSE if they are already M.
#' @param check Run \code{dmsa_check_design()} first (default TRUE).
#' @param seed Optional integer for reproducible permutation.
#' @return An object of class \code{dmsa_fit}.
#' @examples
#' set.seed(1)
#' n <- 80; K <- 6; x <- rnorm(n); d <- sample(c(-1, 1), K, TRUE)
#' M <- matrix(rnorm(n * K), n, K) + outer(x, 0.5 * d)
#' dat <- data.frame(cID = rep(1:40, each = 2), chip = rep(1:10, each = 8),
#'                   x = x, sex_c = rep(c(-0.5, 0.5), 40))
#' probes <- paste0("cg", 1:K); dat[probes] <- 1 / (1 + 2^(-M))
#' al <- dmsa_align(data.frame(cpg = probes, d = d,
#'                             p_plus = ifelse(d > 0, 0.95, 0.05)),
#'                  genes = rep(c("g1", "g2"), each = 3))
#' des <- dmsa_design("x", "sex_c", random = c("chip", "cID"),
#'                    exchangeable = "cID")
#' dmsa_fit(dat, probes, al, des, B = 99, engine = "lm", seed = 1)
#' @export
dmsa_fit <- function(data, probes, alignment, design,
                     method = c("expected", "fixed"),
                     B = 999L, engine = c("auto", "lm", "lmer"),
                     beta_input = TRUE, check = TRUE, seed = NULL) {
  method <- match.arg(method); engine <- match.arg(engine)
  stopifnot(inherits(design, "dmsa_design"))
  al <- as.data.frame(alignment)
  if (length(probes) != nrow(al))
    stop("probes and alignment must cover the same probes, in the same order",
         call. = FALSE)
  if (!is.null(seed)) {
    ## restore the caller's RNG state on exit: a permutation seed is for
    ## reproducing THIS result, not for silently reseeding the user's session.
    ## withr sets the seed and, on exit, restores the caller's state - or
    ## removes it when there was none (E10); the draws are those of set.seed().
    withr::local_seed(seed)
  }
  data <- as.data.frame(data)

  chk <- if (check) dmsa_check_design(design, data) else NULL

  keep <- stats::complete.cases(data[, design$vars, drop = FALSE])
  dat  <- data[keep, , drop = FALSE]
  n    <- nrow(dat)

  Y <- as.matrix(dat[, probes, drop = FALSE]); mode(Y) <- "numeric"
  if (beta_input) Y <- dmsa_mvalues(Y)
  Y <- scale(Y)                                   # per-probe z within this subsample
  ok_probe <- apply(Y, 2, function(z) any(is.finite(z)))

  fixed_part <- c(design$focal, design$fixed)
  X <- stats::model.matrix(stats::reformulate(fixed_part), dat)
  fi <- match(design$focal_test, colnames(X))
  if (is.na(fi)) {
    cand <- grep(gsub(":", ".*:", design$focal_test), colnames(X), value = TRUE)
    if (length(cand) == 1L) fi <- match(cand, colnames(X))
  }
  if (is.na(fi)) {
    ## the common cause is a categorical focal exposure: a k-level factor makes
    ## k - 1 columns, so there is no single signed effect to pool. Say so.
    fv <- design$focal_test
    hits <- grep(paste0("^", gsub("([.\\^$*+?()\\[\\]{}|])", "\\\\\\1", fv)),
                 colnames(X), value = TRUE)
    if (length(hits) > 1L)
      stop("focal term '", fv, "' expands to ", length(hits),
           " model-matrix columns (", paste(hits, collapse = ", "),
           "), so it has no single signed effect to align. DMSA is ",
           "directional: use dmsa_contrasts() to get one directional test per ",
           "contrast, or dmsa_omnibus() for a deliberately directionless test.",
           call. = FALSE)
    stop("focal term '", fv, "' is not a column of the model matrix; ",
         "got: ", paste(colnames(X), collapse = ", "), call. = FALSE)
  }

  ## ---- fast LM for every probe at once (the permutation-matched estimator)
  lm_bse <- function(Xm) {
    XtXi <- solve(crossprod(Xm))
    bh   <- XtXi %*% (t(Xm) %*% Y)
    res  <- Y - Xm %*% bh
    s2   <- colSums(res^2) / (n - ncol(Xm))
    list(b = bh[fi, ], se = sqrt(s2 * XtXi[fi, fi]))
  }
  L <- lm_bse(X)

  ## ---- reported estimator: mixed model per probe, if the design says so ----
  use_lmer <- engine == "lmer" ||
    (engine == "auto" && length(design$random_groups) > 0 &&
       requireNamespace("lme4", quietly = TRUE))
  if (engine == "lmer" && !requireNamespace("lme4", quietly = TRUE))
    stop("engine = 'lmer' needs the lme4 package", call. = FALSE)
  if (use_lmer && !length(design$random_groups))
    stop("engine = 'lmer' but the design declares no random effects", call. = FALSE)

  b_rep <- L$b; se_rep <- L$se; agreement <- NA_real_
  singular <- NA_real_; icc <- NULL
  if (use_lmer) {
    fml <- stats::as.formula(paste(".y ~", paste(c(fixed_part,
             sprintf("(1|%s)", design$random_groups)), collapse = " + ")))
    bl <- sel <- rep(NA_real_, length(probes))
    sing <- 0L; vshare <- list()
    for (j in seq_along(probes)) {
      if (!ok_probe[j]) next
      dd <- dat; dd$.y <- Y[, j]
      fit <- tryCatch(suppressMessages(suppressWarnings(
        lme4::lmer(fml, data = dd))), error = function(e) NULL)
      if (is.null(fit)) next
      if (isTRUE(lme4::isSingular(fit))) sing <- sing + 1L
      cf <- stats::coef(summary(fit))
      ## E1 fix (2026-08-29, PI-approved): a vestigial lookup here used to
      ## REASSIGN `fi` - the focal COLUMN index of the design matrix - to a
      ## ROW index of lmer's coefficient table, inside an unused assignment.
      ## lm_bse(), which builds the whole permutation null after this loop,
      ## reads `fi` from the enclosing environment: whenever lme4 dropped or
      ## reordered a term, every null statistic pooled a DIFFERENT
      ## coefficient than the observed one. `fi` now means "focal column of
      ## X" from start to finish; `k` below does the coefficient-row lookup.
      k <- match(colnames(X)[match(design$focal_test, colnames(X))], rownames(cf))
      k <- if (is.na(k)) match(design$focal_test, rownames(cf)) else k
      if (!is.na(k)) { bl[j] <- cf[k, 1]; sel[j] <- cf[k, 2] }
      vc <- as.data.frame(lme4::VarCorr(fit))
      vshare[[length(vshare) + 1L]] <- stats::setNames(vc$vcov / sum(vc$vcov), vc$grp)
    }
    if (any(!is.na(bl))) {
      agreement <- suppressWarnings(stats::cor(bl, L$b, use = "complete.obs"))
      b_rep <- bl; se_rep <- sel
    }
    singular <- sing / max(sum(ok_probe), 1L)
    if (length(vshare)) {
      nmz <- unique(unlist(lapply(vshare, names)))
      icc <- vapply(nmz, function(g)
        mean(vapply(vshare, function(v) if (g %in% names(v)) v[[g]] else NA_real_,
                    numeric(1)), na.rm = TRUE), numeric(1))
    }
  }

  observed  <- dmsa_test(b_rep, se_rep, al, method = method)
  observed_lm <- dmsa_test(L$b, L$se, al, method = method)

  ## ---- block permutation of the focal columns ------------------------------
  p_perm <- NA_real_; null_z <- numeric(0)
  if (B > 0) {
    blk <- if (!is.null(design$exchangeable)) dat[[design$exchangeable]] else seq_len(n)
    idxs <- dmsa_block_index(blk, B)
    fcols <- unique(unlist(lapply(design$focal_vars, function(v)
      grep(paste0("(^|:)", v, "($|:)"), colnames(X)))))
    fcols <- union(fcols, which(colnames(X) %in% design$focal_vars))
    src <- dat[, design$focal_vars, drop = FALSE]
    null_z <- vapply(idxs, function(i2) {
      dd <- dat
      dd[, design$focal_vars] <- src[i2, , drop = FALSE]
      Xp <- stats::model.matrix(stats::reformulate(fixed_part), dd)
      Lp <- lm_bse(Xp)
      dmsa_test(Lp$b, Lp$se, al, method = method)$z
    }, numeric(1))
    p_perm <- dmsa_perm_pvalue(observed_lm$z, null_z)
  }

  structure(list(
    call = match.call(), design = design, method = method, n = n,
    n_probes = length(probes), n_used = observed$n_used,
    estimate = observed$estimate, se = observed$se, z = observed$z,
    p_perm = p_perm, p_normal = observed$p_normal,
    engine = if (use_lmer) "lmer (reported) / lm (null)" else "lm",
    agreement_lmer_lm = agreement, singular_rate = singular, icc = icc,
    z_lm = observed_lm$z, table = observed$table, null_z = null_z,
    check = chk, B = B
  ), class = "dmsa_fit")
}

#' @export
print.dmsa_fit <- function(x, ...) {
  cat("DMSA set test (", x$method, "-sign)", if (!is.null(x$design$label))
      paste0("  |  ", x$design$label), "\n", sep = "")
  cat(sprintf("  focal %s   n = %d   probes %d used of %d\n",
              x$design$focal_test, x$n, x$n_used, x$n_probes))
  cat(sprintf("  estimate %.4f  se %.4f  z %.2f\n", x$estimate, x$se, x$z))
  if (is.finite(x$p_perm))
    cat(sprintf("  permutation p = %.4g  (B = %d, block = %s, z_null-matched = %.2f)\n",
                x$p_perm, x$B,
                if (is.null(x$design$exchangeable)) "rows" else x$design$exchangeable,
                x$z_lm))
  else cat("  permutation not run (B = 0): descriptive only\n")
  cat("  engine ", x$engine, "\n", sep = "")
  if (is.finite(x$agreement_lmer_lm))
    cat(sprintf("  lmer/lm coefficient agreement r = %.3f   singular fits %.0f%%\n",
                x$agreement_lmer_lm, 100 * x$singular_rate))
  if (length(x$icc))
    cat("  variance shares: ",
        paste(sprintf("%s %.3f", names(x$icc), x$icc), collapse = "   "), "\n", sep = "")
  invisible(x)
}
