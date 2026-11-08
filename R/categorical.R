# ============================================================================
# CATEGORICAL EXPOSURES AT THE SET LEVEL
#
# The set level regresses methylation ON the exposure (M_j ~ X + Z), so the
# exposure's type is almost never a problem: for a SINGLE focal term the t
# statistic of X in M ~ X + Z equals the t statistic of M in X ~ M + Z, so a
# binary exposure is an ordinary two-sample contrast and needs nothing special.
#
# A categorical exposure with k > 2 levels is the one real exception, and the
# reason is DMSA's whole point rather than an implementation limit. DMSA pools a
# SIGNED effect per probe. A k-level factor contributes k - 1 coefficients per
# probe, so there is no single b_j to align and no single direction to report.
# There are exactly two honest ways out:
#
#   contrasts (default, recommended) - run one directional DMSA per contrast
#     against a reference level, then correct across the k - 1 contrasts. Every
#     test keeps a sign and an interpretation ("relative to secure, the anxious
#     group has higher activation tone in the HPA set").
#
#   omnibus - pool a per-probe F on k - 1 df. Answers "does the factor move this
#     set at all" and forfeits the direction, which puts DMSA back among the
#     directionless methods it is meant to beat. Offered, never the default.
#
# Nothing here changes the subject level, where the psychological variable is
# the outcome; that is dmsa_outcome()'s job.
# ============================================================================

#' Expand a categorical focal exposure into one directional design per contrast
#'
#' DMSA pools signed per-probe effects, so a focal factor with \code{k > 2}
#' levels cannot be tested as one directional term. This builds the \code{k - 1}
#' single-contrast designs, each of which \emph{is} directional, and states the
#' multiplicity correction that applies across them.
#'
#' @param data The analysis data frame.
#' @param design A \code{dmsa_design} whose \code{focal_test} names a factor (or
#'   character, or a numeric column with few unique values) with \code{k > 2}
#'   levels.
#' @param reference Level to contrast against. Defaults to the first level,
#'   which for a factor is R's own reference.
#' @param correct Multiplicity correction across the contrasts: \code{"holm"}
#'   (default), \code{"BH"}, \code{"bonferroni"}, or \code{"none"}. Recorded so
#'   the choice is visible rather than assumed.
#' @return An object of class \code{dmsa_contrasts}: a list of
#'   \code{dmsa_design} objects in \code{$designs}, the data frame with the
#'   contrast columns added in \code{$data}, and the correction in
#'   \code{$correct}. Feed each design to \code{dmsa_fit()} and pass the
#'   resulting p-values to \code{dmsa_contrast_adjust()}.
#' @examples
#' set.seed(1)
#' d <- data.frame(style = factor(sample(c("secure", "anxious", "avoidant"), 40,
#'                                       TRUE)),
#'                 age = rnorm(40), cID = rep(1:20, each = 2))
#' des <- dmsa_design(focal = "style", fixed = "age", exchangeable = "cID")
#' dmsa_contrasts(d, des)
#' @export
dmsa_contrasts <- function(data, design, reference = NULL,
                           correct = c("holm", "BH", "bonferroni", "none")) {
  correct <- match.arg(correct)
  if (!inherits(design, "dmsa_design"))
    stop("design must come from dmsa_design() or alpha_design()", call. = FALSE)
  v <- design$focal_test
  if (grepl(":", v, fixed = TRUE))
    stop("dmsa_contrasts() expands a categorical MAIN effect. The focal term ",
         "'", v, "' is an interaction; build the contrast columns first, then ",
         "declare the interaction on one of them.", call. = FALSE)
  if (!v %in% names(data))
    stop("focal term '", v, "' is not a column of data", call. = FALSE)

  x <- data[[v]]
  if (!is.factor(x)) {
    u <- sort(unique(x[!is.na(x)]))
    if (length(u) > 12L)
      stop("'", v, "' has ", length(u), " distinct values, which looks ",
           "continuous rather than categorical. Use it directly as the focal ",
           "term.", call. = FALSE)
    x <- factor(x, levels = u)
  }
  lv <- levels(droplevels(x))
  if (length(lv) < 3L)
    stop("'", v, "' has ", length(lv), " level(s). A 2-level exposure is ",
         "already a single directional contrast - pass it to dmsa_fit() as is.",
         call. = FALSE)

  ref <- if (is.null(reference)) lv[1] else as.character(reference)
  if (!ref %in% lv)
    stop("reference level '", ref, "' is not among the levels of '", v, "': ",
         paste(lv, collapse = ", "), call. = FALSE)
  others <- setdiff(lv, ref)

  n_ref <- sum(x == ref, na.rm = TRUE)
  small <- c(ref, others)[c(n_ref, vapply(others, function(l)
    sum(x == l, na.rm = TRUE), numeric(1))) < 20L]
  if (length(small))
    warning("level(s) with fewer than 20 observations: ",
            paste(small, collapse = ", "),
            ". A contrast on a thin cell is unstable and its permutation ",
            "null is coarse.", call. = FALSE)

  dd <- data
  designs <- vector("list", length(others))
  names(designs) <- others
  for (i in seq_along(others)) {
    l <- others[i]
    cn <- .dmsa_cname(v, l, ref, names(dd))
    ## a two-level factor, not 0/1 numeric, so the contrast is only estimated on
    ## the two levels involved and the other groups drop out rather than being
    ## silently folded into the reference.
    keep <- as.character(x)
    keep[!keep %in% c(ref, l)] <- NA
    dd[[cn]] <- factor(keep, levels = c(ref, l))
    d2 <- design
    d2$focal <- c(cn, setdiff(design$focal, v))
    d2$focal_test <- paste0(cn, l)
    d2$focal_vars <- unique(c(cn, setdiff(design$focal_vars, v)))
    d2$vars <- unique(c(d2$focal_vars, design$fixed, design$random_groups,
                        design$exchangeable))
    d2$label <- paste0(if (!is.null(design$label)) paste0(design$label, ": "),
                       l, " vs ", ref)
    designs[[i]] <- d2
  }

  structure(list(designs = designs, data = dd, variable = v, reference = ref,
                 levels = lv, contrasts = others, correct = correct,
                 n_per_level = table(x)),
            class = "dmsa_contrasts")
}

.dmsa_cname <- function(v, l, ref, taken) {
  base <- paste0(v, "__", make.names(l), "_vs_", make.names(ref))
  cn <- base; k <- 1L
  while (cn %in% taken) { k <- k + 1L; cn <- paste0(base, "_", k) }
  cn
}

#' @export
print.dmsa_contrasts <- function(x, ...) {
  cat("DMSA categorical focal exposure: ", x$variable, "\n", sep = "")
  cat("  levels:    ", paste(x$levels, collapse = ", "), "\n", sep = "")
  cat("  reference: ", x$reference, "   (n = ", x$n_per_level[[x$reference]],
      ")\n", sep = "")
  cat("  ", length(x$contrasts), " directional contrast(s), corrected by ",
      x$correct, ":\n", sep = "")
  for (i in seq_along(x$designs))
    cat("    ", i, ". ", x$designs[[i]]$focal_test, "   (n = ",
        x$n_per_level[[x$contrasts[i]]], " vs ",
        x$n_per_level[[x$reference]], ")\n", sep = "")
  cat("  Each contrast keeps a sign. A single omnibus test over ",
      length(x$levels) - 1L, " df would not.\n", sep = "")
  invisible(x)
}

#' Correct p-values across the contrasts of one categorical exposure
#'
#' @param x A \code{dmsa_contrasts} object.
#' @param p Named or ordered numeric vector of p-values, one per contrast, in
#'   the order of \code{x$contrasts}.
#' @param direction Optional character or numeric vector of the pooled
#'   directions, carried through so the table reports sign next to significance.
#' @return A data frame with raw and adjusted p-values.
#' @examples
#' set.seed(1)
#' d <- data.frame(style = factor(sample(c("secure", "anxious", "avoidant"), 90,
#'                                       TRUE)),
#'                 age = rnorm(90), cID = rep(1:45, each = 2))
#' des <- dmsa_design(focal = "style", fixed = "age", exchangeable = "cID")
#' cc <- dmsa_contrasts(d, des, reference = "secure")
#' # one p-value per directional contrast, corrected within the family
#' dmsa_contrast_adjust(cc, c(anxious = 0.01, avoidant = 0.04),
#'                      direction = c("+", "-"))
#' @export
dmsa_contrast_adjust <- function(x, p, direction = NULL) {
  if (!inherits(x, "dmsa_contrasts"))
    stop("x must come from dmsa_contrasts()", call. = FALSE)
  if (length(p) != length(x$contrasts))
    stop("need one p-value per contrast: ", length(x$contrasts), " expected, ",
         length(p), " given", call. = FALSE)
  if (!is.null(names(p)) && all(x$contrasts %in% names(p)))
    p <- p[x$contrasts]
  meth <- if (x$correct == "none") "none" else x$correct
  out <- data.frame(contrast = paste0(x$contrasts, " vs ", x$reference),
                    p = as.numeric(p),
                    p_adj = if (meth == "none") as.numeric(p) else
                      stats::p.adjust(as.numeric(p), method = meth),
                    row.names = NULL)
  if (!is.null(direction)) out$direction <- direction[seq_len(nrow(out))]
  attr(out, "correct") <- meth
  out
}

#' Directionless omnibus test of a categorical focal exposure over a set
#'
#' Pools the per-probe F statistic for the \code{k - 1} contrast columns and
#' calibrates it by the same block permutation the directional tests use. This
#' answers "does this factor move the set at all" and \strong{gives up the
#' direction}, which is the property that separates DMSA from the methods it is
#' compared against. Use it as a screen, or when the levels have no ordering and
#' no reference level is defensible - and report \code{dmsa_contrasts()}
#' alongside it.
#'
#' @param data Analysis data frame.
#' @param Y Probe matrix (subjects x probes), betas or M-values.
#' @param design A \code{dmsa_design} whose \code{focal_test} names the factor.
#' @param beta_input Passed through: are the columns of \code{Y} betas?
#' @param B Permutations.
#' @param seed RNG seed.
#' @return A list with the pooled statistic, the permutation p-value, and an
#'   explicit \code{directional = FALSE}.
#' @examples
#' set.seed(1)
#' n <- 90
#' d <- data.frame(style = factor(sample(c("secure", "anxious", "avoidant"),
#'                                       n, TRUE)),
#'                 age = rnorm(n), cID = rep(seq_len(n / 2), each = 2))
#' Y <- matrix(plogis(rnorm(n * 8)), n, 8)     # beta values for 8 probes
#' des <- dmsa_design(focal = "style", fixed = "age", exchangeable = "cID")
#' ## three attachment styles with no defensible reference level: the pooled
#' ## F asks only whether the factor moves the set, and gives up the direction
#' o <- dmsa_omnibus(d, Y, des, B = 99, seed = 5)
#' round(o$p_perm, 3)
#' o$directional    # always FALSE: a pooled F on 2 df has no sign
#' @export
dmsa_omnibus <- function(data, Y, design, beta_input = TRUE, B = 999,
                         seed = 1L) {
  if (!inherits(design, "dmsa_design"))
    stop("design must come from dmsa_design() or alpha_design()", call. = FALSE)
  v <- design$focal_test
  if (!v %in% names(data)) stop("'", v, "' is not a column of data",
                                call. = FALSE)
  x <- data[[v]]
  if (!is.factor(x)) x <- factor(x)
  if (nlevels(droplevels(x)) < 3L)
    stop("an omnibus over fewer than 3 levels is just the directional test; ",
         "use dmsa_fit()", call. = FALSE)

  vars <- unique(c(v, design$fixed, design$exchangeable))
  vars <- vars[!is.na(vars)]
  cc <- stats::complete.cases(data[, vars, drop = FALSE]) &
    stats::complete.cases(Y)
  dat <- data[cc, , drop = FALSE]; Ym <- as.matrix(Y)[cc, , drop = FALSE]
  dat[[v]] <- droplevels(factor(dat[[v]]))
  if (beta_input) Ym <- dmsa_mvalues(Ym)
  Ym <- scale(Ym)
  Ym <- Ym[, apply(Ym, 2, function(z) all(is.finite(z))), drop = FALSE]
  n <- nrow(Ym); df1 <- nlevels(dat[[v]]) - 1L

  Z <- stats::model.matrix(stats::reformulate(c("1", design$fixed)), dat)
  full_rhs <- c(v, design$fixed)
  Xf <- stats::model.matrix(stats::reformulate(full_rhs), dat)
  df2 <- n - ncol(Xf)

  sumF <- function(Xfull) {
    r0 <- Ym - Z %*% solve(crossprod(Z), t(Z) %*% Ym)
    r1 <- Ym - Xfull %*% solve(crossprod(Xfull), t(Xfull) %*% Ym)
    ss0 <- colSums(r0^2); ss1 <- colSums(r1^2)
    sum(((ss0 - ss1) / df1) / (ss1 / df2))
  }
  obs <- sumF(Xf)

  blk <- if (!is.null(design$exchangeable)) dat[[design$exchangeable]] else
    seq_len(n)
  ## restore the caller's RNG state on exit: a permutation seed is for
  ## reproducing THIS result, not for silently reseeding the user's session.
  ## withr sets the seed and, on exit, restores the caller's state - or
  ## removes it when there was none (E10); the draws are those of set.seed().
  if (is.null(seed)) withr::local_preserve_seed() else withr::local_seed(seed)
  ## E6 fix (2026-08-29, PI-approved): the hand-rolled shuffle here wrote
  ## block-sorted values back into ORIGINAL row order, scrambling values
  ## across families whenever rows were not block-contiguous - a null that
  ## was neither block-respecting nor cleanly unrestricted. The shared
  ## engine routine (size-stratified whole-block swaps, rows kept in place)
  ## is used instead, one permutation scheme package-wide.
  idxs <- dmsa_block_index(blk, B)
  null <- vapply(idxs, function(o) {
    d2 <- dat; d2[[v]] <- dat[[v]][o]
    sumF(stats::model.matrix(stats::reformulate(full_rhs), d2))
  }, numeric(1))

  list(variable = v, statistic = obs, df1 = df1, n = n,
       n_probes = ncol(Ym), B = B,
       p_perm = (1 + sum(null >= obs)) / (B + 1),
       directional = FALSE,
       note = paste0("Directionless: pooled F on ", df1, " df has no sign. ",
                     "Report dmsa_contrasts() for the direction."))
}
