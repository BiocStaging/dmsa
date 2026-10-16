# dmsa_change() paired-mode null calibration
# ------------------------------------------
# Slow validation harness; not part of routine R CMD check.
# Run after installing the package version being validated.
#
# It targets the failure modes that motivated 0.99.18:
#   * fixed mother/father role within every pair;
#   * pair dependence;
#   * role-specific heteroscedasticity;
#   * change-scale nuisance covariate shifts;
#   * correlated genes / shared null stream;
#   * global-null FWER and partial-null strong-FWER behaviour.
#
# Example:
#   source(system.file("validation", "dmsa_change_null_calibration.R",
#                      package = "dmsa"))
#   ans <- dmsa_change_null_calibration(R = 500, B = 499, n_pairs = 60)
#   ans$summary

.dcc_make <- function(n_pairs = 60, planted = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- 2L * n_pairs
  genes <- rep(paste0("g", 1:5), each = 4)
  P <- length(genes)
  d <- rep(c(1, -1, 1, -1), 5)
  role <- rep(c(0, 1), n_pairs)
  cid <- rep(seq_len(n_pairs), each = 2)

  # Pair dependence plus person-specific noise.
  pair0 <- matrix(rnorm(n_pairs * P, sd = .55), n_pairs, P)
  M0 <- pair0[cid, , drop = FALSE] + matrix(rnorm(n * P, sd = .75), n, P)

  # Cell/technical change differs by role. It is a nuisance here by design.
  dcell <- rnorm(n, mean = .20 * role, sd = .20)
  pair_delta <- matrix(rnorm(n_pairs * P, sd = .25), n_pairs, P)
  eps <- matrix(rnorm(n * P), n, P)
  eps[role == 1, ] <- 1.55 * eps[role == 1, , drop = FALSE] # heteroscedastic
  M1 <- .65 * M0 + pair_delta[cid, , drop = FALSE] + .55 * dcell + .65 * eps

  if (planted) {
    j <- which(genes == "g1")
    M1[role == 1, j] <- M1[role == 1, j, drop = FALSE] +
      matrix(rep(.70 * d[j], each = n_pairs), nrow = n_pairs)
  }

  al <- dmsa_align(data.frame(cpg = paste0("cg", seq_len(P)), d = d,
                              p_plus = ifelse(d > 0, .95, .05)),
                   genes = genes, level = "gene")
  list(M0 = M0, M1 = M1, genes = genes, al = al,
       dat = data.frame(mother = role, cid = cid, dcell = dcell))
}

dmsa_change_null_calibration <- function(R = 500L, B = 499L,
                                         n_pairs = 60L, seed = 20260912L,
                                         alpha = .05) {
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, 2L * R)

  global_any <- logical(R)
  global_point <- matrix(FALSE, R, 5,
                         dimnames = list(NULL, paste0("g", 1:5)))
  partial_any_null <- logical(R)
  partial_power <- logical(R)

  for (r in seq_len(R)) {
    z <- .dcc_make(n_pairs = n_pairs, planted = FALSE, seed = seeds[r])
    fit <- suppressMessages(dmsa_change(
      z$M0, z$M1, z$dat, exposure = "mother", role = "outcome",
      units = z$genes, alignment = z$al, covariates = "dcell",
      block = z$dat$cid, paired = TRUE, pair_focal = 1,
      baseline = FALSE, weighting = "flat", B = B, seed = seeds[R + r]))
    pa <- setNames(fit$p_adj, fit$unit)
    pp <- setNames(fit$p, fit$unit)
    global_any[r] <- any(pa < alpha, na.rm = TRUE)
    global_point[r, names(pp)] <- pp < alpha

    z <- .dcc_make(n_pairs = n_pairs, planted = TRUE,
                   seed = seeds[R + r])
    fit <- suppressMessages(dmsa_change(
      z$M0, z$M1, z$dat, exposure = "mother", role = "outcome",
      units = z$genes, alignment = z$al, covariates = "dcell",
      block = z$dat$cid, paired = TRUE, pair_focal = 1,
      baseline = FALSE, weighting = "flat", B = B, seed = seeds[r]))
    pa <- setNames(fit$p_adj, fit$unit)
    partial_power[r] <- isTRUE(pa["g1"] < alpha)
    partial_any_null[r] <- any(pa[setdiff(names(pa), "g1")] < alpha,
                               na.rm = TRUE)

    if (r %% 25L == 0L)
      message("replicate ", r, "/", R)
  }

  se <- function(p) sqrt(p * (1 - p) / R)
  gfw <- mean(global_any)
  sfw <- mean(partial_any_null)
  pow <- mean(partial_power)
  point <- colMeans(global_point)

  summary <- data.frame(
    quantity = c("global-null FWER", "partial-null FWER", "planted g1 power",
                 paste0("pointwise null ", names(point))),
    estimate = c(gfw, sfw, pow, point),
    mcse = c(se(gfw), se(sfw), se(pow), vapply(point, se, numeric(1))),
    target = c(alpha, alpha, NA, rep(alpha, length(point))),
    stringsAsFactors = FALSE
  )

  list(summary = summary,
       global_null_any = global_any,
       partial_null_any = partial_any_null,
       planted_power = partial_power,
       pointwise_null = global_point,
       settings = list(R = R, B = B, n_pairs = n_pairs,
                       seed = seed, alpha = alpha))
}
