## ---------------------------------------------------------------------------
## The vector core of dmsa_modmed(): any number of predictors X and of
## parallel mediators M (one moderator W or none). One M-model per mediator
## (M_j ~ all X + W + X:W ...), one Y-model (Y ~ all M + all X + W + M:W ...),
## the same corrected moments and the same Freedman-Lane / HC pivot as the
## scalar core (.mm_core), which it reproduces when x and m are scalars
## (tested). Estimands: every X_i -> M_j -> Y product by joint significance
## with Holm over the products (the declared family), the direct effects,
## and with two or more LATENT mediators the connected measurement of the
## mediators (MATH 1.7: the correlated-error covariance theta_MM of the two
## Bartlett scores from a joint two-factor fit, subtracted from the (M_i, M_j)
## cell of the corrected moments when |theta| / sqrt(V_i V_j) > .10).
## Engine untouched: .med_moments / .mm_moments_sq / .med_reg / .med_fl.
## ---------------------------------------------------------------------------

.mm_core_v <- function(meas, x, m, w, y, w_on, direct_moderation, quadratic, covM, covn,
                       B, blocks, block_mode, conditional_at, alpha_ssc, pivot = "hc", strata = NULL,
                       theta_mm = NULL) {
  n <- length(meas[[y]]$f)
  nm <- c(x, m, if (!is.null(w)) w, y)
  F <- vapply(nm, function(v) meas[[v]]$f - mean(meas[[v]]$f), numeric(n)); colnames(F) <- nm
  voff <- vapply(nm, function(v) meas[[v]]$V, numeric(1)); names(voff) <- nm
  voff <- voff[voff > 0]
  prods <- list(); D <- F
  addcol <- function(D, v, pn) { D <- cbind(D, v); colnames(D)[ncol(D)] <- pn; D }
  xw <- character(); mw_ <- character()
  if (!is.null(w)) {
    if (w_on %in% c("a", "both") || direct_moderation) for (xi in x) {
      pn <- paste0(xi, ":", w); D <- addcol(D, F[, xi] * F[, w], pn); prods[[pn]] <- c(xi, w); xw <- c(xw, pn) }
    if (w_on %in% c("b", "both")) for (mj in m) {
      pn <- paste0(mj, ":", w); D <- addcol(D, F[, mj] * F[, w], pn); prods[[pn]] <- c(mj, w); mw_ <- c(mw_, pn) }
  }
  if (quadratic) for (v in unique(c(x, w, m))) {
    pn <- paste0(v, "^2"); D <- addcol(D, F[, v]^2, pn); prods[[pn]] <- c(v, v) }
  if (!is.null(covM)) { C <- scale(as.matrix(covM), scale = FALSE); D <- cbind(D, C)
    colnames(D)[seq(ncol(D) - ncol(C) + 1, ncol(D))] <- covn }
  D <- scale(D, scale = FALSE)
  mm <- if (quadratic) .mm_moments_sq(D, voffset = voff, products = prods, n = n, alpha = alpha_ssc) else
    .med_moments(D, voffset = voff, products = prods, n = n, alpha = alpha_ssc)
  ## connected measurement of the mediators (MATH 1.7)
  theta_applied <- list()
  if (!is.null(theta_mm)) for (k in seq_along(theta_mm)) {
    th <- theta_mm[[k]]; pr <- th$pair
    if (!isTRUE(th$apply)) next
    S <- mm$S; S[pr[1], pr[2]] <- S[pr[1], pr[2]] - th$theta; S[pr[2], pr[1]] <- S[pr[1], pr[2]]
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) > 1e-6 * max(ev)) { mm$S <- S; theta_applied[[paste(pr, collapse = ",")]] <- th$theta }
    else theta_applied[[paste(pr, collapse = ",")]] <- NA_real_
  }
  Vof <- function(v) if (v %in% names(voff)) voff[[v]] else 0
  a_terms <- c(x, if (!is.null(w)) w, if (w_on %in% c("a", "both")) xw)
  y_terms <- c(m, x, if (!is.null(w)) w, if (w_on %in% c("b", "both")) mw_, if (direct_moderation) xw)
  sq_m <- if (quadratic) paste0(unique(c(x, w)), "^2") else character()
  sq_y <- if (quadratic) paste0(unique(c(m, x, w)), "^2") else character()
  mx <- c(a_terms, sq_m, covn); yx <- c(y_terms, sq_y, covn)
  fl <- function(xn, yn, tm, Vy) .med_fl(D, mm$S, xn, yn, tm, n, B, blocks, Vy = Vy, block_mode = block_mode, Sobs = mm$Sobs, pivot = pivot, strata = strata)
  paths <- list()
  ## the scalar call order (a-terms, b[m], the direct terms, b[m:w], c'[x:w], squares) per mediator
  for (mj in m) for (tm in a_terms) paths[[paste0("a[", tm, if (length(m) > 1L) paste0("->", mj) else "", "]")]] <- fl(mx, mj, tm, Vof(mj))
  for (mj in m) paths[[paste0("b[", mj, "]")]] <- fl(yx, y, mj, Vof(y))
  for (tm in c(x, if (!is.null(w)) w)) paths[[paste0("c'[", tm, "]")]] <- fl(yx, y, tm, Vof(y))
  if (w_on %in% c("b", "both")) for (pn in mw_) paths[[paste0("b[", pn, "]")]] <- fl(yx, y, pn, Vof(y))
  if (direct_moderation) for (pn in xw) paths[[paste0("c'[", pn, "]")]] <- fl(yx, y, pn, Vof(y))
  if (quadratic) for (mj in m) for (tm in sq_m) paths[[paste0("a[", tm, if (length(m) > 1L) paste0("->", mj) else "", "]")]] <- fl(mx, mj, tm, Vof(mj))
  g <- function(k) paths[[k]]$b; gp <- function(k) paths[[k]]$p
  ak <- function(tm, mj) paste0("a[", tm, if (length(m) > 1L) paste0("->", mj) else "", "]")
  sdw <- if (!is.null(w)) sqrt(max(mm$S[w, w], 0)) else NA_real_
  ## estimands per X_i -> M_j pair
  idx <- list(); cond <- list(); est <- c()
  for (xi in x) for (mj in m) {
    tag <- if (length(x) > 1L || length(m) > 1L) paste0("[", xi, "->", mj, "]") else ""
    a1 <- g(ak(xi, mj)); b1 <- g(paste0("b[", mj, "]"))
    a3 <- if (!is.null(w) && w_on %in% c("a", "both")) g(ak(paste0(xi, ":", w), mj)) else 0
    b3 <- if (!is.null(w) && w_on %in% c("b", "both")) g(paste0("b[", mj, ":", w, "]")) else 0
    theta <- function(Wsd) if (is.null(w)) rep(a1 * b1, length(Wsd)) else (a1 + a3 * Wsd * sdw) * (b1 + b3 * Wsd * sdw)
    p_a1 <- gp(ak(xi, mj)); p_b1 <- gp(paste0("b[", mj, "]"))
    idx[[paste0("IE", tag)]] <- list(est = a1 * b1, p_js = max(p_a1, p_b1), kind = "a1*b1 (indirect effect at W = 0)", parts = c(ak(xi, mj), paste0("b[", mj, "]")), pair = c(xi, mj))
    if (!is.null(w)) {
      pa <- if (w_on %in% c("a", "both")) gp(ak(paste0(xi, ":", w), mj)) else NA_real_
      pb <- if (w_on %in% c("b", "both")) gp(paste0("b[", mj, ":", w, "]")) else NA_real_
      if (w_on == "a") idx[[paste0("IMM", tag)]] <- list(est = a3 * b1, p_js = max(pa, p_b1), kind = "index of moderated mediation a3*b1 (first stage)", parts = c(ak(paste0(xi, ":", w), mj), paste0("b[", mj, "]")), pair = c(xi, mj))
      if (w_on == "b") idx[[paste0("IMM", tag)]] <- list(est = a1 * b3, p_js = max(pb, p_a1), kind = "index of moderated mediation a1*b3 (second stage)", parts = c(ak(xi, mj), paste0("b[", mj, ":", w, "]")), pair = c(xi, mj))
      if (w_on == "both") {
        pj <- c(I_a = max(pa, p_b1), I_b = max(pb, p_a1), I_ab = max(pa, pb)); ph <- stats::p.adjust(pj, "holm")
        idx[[paste0("IMM", tag)]] <- list(est = theta(1) - theta(-1), p_js = min(ph), products = c(I_a = a3 * b1, I_b = a1 * b3, I_ab = a3 * b3), p_products = pj, p_holm = ph,
                                          kind = "dual stage: delta theta = theta(+1 SD) - theta(-1 SD); test = Holm union over the joint-significance p of I_a, I_b, I_ab", pair = c(xi, mj))
      }
    }
    cond[[length(cond) + 1L]] <- if (!is.null(w)) data.frame(x = xi, m = mj, w_sd = conditional_at, indirect = theta(conditional_at)) else data.frame(x = xi, m = mj, w_sd = 0, indirect = a1 * b1)
    sfx <- if (nzchar(tag)) tag else ""
    est <- c(est, stats::setNames(c(a1, a3, b1, b3, a1 * b1, theta(-1), theta(1), theta(1) - theta(-1)), paste0(c("a1", "a3", "b1", "b3", "IE", "theta_lo", "theta_hi", "dtheta"), sfx)))
  }
  ## Holm over the declared family: the IE products across pairs, the IMM indices across pairs
  for (fam in c("IE", "IMM")) { k <- grep(paste0("^", fam), names(idx)); if (length(k) > 1L) {
    ph <- stats::p.adjust(vapply(idx[k], function(i) i$p_js, numeric(1)), "holm")
    for (j in seq_along(k)) idx[[k[j]]]$p_holm_family <- ph[[j]] } }
  for (i in seq_along(idx)) { if (is.null(idx[[i]]$p_holm)) idx[[i]]$p_holm <- idx[[i]]$p_js; if (is.null(idx[[i]]$p_holm_family)) idx[[i]]$p_holm_family <- idx[[i]]$p_js }
  cond <- do.call(rbind, cond)
  if (length(x) == 1L && length(m) == 1L) cond <- cond[, c("w_sd", "indirect")]
  list(paths = paths, indices = idx, conditional = cond, moments = mm, n = n,
       terms = list(M_model = mx, Y_model = yx),
       collinearity = list(M_model = .med_collin(mm$S, mx), Y_model = .med_collin(mm$S, yx)),
       estimates = est, theta_applied = theta_applied)
}

## MATH 1.7: theta between two Bartlett scores from a CONNECTED two-factor
## fit (lavaan, Suggests): w_a' R_ab w_b in the unit-latent-variance metric.
.mm_theta_connected <- function(Ta, Tb, ma, mb, pair) {
  out <- list(pair = pair, theta = NA_real_, ratio = NA_real_, phi = NA_real_, apply = FALSE, note = "")
  if (ma$model == "single" || mb$model == "single") { out$note <- "a mediator is at the single rung: no factor model to connect"; return(out) }
  if (!requireNamespace("lavaan", quietly = TRUE)) { out$note <- "lavaan not installed: connected measurement not fitted"; return(out) }
  Ta <- as.matrix(Ta); Tb <- as.matrix(Tb); colnames(Ta) <- paste0("a", seq_len(ncol(Ta))); colnames(Tb) <- paste0("b", seq_len(ncol(Tb)))
  dat <- as.data.frame(cbind(Ta, Tb))
  mod <- sprintf("Fa =~ %s\nFb =~ %s", paste(colnames(Ta), collapse = " + "), paste(colnames(Tb), collapse = " + "))
  fit <- tryCatch(suppressWarnings(lavaan::cfa(mod, data = dat, std.lv = TRUE)), error = function(e) NULL)
  if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit, "converged"))) { out$note <- "joint CFA did not converge"; return(out) }
  est <- lavaan::lavInspect(fit, "est"); lam <- est$lambda; phi <- est$psi["Fb", "Fa"]
  S <- stats::cov(dat); ia <- colnames(Ta); ib <- colnames(Tb)
  R_ab <- S[ia, ib] - lam[ia, "Fa"] %*% t(lam[ib, "Fb"]) * phi
  w_of <- function(mm, Z) { a <- mm$lambda / mm$theta; w <- a / sum(mm$lambda * a); f0 <- as.numeric(scale(Z, scale = FALSE) %*% w); V0 <- 1 / sum(mm$lambda * a)
    s2 <- max(stats::var(f0) - V0, 0.05 * stats::var(f0)); w / sqrt(s2) }
  wa <- w_of(ma, Ta); wb <- w_of(mb, Tb)
  th <- as.numeric(t(wa) %*% R_ab %*% wb)
  out$theta <- th; out$ratio <- th / sqrt(ma$V * mb$V); out$phi <- phi
  out$fit <- c(cfi = unname(lavaan::fitMeasures(fit, "cfi")), rmsea = unname(lavaan::fitMeasures(fit, "rmsea")))
  out$apply <- is.finite(out$ratio) && abs(out$ratio) > 0.10
  out$note <- if (out$apply) "theta_MM applied (|theta| / sqrt(V_a V_b) > .10)" else "theta_MM below the .10 rule: not applied"
  out
}
