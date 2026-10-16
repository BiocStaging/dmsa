# ============================================================================
# THE TWO-CALL INTERFACE, CALL 2: dmsa_report()
#
# Runs the frame and delivers: figures/ tables/ summary.md under
# frame$outdir. Statistics are the ones locked elsewhere in the package -
# dmsa_triangulate() at every level (three lenses, one permutation stream,
# maxT/minP within the LEVEL-LOCAL family only), the composite one-model for
# moderation (the one lens that carries to products), the Westfall-Young
# family logic throughout. No new statistics live in this file.
#
# Level-local families, exactly as specified:
#   system : the named systems are the family (1 system -> no correction)
#   module : the modules of that system
#   gene   : the genes of that system
#   probe  : the probes of that gene (run for surviving genes only)
#
# Brain-bridge post-hoc (always on, never fatal): after selection, the union
# of surviving units' probes goes through cpgdirection::cpg_brain_bridge();
# hit rule FIXED as: bridge_usable | has_t2_direction | !is.na(bridge_grade).
# ensemble_pctile is reported as context, never as a hit.
# ============================================================================

## Spec 50: `frame$outcome_type` carries one family per outcome when the frame
## declares several. A single-outcome frame still stores the bare string, so
## this accessor reads both shapes and always answers for the outcome in hand.
.rp_otype <- function(frame, oc) {
  ot <- frame$outcome_type
  if (is.null(ot) || !length(ot)) return("gaussian")
  nm <- names(ot)
  if (!is.null(nm) && oc %in% nm) return(unname(ot[[oc]]))
  unname(ot[[1L]])
}

.rp_pal <- function(frame, n = 6) {
  pal <- frame$palette
  if (requireNamespace("viridisLite", quietly = TRUE) &&
      pal %in% c("viridis", "magma", "plasma", "cividis", "inferno", "mako",
                 "rocket", "turbo"))
    return(do.call(getExportedValue("viridisLite", pal), list(n, end = .92)))
  grDevices::hcl.colors(n, "viridis")
}

.rp_dev <- function(frame, file, width = 7, height = 5) {
  if (frame$plot_type == "pdf")
    grDevices::pdf(paste0(file, ".pdf"), width = width, height = height)
  else grDevices::png(paste0(file, ".png"), width = width * 300,
                      height = height * 300, res = 300)
  paste0(file, ".", frame$plot_type)
}

## ---- design pieces ---------------------------------------------------------
.rp_rhs <- function(frame, oc) {
  others <- setdiff(frame$outcome, oc)
  Filter(nzchar, c(oc, others, frame$covariates, frame$chip))
}

## per-probe slopes for one unit (for the locus panel)
.rp_probe_fits <- function(frame, oc, cols) {
  ff <- stats::as.formula(paste("~", paste(.rp_rhs(frame, oc), collapse = "+")))
  X <- stats::model.matrix(ff, frame$data)
  fi <- which(colnames(X) == oc)
  XtXi <- solve(crossprod(X)); H <- XtXi %*% t(X)
  Y <- frame$M[, cols, drop = FALSE]
  Bh <- H %*% Y; R <- Y - X %*% Bh
  s2 <- colSums(R^2) / (nrow(X) - ncol(X))
  b <- Bh[fi, ]; se <- sqrt(s2 * XtXi[fi, fi])
  ## nominal two-sided per-probe p from the SAME fit the panel's b/se come
  ## from. The locus figure greys probes above alpha so a reader can see at a
  ## glance that a gene-level hit need not rest on any single CpG (PI,
  ## 2026-08-29). Nominal, not adjusted: an adjusted probe p almost never
  ## fires, which would grey every probe of a gene selected precisely because
  ## its CpGs agree.
  list(b = b, se = se,
       p = 2 * stats::pt(-abs(b / se), df = nrow(X) - ncol(X)))
}

## one gene-level (or module/probe/system-level) triangulation for one system
.report_gene_level <- function(frame, oc, sid, B = frame$B, quiet = FALSE,
                               units = NULL, alignment = NULL, Ms = NULL) {
  s <- frame$sets[[as.character(sid)]]
  if (is.null(Ms)) Ms <- frame$M[, s$columns, drop = FALSE]
  if (is.null(units)) units <- s$map$gene
  if (is.null(alignment)) alignment <- s$alignment
  r <- dmsa_triangulate(Ms, frame$data, .rp_rhs(frame, oc), oc, units,
                        alignment, block = frame$block, B = B,
                        correction = frame$correction,
                        weighting = frame$weighting %||% "combined",
                        w_floor = frame$w_floor %||% 1.5, seed = frame$seed,
                        ri_group = if (.nz(frame$chip_random))
                                     frame$data[[frame$chip_random]] else NULL)
  r$system_id <- s$system_id; r$system <- s$system; r$outcome <- oc
  ## realized family-wise error of the pre-registered any-lens rule IN THIS FAMILY,
  ## measured from the run's own permutation stream (share of null draws
  ## whose best adjusted p anywhere in the family clears alpha) - the
  ## design-specific number the MS's ".04-.12" generalises
  .mb <- attr(r, "union_null_min")
  r$fwer_realized <- if (length(.mb) && any(is.finite(.mb)))
    mean(.mb < frame$alpha, na.rm = TRUE) else NA_real_
  ## the null false-alarm rate of ordinary BH at the naming level in THIS
  ## family (PI ruling 2026-10-15): the naming rule itself is calibrated to
  ## fire in a fraction fdr of null families; this shows what the uncalibrated
  ## BH would have done under the family's measured dependence
  .qm <- attr(r, "unit_null_qmin")
  r$null_alarm_bh <- if (length(.qm) && any(is.finite(.qm)))
    mean(.qm <= (frame$fdr %||% .10), na.rm = TRUE) else NA_real_
  r
}

## composite-only moderated test for a set of units within one family
.report_moderated <- function(frame, oc, units_list, family_label) {
  d0 <- frame$data
  mc <- scale(d0[[frame$mod]])[, 1]
  m2c <- if (nzchar(frame$mod2)) scale(d0[[frame$mod2]])[, 1] else NULL
  res <- list(); nulls <- list(); nulls_q <- list()
  for (u in names(units_list)) {
    cols <- units_list[[u]]$columns
    al <- units_list[[u]]$alignment
    S <- dmsa_scores(frame$M[, cols, drop = FALSE], al,
                     weighting = frame$weighting %||% "combined",
                     w_floor = frame$w_floor %||% 1.5)$aligned
    d2 <- d0; d2$S <- as.numeric(scale(S)); d2$.mc <- mc
    if (!is.null(m2c)) d2$.m2c <- m2c
    d2$.sq <- d2$S^2 - mean(d2$S^2, na.rm = TRUE)
    covs <- Filter(nzchar, c(setdiff(frame$outcome, oc), frame$covariates, frame$chip))
    if (frame$frame_role == "predictor") {
      lhs <- oc
      inter <- if (is.null(m2c)) "S * .mc" else "S * .mc * .m2c"
      term <- if (is.null(m2c)) "S:.mc" else "S:.mc:.m2c"
    } else {
      lhs <- "S"
      inter <- if (is.null(m2c)) paste(oc, "* .mc") else
        paste(oc, "* .mc * .m2c")
      term <- if (is.null(m2c)) paste0(oc, ":.mc") else
        paste0(oc, ":.mc:.m2c")
    }
    ff <- stats::as.formula(paste(lhs, "~", paste(c(inter, covs),
                                                  collapse = " + ")))
    f <- dmsa_model(ff, d2, term, block = frame$block, B = frame$B,
                    seed = frame$seed, nulls = TRUE)
    ## MODERATED NON-LINEAR. type != "linear" under moderation = TRUE used to be
    ## silently ignored: the moderation branch returns before the shape scan
    ## runs, so the user asked for a moderated non-linear model and received a
    ## moderated LINEAR one with no warning. The question it should answer is
    ## whether the CURVATURE itself depends on the moderator, so the tested term
    ## is the moderator by squared-score product, with every lower-order term
    ## present (S, S^2, moderator, S x moderator) for marginality.
    fq <- NULL
    if (frame$type != "linear") {
      iq <- if (is.null(m2c)) "S * .mc + .sq * .mc"
            else "S * .mc * .m2c + .sq * .mc * .m2c"
      ffq <- stats::as.formula(paste(lhs, "~", paste(c(iq, covs),
                                                     collapse = " + ")))
      ## Resolve the tested term from the DESIGN rather than assuming how R
      ## orders an interaction's name: for `S * .mc + .sq * .mc` the column is
      ## `.mc:.sq`, not `.sq:.mc`, and guessing wrong makes dmsa_model() throw -
      ## which, swallowed by a tryCatch, would silently produce NA.
      want <- if (is.null(m2c)) c(".sq", ".mc") else c(".sq", ".mc", ".m2c")
      cn <- .rp_try("moderation: quadratic design", colnames(stats::model.matrix(
              ffq, stats::model.frame(ffq, d2, na.action = stats::na.omit))),
              outcome = oc, fallback = "quadratic moderation term not tested")
      tq <- NULL
      if (!is.null(cn)) {
        parts <- strsplit(cn, ":", fixed = TRUE)
        ok <- vapply(parts, function(z) length(z) == length(want) &&
                       setequal(z, want), logical(1))
        if (any(ok)) tq <- cn[which(ok)[1L]]
      }
      if (is.null(tq)) {
        warning("moderated non-linear term (", paste(want, collapse = " x "),
                ") is not estimable for unit '", u, "' - reported as NA",
                call. = FALSE)
      } else {
        fq <- tryCatch(dmsa_model(ffq, d2, tq, block = frame$block,
                                  B = frame$B, seed = frame$seed, nulls = TRUE),
                       error = function(e) {
                         warning("moderated non-linear test failed for unit '",
                                 u, "': ", conditionMessage(e), call. = FALSE)
                         NULL })
      }
    }
    res[[u]] <- data.frame(unit = u, n_probes = sum(al$s != 0, na.rm = TRUE),
                           b = f$b, t = f$t, p_composite = f$p_perm,
                           b_curv = if (is.null(fq)) NA_real_ else fq$b,
                           t_curv = if (is.null(fq)) NA_real_ else fq$t,
                           p_curv = if (is.null(fq)) NA_real_ else fq$p_perm,
                           stringsAsFactors = FALSE)
    nulls[[u]] <- f$null_t
    nulls_q[[u]] <- if (is.null(fq)) NULL else fq$null_t
  }
  out <- do.call(rbind, res)
  ## family maxT on |t| across units (shared permutation stream: same seed,
  ## same blocks, same rows -> the b-th draw is the same reordering for all)
  TN <- do.call(cbind, nulls)
  if (!is.null(TN) && ncol(TN) > 1) {
    mx <- apply(abs(TN), 1, max, na.rm = TRUE)
    out$p_composite_adj <- vapply(abs(out$t), function(o)
      (1 + sum(mx >= o, na.rm = TRUE)) / (length(mx) + 1), numeric(1))
    o <- order(out$p_composite); out$p_composite_adj[o] <-
      cummax(out$p_composite_adj[o])
  } else out$p_composite_adj <- out$p_composite
  ## A family of one has nothing to correct against, so its "adjusted" p is just
  ## its raw p. Presenting that as family-adjusted makes an uncorrected number
  ## look like the strongest result in the run.
  out$n_family <- length(units_list)
  out$corrected <- length(units_list) > 1L
  ## the moderated-curvature term gets the same level-local maxT
  TQ <- if (length(nulls_q) == nrow(out) &&
            !any(vapply(nulls_q, is.null, logical(1))))
    do.call(cbind, nulls_q) else NULL
  if (!is.null(TQ) && ncol(TQ) > 1L && any(is.finite(out$t_curv))) {
    mq <- apply(abs(TQ), 1, max, na.rm = TRUE)
    out$p_curv_adj <- vapply(abs(out$t_curv), function(o)
      if (!is.finite(o)) NA_real_ else
        (1 + sum(mq >= o, na.rm = TRUE)) / (length(mq) + 1), numeric(1))
    ok <- which(is.finite(out$p_curv_adj))
    if (length(ok)) { oq <- ok[order(out$p_curv[ok])]
                      out$p_curv_adj[oq] <- cummax(out$p_curv_adj[oq]) }
  } else out$p_curv_adj <- out$p_curv
  out$family <- family_label
  out
}

## shape scan: non-linear / exponential arms on the aligned score (ACAT)
.report_shape <- function(frame, oc, cols, al, nulls = FALSE) {
  S <- as.numeric(scale(dmsa_scores(frame$M[, cols, drop = FALSE], al,
                        weighting = frame$weighting %||% "combined",
                        w_floor = frame$w_floor %||% 1.5)$aligned))
  covs <- Filter(nzchar, c(setdiff(frame$outcome, oc), frame$covariates, frame$chip))
  d <- frame$data
  d$S <- S
  d$.sq  <- S^2 - mean(S^2, na.rm = TRUE)
  d$.thr <- as.numeric(S > 0)
  if (frame$type == "exponential") d$.ex <- as.numeric(scale(exp(pmin(S, 5))))
  rhs <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""
  B <- max(499L, frame$B %/% 4L)

  ## MARGINALITY. Every departure term is tested with the LINEAR term in the
  ## model, so the test is the incremental one - "is there a quadratic component
  ## over and above the linear component". Substituting the squared score FOR
  ## the linear score tests a curve whose turning point is pinned at the
  ## centring value, which is not evidence of curvature at all.
  fit1 <- function(extra, term) {
    ff <- stats::as.formula(paste(oc, "~ S",
                                  if (nzchar(extra)) paste("+", extra) else "",
                                  rhs))
    .rp_try("moderation: model", dmsa_model(ff, d, term, block = frame$block,
                        B = B, seed = frame$seed, nulls = nulls),
            outcome = oc, fallback = "unit dropped from the moderation table")
  }
  gp <- function(f) if (is.null(f)) NA_real_ else f$p_perm
  f_lin  <- fit1("", "S")
  f_quad <- fit1(".sq",  ".sq")
  f_thr  <- fit1(".thr", ".thr")
  f_expo <- if (frame$type == "exponential") fit1(".ex", ".ex") else NULL
  p_lin <- gp(f_lin); p_quad <- gp(f_quad); p_thr <- gp(f_thr)
  p_expo <- gp(f_expo)

  ## ACAT over the DEPARTURE terms only. Pooling the linear arm in would make
  ## the omnibus "is there any association", which a significant linear term
  ## alone can drive and which therefore cannot support a claim about shape.
  dep <- c(p_quad, p_thr, p_expo); dep <- dep[is.finite(dep)]
  p_nl <- if (!length(dep)) NA_real_ else {
    q <- pmin(pmax(dep, 1e-15), 1 - 1e-15)
    0.5 - atan(mean(tan((0.5 - q) * pi))) / pi
  }

  ## SHAPE, not curvature. Lind & Mehlum (2010): H1 is slope(xl) < 0 AND
  ## slope(xh) > 0 - opposite signs at the two ENDPOINTS. It is an
  ## intersection-union (Sasabuchi 1980) test, so the p-value is the LARGER of
  ## the two one-sided p-values and both ends must clear alpha. A significant
  ## squared term alone is not enough: a convex but monotone relationship
  ## produces one, and a quadratic approximation then invents a turning point.
  u <- list(shape = NA_character_, p_ushape = NA_real_,
            p_ushape_trim = NA_real_, turn = NA_real_, turn_lo = NA_real_,
            turn_hi = NA_real_, turn_inside = NA, turn_ci_inside = NA,
            slope_lo = NA_real_, slope_hi = NA_real_)
  m2 <- try(stats::lm(stats::as.formula(paste(oc, "~ S + .sq", rhs)), d),
            silent = TRUE)
  if (!inherits(m2, "try-error")) {
    V <- stats::vcov(m2); cf <- stats::coef(m2)[rownames(V)]
    if (all(c("S", ".sq") %in% names(cf))) {
      Sk <- m2$model$S; xL <- min(Sk); xH <- max(Sk)
      cv <- function(x) { v <- stats::setNames(numeric(length(cf)), names(cf))
                          v["S"] <- 1; v[".sq"] <- 2 * x; v }
      sl  <- function(x) sum(cv(x) * cf)
      sse <- function(x) { v <- cv(x); sqrt(drop(v %*% V %*% v)) }
      dfr <- stats::df.residual(m2)
      sasa <- function(lo, hi) {
        tl <- sl(lo) / sse(lo); th <- sl(hi) / sse(hi)
        pu  <- max(stats::pt(tl, dfr),  stats::pt(-th, dfr))
        piu <- max(stats::pt(-tl, dfr), stats::pt(th, dfr))
        list(shape = if (pu <= piu) "U" else "inverted-U", p = min(pu, piu))
      }
      full <- sasa(xL, xH)
      ## Lind & Mehlum note the interval may also be taken in the interior of
      ## the domain "to make sure the shape is not only a marginal phenomenon".
      ## In a skewed score the endpoints ARE its sparsest points, so repeat on
      ## the central 90% and let a tail-borne shape show itself.
      qq <- stats::quantile(Sk, c(.05, .95), na.rm = TRUE)
      trim <- sasa(qq[1], qq[2])
      b1 <- unname(cf["S"]); b2 <- unname(cf[".sq"])
      xs <- if (is.finite(b2) && b2 != 0) -b1 / (2 * b2) else NA_real_
      ## Fieller (1943) interval for the extremum. Lind & Mehlum present the
      ## containment of this interval in [xl, xh] as the equivalent test, and
      ## warn the delta method is severely biased in finite samples. A <= 0
      ## means the set is unbounded - the extremum is simply not determined,
      ## and saying so is the point.
      tc <- stats::qt(.975, dfr)
      av <- -b1; bv <- 2 * b2
      Va <- V["S", "S"]; Vb <- 4 * V[".sq", ".sq"]; Cab <- -2 * V["S", ".sq"]
      A <- bv^2 - tc^2 * Vb
      Bq <- -2 * (av * bv - tc^2 * Cab)
      Cq <- av^2 - tc^2 * Va
      disc <- Bq^2 - 4 * A * Cq
      fl <- c(NA_real_, NA_real_)
      if (is.finite(disc) && disc >= 0 && A > 0)
        fl <- sort(c((-Bq - sqrt(disc)) / (2 * A), (-Bq + sqrt(disc)) / (2 * A)))
      else if (is.finite(A) && A <= 0) fl <- c(-Inf, Inf)
      u <- list(shape = full$shape, p_ushape = full$p,
                p_ushape_trim = trim$p, turn = xs,
                turn_lo = fl[1], turn_hi = fl[2],
                turn_inside = is.finite(xs) && xs > xL && xs < xH,
                turn_ci_inside = all(is.finite(fl)) && fl[1] > xL && fl[2] < xH,
                slope_lo = sl(xL), slope_hi = sl(xH))
    }
  }
  c(list(p_lin = p_lin, p_quad = p_quad, p_thr = p_thr, p_expo = p_expo,
         p_nonlin = p_nl,
         null_quad = if (nulls && !is.null(f_quad)) f_quad$null_t else NULL,
         null_thr  = if (nulls && !is.null(f_thr))  f_thr$null_t  else NULL,
         null_expo = if (nulls && !is.null(f_expo)) f_expo$null_t else NULL,
         t_quad = if (!is.null(f_quad)) f_quad$t else NA_real_,
         t_thr  = if (!is.null(f_thr))  f_thr$t  else NA_real_,
         t_expo = if (!is.null(f_expo)) f_expo$t else NA_real_), u)
}

## One level-local family of shape scans, with Westfall-Young maxT on the
## quadratic departure statistic. Everything else in DMSA is corrected inside
## its level-local family; the shape scan used to sit outside that entirely,
## reporting raw p-values, which is exactly the inconsistency that makes a
## scan across many units unreportable.
.report_shape_family <- function(frame, oc, units_list, family_label, level) {
  rows <- list(); nulls <- list()
  for (u in names(units_list)) {
    o <- .rp_try("shape scan: unit", .report_shape(frame, oc, units_list[[u]]$columns,
                                units_list[[u]]$alignment, nulls = TRUE),
                 outcome = oc, unit = u, fallback = "unit absent from the shape table")
    if (is.null(o)) next
    g <- function(nm) if (is.null(o[[nm]])) NA else unname(o[[nm]])
    rows[[u]] <- data.frame(
      outcome = oc, level = level, family = family_label, unit = u,
      n_probes = length(units_list[[u]]$columns),
      p_lin = as.numeric(g("p_lin")), p_quad = as.numeric(g("p_quad")),
      p_thr = as.numeric(g("p_thr")), p_expo = as.numeric(g("p_expo")),
      p_nonlin = as.numeric(g("p_nonlin")), t_quad = as.numeric(g("t_quad")),
      t_thr = as.numeric(g("t_thr")), t_expo = as.numeric(g("t_expo")),
      shape = as.character(g("shape")),
      p_ushape = as.numeric(g("p_ushape")),
      p_ushape_trim = as.numeric(g("p_ushape_trim")),
      turn = as.numeric(g("turn")), turn_lo = as.numeric(g("turn_lo")),
      turn_hi = as.numeric(g("turn_hi")),
      turn_inside = as.logical(g("turn_inside")),
      turn_ci_inside = as.logical(g("turn_ci_inside")),
      n_family = length(units_list),
      stringsAsFactors = FALSE)
    nulls[[u]] <- list(quad = o$null_quad, thr = o$null_thr, expo = o$null_expo)
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  ## maxT within this level-local family, for EVERY departure arm - not just
  ## the quadratic. An exponential run's own arm was being reported raw while
  ## the quadratic beside it was corrected, which is not a defensible pair of
  ## numbers to put in one table.
  for (arm in c("quad", "thr", "expo")) {
    praw <- out[[paste0("p_", arm)]]; traw <- out[[paste0("t_", arm)]]
    nl <- lapply(nulls, function(z) z[[arm]])
    TN <- if (length(nl) && !any(vapply(nl, is.null, logical(1))))
      do.call(cbind, nl) else NULL
    if (!is.null(TN) && ncol(TN) > 1L && any(is.finite(traw))) {
      mx <- apply(abs(TN), 1, max, na.rm = TRUE)
      adj <- vapply(abs(traw), function(o)
        if (!is.finite(o)) NA_real_ else
          (1 + sum(mx >= o, na.rm = TRUE)) / (length(mx) + 1), numeric(1))
      ok <- which(is.finite(adj))
      if (length(ok)) {
        o2 <- ok[order(praw[ok])]; adj[o2] <- cummax(adj[o2])
      }
      out[[paste0("p_", arm, "_adj")]] <- adj
    } else out[[paste0("p_", arm, "_adj")]] <- praw
  }
  out
}

## Shape figure: the aligned system score against the outcome it was tested on,
## with each arm's fit drawn over the partial residuals. `type = "non-linear"`
## is a claim about functional form, and a claim about functional form that
## ships as a p-value with no picture is not checkable.
.rp_fig_shape <- function(frame, oc, level, family, unit, sh, file) {
  ul <- try(.rp_units(frame, level), silent = TRUE)
  if (inherits(ul, "try-error") || is.null(ul[[family]])) return(invisible(NULL))
  u <- ul[[family]][[unit]]
  if (is.null(u)) return(invisible(NULL))
  s <- list(system = .rp_unit_label(frame, unit, level),
            columns = u$columns, alignment = u$alignment)
  S <- as.numeric(scale(dmsa_scores(frame$M[, s$columns, drop = FALSE],
                        s$alignment, weighting = frame$weighting %||% "combined",
                        w_floor = frame$w_floor %||% 1.5)$aligned))
  d <- frame$data; d$S <- S
  covs <- Filter(nzchar, c(setdiff(frame$outcome, oc), frame$covariates, frame$chip))
  rhs <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""
  ## Adjust the outcome for covariates once, then show every arm against the
  ## same adjusted outcome so the curves are comparable by eye.
  f0 <- try(stats::lm(stats::as.formula(paste(oc, "~ 1", rhs)), d), silent = TRUE)
  if (inherits(f0, "try-error")) return(invisible(NULL))
  keep <- match(rownames(f0$model), rownames(d))
  yadj <- as.numeric(stats::residuals(f0))
  Sk <- S[keep]
  ok <- is.finite(yadj) & is.finite(Sk)
  yadj <- yadj[ok]; Sk <- Sk[ok]
  if (length(Sk) < 20) return(invisible(NULL))

  pal <- .rp_pal(frame, 6); CL <- pal[1]; CQ <- pal[3]; CT <- pal[5]
  CE <- pal[2]
  fp <- .rp_dev(frame, file, width = 7.6, height = 5.4)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(4.4, 4.6, 4.4, 1.6))
  graphics::plot(Sk, yadj, pch = 19, cex = .55,
                 col = grDevices::adjustcolor("grey25", .40),
                 xlab = sprintf("%s aligned tone score (z)", s$system),
                 ylab = sprintf("%s | adjusted for covariates", .lab(frame, oc)),
                 cex.lab = .8, cex.axis = .75, col.axis = "grey30", bty = "n")
  graphics::abline(h = 0, col = "grey88")
  gx <- seq(min(Sk), max(Sk), length.out = 200)
  ## A quadratic's tails are set by whatever few points sit out there. Draw the
  ## fits SOLID across the central 90% of the score and DASHED beyond it, so
  ## curvature carried by four observations cannot be mistaken for the finding.
  qlim <- stats::quantile(Sk, c(.05, .95), na.rm = TRUE)
  drawfit <- function(fx, col) {
    inr <- gx >= qlim[1] & gx <= qlim[2]
    graphics::lines(gx[inr], fx[inr], col = col, lwd = 2.4)
    for (side in list(gx <= qlim[1], gx >= qlim[2]))
      if (sum(side) > 1) graphics::lines(gx[side], fx[side], col = col,
                                         lwd = 1.6, lty = 2)
  }
  ## Every drawn fit is the SAME hierarchical model its p-value came from -
  ## linear term always present. Drawing a curve from one model and labelling it
  ## with another model's p is exactly the mismatch this panel has to avoid.
  bl <- stats::lm(yadj ~ Sk)
  drawfit(stats::predict(bl, data.frame(Sk = gx)), CL)
  bq <- stats::lm(yadj ~ Sk + I(Sk^2))
  drawfit(stats::predict(bq, data.frame(Sk = gx)), CQ)
  bt <- stats::lm(yadj ~ Sk + I(Sk > 0))
  drawfit(stats::predict(bt, data.frame(Sk = gx)), CT)
  if (isTRUE(is.finite(sh$p_expo))) {
    be <- stats::lm(yadj ~ Sk + I(as.numeric(scale(exp(pmin(Sk, 5))))))
    ge <- (exp(pmin(gx, 5)) - mean(exp(pmin(Sk, 5)))) / stats::sd(exp(pmin(Sk, 5)))
    drawfit(as.numeric(stats::coef(be)[1] + stats::coef(be)[2] * gx +
                       stats::coef(be)[3] * ge), CE)
  }
  ## turning point of the quadratic fit, and whether it is inside the data
  if (is.finite(sh$turn) && sh$turn > min(Sk) && sh$turn < max(Sk))
    graphics::abline(v = sh$turn, lty = 3, col = "firebrick", lwd = 1.6)
  graphics::abline(v = 0, lty = 3, col = "grey60")
  usr <- graphics::par("usr")
  graphics::segments(Sk, usr[3], Sk, usr[3] + .035 * diff(usr[3:4]),
                     col = grDevices::adjustcolor("grey15", .35), lwd = 1)
  pf <- function(x) if (is.null(x) || !is.finite(x)) "n/a" else sprintf("%.4f", x)
  ## One legend, not two: a second box lands on the density strip.
  leg_txt <- c(sprintf("linear               p = %s", pf(sh$p_lin)),
               sprintf("+ quadratic (incr.)  p = %s  adj %s",
                       pf(sh$p_quad), pf(sh$p_quad_adj)),
               sprintf("+ threshold (incr.)  p = %s  adj %s",
                       pf(sh$p_thr), pf(sh$p_thr_adj)))
  leg_col <- c(CL, CQ, CT)
  if (isTRUE(is.finite(sh$p_expo))) {
    leg_txt <- c(leg_txt, sprintf("+ exponential (incr.) p = %s  adj %s",
                                  pf(sh$p_expo), pf(sh$p_expo_adj)))
    leg_col <- c(leg_col, CE)
  }
  graphics::legend("topleft", bty = "n", cex = .62, lwd = 2.4,
                   legend = leg_txt, col = leg_col)
  fci <- if (!is.finite(sh$turn)) "" else
    if (!is.finite(sh$turn_lo) || !is.finite(sh$turn_hi))
      sprintf("; extremum %.2f, Fieller interval unbounded", sh$turn)
    else sprintf("; extremum %.2f, Fieller 95%% [%.2f, %.2f]%s", sh$turn,
                 sh$turn_lo, sh$turn_hi,
                 if (isTRUE(sh$turn_ci_inside)) " (inside the data)"
                 else " (NOT contained in the data range)")
  ## The panel carries the NUMBERS; the explanation of what they mean lives in
  ## summary.md. A figure that has to be read as a paragraph is not a figure.
  shp <- if (!is.finite(sh$p_ushape)) "shape test n/a" else
    if (sh$p_ushape < 0.05 && isTRUE(sh$turn_inside))
      sprintf("%s: Sasabuchi p = %s%s", sh$shape, pf(sh$p_ushape),
              if (is.finite(sh$p_ushape_trim) && sh$p_ushape_trim >= .05)
                sprintf(" (%s on central 90%%)", pf(sh$p_ushape_trim)) else "")
    else sprintf("no %s (Sasabuchi p = %s)", sh$shape, pf(sh$p_ushape))
  ext <- if (!is.finite(sh$turn)) "" else
    if (!is.finite(sh$turn_lo) || !is.finite(sh$turn_hi))
      sprintf("   extremum %.2f, Fieller unbounded", sh$turn)
    else sprintf("   extremum %.2f, Fieller [%.2f, %.2f]", sh$turn,
                 sh$turn_lo, sh$turn_hi)
  .ttl <- sprintf("%s (%s) -> %s", .rp_unit_label(frame, unit, level), level,
                  .lab(frame, oc))
  .ft <- .rp_fit(.ttl, cex = .78, width = graphics::par("pin")[1])
  for (k in seq_along(.ft$lines))
    graphics::mtext(.ft$lines[k], side = 3,
                    line = 2.6 + (length(.ft$lines) - k) * .9,
                    adj = 0, font = 2, cex = .ft$cex)
  .fs <- .rp_fit(sprintf("%s%s", shp, ext), cex = .6,
                 width = graphics::par("pin")[1])
  for (k in seq_along(.fs$lines))
    graphics::mtext(.fs$lines[k], side = 3, line = 1.5 - (k - 1) * .7,
                    adj = 0, cex = .fs$cex, col = "grey25")
  graphics::mtext(sprintf(
    "solid = central 90%% of the score, dashed = extrapolation (%d of %d outside); ticks = data",
    sum(Sk < qlim[1] | Sk > qlim[2]), length(Sk)),
    side = 3, line = .55, adj = 0, cex = .55, col = "grey45")
  invisible(fp)
}

## ---- brain bridge (always on, never fatal) --------------------------------
.report_bridge <- function(probes) {
  if (!length(probes)) return(list(status = "no surviving probes", table = NULL))
  br <- .rp_try("brain bridge", {
    .bridge <- .cpgd("cpg_brain_bridge")
    if (is.null(.bridge)) stop("cpgdirection not installed")
    as.data.frame(.bridge(unique(probes), tissue = "blood"))
  }, fallback = "bridge line reads 'bridge check unavailable'")
  if (is.null(br)) return(list(status = "bridge check unavailable", table = NULL))
  hit <- vapply(seq_len(nrow(br)), function(i)
    isTRUE(br$bridge_usable[i]) || isTRUE(br$has_t2_direction[i]) ||
      !is.na(br$bridge_grade[i]), logical(1))
  br$bridge_hit <- hit
  list(status = if (any(hit))
    sprintf("%d of %d probes carry a formal brain bridge", sum(hit), nrow(br))
    else sprintf(
      "checked %d probes: no formal bridge (ensemble percentiles %s)",
      nrow(br), if (!is.null(br$ensemble_pctile) &&
                    any(is.finite(br$ensemble_pctile)))
        paste(range(round(br$ensemble_pctile), na.rm = TRUE),
              collapse = "-") else "none recorded"),
    table = br, any_hit = any(hit))
}

## ---- tables: gt when present, self-contained fallback otherwise -----------
## Every table is written in THREE forms, because they answer different needs
## and a user should not have to re-run to get the one they wanted: .csv to
## compute on, .html to read, .docx to paste into a manuscript. The docx is the
## one that saves real work - a results table retyped by hand is a table with
## typos in it.
## Display label for a column, falling back to the column name. Cosmetic only.
.lab <- function(frame, x) {
  L <- frame$labels
  if (is.null(L) || !length(L)) return(x)
  vapply(as.character(x), function(z)
    if (!is.na(z) && z %in% names(L) && nzchar(L[[z]][1]))
      unname(L[[z]][1]) else z,
    character(1), USE.NAMES = FALSE)
}




## Add a display label column next to `unit` without touching `unit` itself -
## the raw name is the key a reader joins on, the label is what they read.
.rp_add_unit_label <- function(frame, df) {
  if (is.null(df) || !nrow(df) || !"unit" %in% names(df)) return(df)
  lv <- if ("level" %in% names(df)) df$level else NULL
  lab <- .rp_unit_label(frame, df$unit, lv)
  if (identical(as.character(lab), as.character(df$unit))) return(df)
  j <- match("unit", names(df))
  data.frame(df[, seq_len(j), drop = FALSE], unit_label = lab,
             df[, setdiff(seq_along(df), seq_len(j)), drop = FALSE],
             check.names = FALSE, stringsAsFactors = FALSE)
}

## ---- probe naming --------------------------------------------------------
## A bare "cg25738176" does not say which gene it belongs to, and the gene IS
## the level-local family it was corrected inside - so a reader cannot even
## tell what the correction was over. Render probes as GENE's cgXXXXXXX.
.rp_probe_gene <- function(frame, probe) {
  probe <- as.character(probe)
  mp <- frame$map
  if (is.null(mp) || !all(c("probe", "gene") %in% names(mp)))
    return(rep(NA_character_, length(probe)))
  as.character(mp$gene)[match(probe, as.character(mp$probe))]
}

.rp_unit_label <- function(frame, unit, level = NULL) {
  u <- as.character(unit)
  if (!length(u)) return(u)
  is_probe <- if (!is.null(level)) as.character(level) == "probe"
              else grepl("^cg[0-9]+$", u)
  g <- .rp_probe_gene(frame, u)
  out <- ifelse(is_probe & !is.na(g) & nzchar(g), paste0(g, "'s ", u), u)
  ## an apostrophe is fine in prose and titles but not in a file name
  out
}

## file-name-safe variant: GENE_cgXXXXXXX
.rp_unit_file <- function(frame, unit, level = NULL) {
  u <- as.character(unit)
  is_probe <- if (!is.null(level)) as.character(level) == "probe"
              else grepl("^cg[0-9]+$", u)
  g <- .rp_probe_gene(frame, u)
  ifelse(is_probe & !is.na(g) & nzchar(g), paste0(g, "_", u), u)
}

## ---- text fitting -------------------------------------------------------
## A 202-character moderation title measured 15.25 in into 6.6 in of usable
## width - a 2.31x overflow, centred, so it was cut at BOTH ends. Wrapping by
## CHARACTER count (strwrap(width = 62)) is not a fix either: the width that
## actually fits depends on cex and on the glyphs, so a fixed count overflows
## for wide strings and wastes space for narrow ones. Measure instead.
##
## Returns the wrapped lines and the cex they fit at. Wrapping is tried first
## (it preserves size), then shrinking, down to `min_cex`.
.rp_fit <- function(txt, cex = 0.8, width = NULL, min_cex = 0.45,
                    max_lines = 4L) {
  txt <- paste(as.character(txt), collapse = " ")
  if (!nzchar(txt)) return(list(lines = character(0), cex = cex))
  if (is.null(width)) width <- graphics::par("din")[1] - 0.5
  wrap_at <- function(nl)
    strwrap(txt, width = max(8L, ceiling(nchar(txt) / nl)))
  fits <- function(ls, cx)
    all(is.finite(graphics::strwidth(ls, units = "inches", cex = cx))) &&
    max(graphics::strwidth(ls, units = "inches", cex = cx)) <= width
  for (nl in seq_len(max_lines)) {
    ls <- wrap_at(nl)
    if (fits(ls, cex)) return(list(lines = ls, cex = cex))
  }
  ls <- wrap_at(max_lines); cx <- cex
  while (cx > min_cex) {
    cx <- round(cx - 0.05, 2)
    if (fits(ls, cx)) return(list(lines = ls, cex = cx))
  }
  ## Nothing fits even at min_cex (a single unbreakable token longer than the
  ## device). Truncate visibly rather than let the device clip it silently.
  ## Shrink the BASE string, not the string-with-ellipsis: appending "..." to
  ## (n - 3) characters gives n characters again, so the obvious loop never
  ## terminates. Truncate the base and re-append each time.
  ls <- vapply(ls, function(l) {
    if (graphics::strwidth(l, units = "inches", cex = min_cex) <= width)
      return(l)
    k <- nchar(l)
    while (k > 1L && graphics::strwidth(paste0(substr(l, 1L, k), "..."),
                                        units = "inches", cex = min_cex) > width)
      k <- k - 1L
    paste0(substr(l, 1L, k), "...")
  }, character(1), USE.NAMES = FALSE)
  list(lines = ls, cex = min_cex)
}

## ---- the title block -------------------------------------------------------
## Size the top margin from the text that will go in it, instead of guessing.
##
## The module figure fixed mar[3] at 4.0 lines and drew the first of its
## wrapped title lines at line 4.2. So the moment a title needed three lines,
## line ONE was drawn off the top of the device. The figure still looked
## finished, because what survived - "... level, three lenses, maxT within each
## system" - reads like a title. What went missing was the outcome name, the
## one part a reader cannot reconstruct from the rest of the panel.
##
## Measuring needs a graphics context, and setting the margin has to happen
## before the plot. Both are satisfied by opening the frame with plot.new(),
## measuring, then re-setting par() with new = TRUE so the caller's plot()
## reuses that same page rather than starting a second one.
##
## Only mar[3] moves, and the measurement depends on mar[2] and mar[4] alone,
## so there is no circularity between the two steps.
.rp_head <- function(mar, ttl, sub, ttl_cex = .85, sub_cex = .55, pad = 1.0) {
  graphics::par(mar = mar)
  graphics::plot.new()
  ## mtext(side = 3, adj = 0) starts at the left edge of the plot region and is
  ## free to run on into the right margin, so the usable width is the plot
  ## region PLUS that margin, less a hair so it never touches the device edge.
  w <- graphics::par("pin")[1] + graphics::par("mai")[4] - 0.08
  tf <- .rp_fit(ttl, cex = ttl_cex, width = w)
  sf <- .rp_fit(sub, cex = sub_cex, width = w, max_lines = 6L)
  sub_at <- if (length(sf$lines))
    0.5 + (length(sf$lines) - seq_along(sf$lines)) * 0.75 else numeric(0)
  base <- if (length(sub_at)) max(sub_at) + 1.1 else 0.6
  ttl_at <- base + (length(tf$lines) - seq_along(tf$lines)) * 1.0
  top <- if (length(ttl_at)) max(ttl_at) + pad else mar[3]
  graphics::par(mar = c(mar[1], mar[2], max(mar[3], top), mar[4]), new = TRUE)
  ## The caller draws its plot, then calls draw(). cex comes from .rp_fit, not
  ## from the caller: a title that had to shrink to fit must be DRAWN at the
  ## size it was measured at, or the fitting was for nothing. Both callers used
  ## to take $lines and then hardcode cex = .85, which is exactly that bug.
  function() {
    for (k in seq_along(tf$lines))
      graphics::mtext(tf$lines[k], side = 3, line = ttl_at[k], adj = 0,
                      font = 2, cex = tf$cex)
    for (k in seq_along(sf$lines))
      graphics::mtext(sf$lines[k], side = 3, line = sub_at[k], adj = 0,
                      cex = sf$cex, col = "grey35")
  }
}

## TRUE if every line fits - used by the battery to assert no figure overflows.
.rp_text_fits <- function(txt, cex, width = NULL) {
  if (is.null(width)) width <- graphics::par("din")[1] - 0.5
  ls <- unlist(strsplit(paste(txt, collapse = "\n"), "\n", fixed = TRUE))
  if (!length(ls)) return(TRUE)
  max(graphics::strwidth(ls, units = "inches", cex = cex)) <= width
}

## Draw a fitted title (and optional subtitle) in the outer margin. Returns the
## number of title lines so the caller can size the margin.
.rp_draw_title <- function(main, sub = NULL, cex = 0.8, sub_cex = 0.62,
                           outer = TRUE, top = 1.4) {
  ## strwidth() measures against the CURRENT panel's par("cex"), but mtext(outer
  ## = TRUE) draws at device scale. In a 2x2 layout par("cex") is 0.83, so a
  ## title measured at cex 0.8 came out 1.20x wider than measured and was cut at
  ## both ends even though the fit said it passed. Measure at the size the outer
  ## text will actually have, then convert back when drawing.
  pc <- graphics::par("cex")
  if (!is.finite(pc) || pc <= 0) pc <- 1
  k <- if (isTRUE(outer)) 1 / pc else 1
  ft <- .rp_fit(main, cex = cex * k)
  n <- length(ft$lines)
  for (i in seq_len(n))
    graphics::mtext(ft$lines[i], side = 3, outer = outer,
                    line = top + (n - i) * (0.95 * (ft$cex / k) / 0.8),
                    font = 2, cex = ft$cex / k)
  if (!is.null(sub) && nzchar(paste(sub, collapse = ""))) {
    fs <- .rp_fit(sub, cex = sub_cex * k)
    for (i in seq_along(fs$lines))
      graphics::mtext(fs$lines[i], side = 3, outer = outer,
                      line = top - 0.75 - (i - 1) * 0.7,
                      cex = fs$cex / k, col = "grey35")
  }
  invisible(n)
}

## Report-scoped state. Only use: say the "gt is missing" line once per
## dmsa_report() call rather than once per table.
.rp_env <- new.env(parent = emptyenv())

## ---- spec 53: report QC ledger ---------------------------------------------
## A report component that fails must not become a hole. Every place where a
## component used to be swallowed by tryCatch(..., error = function(e) NULL)
## now goes through .rp_try(): the failure is RECORDED (component, where,
## what it fell back to, the message), the analysis continues exactly as
## before, and at the end the ledger is written to tables/qc_report.csv,
## listed in summary.md and counted once on the console. A user can then
## tell "absent because not applicable" from "absent because it broke".
.rp_qc <- new.env(parent = emptyenv())
.rp_qc_reset <- function() assign("rows", list(), envir = .rp_qc)
.rp_qc_note <- function(component, message, outcome = NA_character_,
                        level = NA_character_, unit = NA_character_,
                        fallback = NA_character_, severity = "failed") {
  rows <- if (exists("rows", envir = .rp_qc, inherits = FALSE))
    get("rows", envir = .rp_qc) else list()
  rows[[length(rows) + 1L]] <- data.frame(
    component = as.character(component), severity = as.character(severity),
    outcome = as.character(outcome), level = as.character(level),
    unit = as.character(unit), fallback = as.character(fallback),
    message = gsub("[\r\n]+", " ", as.character(message)),
    stringsAsFactors = FALSE)
  assign("rows", rows, envir = .rp_qc)
  invisible(NULL)
}
.rp_qc_rows <- function() {
  rows <- if (exists("rows", envir = .rp_qc, inherits = FALSE))
    get("rows", envir = .rp_qc) else list()
  if (!length(rows))
    return(data.frame(component = character(0), severity = character(0),
                      outcome = character(0), level = character(0),
                      unit = character(0), fallback = character(0),
                      message = character(0), stringsAsFactors = FALSE))
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}
## evaluate `expr`; on error record it in the ledger and return NULL
.rp_try <- function(component, expr, outcome = NA_character_,
                    level = NA_character_, unit = NA_character_,
                    fallback = NA_character_) {
  tryCatch(expr, error = function(e) {
    .rp_qc_note(component, conditionMessage(e), outcome = outcome,
                level = level, unit = unit, fallback = fallback)
    NULL
  })
}

.rp_write_table <- function(df, file, frame, title = "") {
  utils::write.csv(df, paste0(file, ".csv"), row.names = FALSE)
  ok_gt <- requireNamespace("gt", quietly = TRUE)
  ext <- frame$table_type
  wants <- unique(c("html", if (ext %in% c("docx", "rtf")) ext else NULL,
                    "docx"))
  if (ok_gt) {
    tb <- gt::gt(df)
    if (nzchar(title)) tb <- gt::tab_header(tb, title = title)
    for (w in wants)
      tryCatch(gt::gtsave(tb, paste0(file, ".", w)), error = function(e) {
        ## html has a dependency-free fallback; the others do not, so a failure
        ## there has to be said out loud rather than leaving a missing file the
        ## user only notices when pasting into a manuscript.
        if (w == "html") .rp_html_fallback(df, paste0(file, ".html"), title)
        else message("could not write .", w, " for ", basename(file), ": ",
                     conditionMessage(e), " - the .csv and .html are still there.")
      })
  } else {
    ## No gt. Say so unconditionally: .docx is now written for every run, not
    ## only when table_type asks for it, so staying quiet here leaves the user
    ## looking for a Word table that was never attempted.
    if (!isTRUE(.rp_env$gt_said)) {
      .rp_env$gt_said <- TRUE
      no_ext <- if (!ext %in% c("html", "docx")) paste0(" and no .", ext)
      else ""
      message("gt is not installed, so the tables are written as .csv and ",
              ".html only - no .docx", no_ext,
              ". install.packages(\"gt\") adds the Word table.")
    }
    .rp_html_fallback(df, paste0(file, ".html"), title)
  }
  invisible(file)
}

.rp_html_fallback <- function(df, out, title = "") {
  esc <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", as.character(x)))
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(z) ifelse(is.finite(z),
                                                sprintf("%.4g", z), ""))
  rows <- apply(df, 1, function(r)
    paste0("<tr><td>", paste(esc(r), collapse = "</td><td>"), "</td></tr>"))
  html <- c("<!DOCTYPE html><html><head><meta charset='utf-8'><style>",
    "body{font-family:system-ui,sans-serif;margin:24px;color:#222}",
    "h3{font-weight:600}table{border-collapse:collapse;font-size:13px}",
    "th{border-bottom:2px solid #444;text-align:left;padding:4px 10px}",
    "td{border-bottom:1px solid #ddd;padding:4px 10px;font-variant-numeric:tabular-nums}",
    "tr:hover{background:#f6f6f6}</style></head><body>",
    if (nzchar(title)) paste0("<h3>", esc(title), "</h3>") else "",
    "<table><tr><th>", paste(esc(names(df)), collapse = "</th><th>"),
    "</th></tr>", rows, "</table></body></html>")
  writeLines(html, out)
  invisible(out)
}

## ---- the module figure per outcome ----------------------------------------
## The module level had no figure at all, which left the one thing DMSA can say
## that a gene-level engine cannot - WHERE IN THE SYSTEM the signal sits -
## visible only in a CSV. It also strands the module evidence audit: a module
## result is only as good as the module's definition, and the tier belongs
## beside the p, not in an appendix a reader will not open.
.rp_fig_modules <- function(frame, res, oc, file) {
  pal <- .rp_pal(frame); COH <- pal[1]; CMP <- pal[3]; DIF <- pal[4]
  r <- res[res$outcome == oc & res$level == "module" & res$n_probes > 0, ,
           drop = FALSE]
  if (!nrow(r)) return(invisible(NULL))
  r <- r[order(r$system_id, r$p_omnibus), , drop = FALSE]

  ## evidence tier per module, if the run carries the audit
  ev <- frame$module_evidence
  tier <- rep(NA_character_, nrow(r))
  if (!is.null(ev) && all(c("module", "evidence_strength") %in% names(ev)))
    tier <- ev$evidence_strength[match(r$unit, ev$module)]

  fp <- .rp_dev(frame, file, width = 8.4,
                height = max(3.4, 0.30 * nrow(r) + 1.9))
  on.exit(grDevices::dev.off(), add = TRUE)
  n <- nrow(r); yy <- rev(seq_len(n))
  ## Title and subtitle are measured first and the top margin sized to hold
  ## them; head() draws them once the panel is up.
  ttl <- sprintf("%s - module level, three lenses, %s within each system",
                 .lab(frame, oc), frame$correction)
  sub <- paste("Each block is one family: a module is corrected only against",
               "the other modules of its own system. Right margin: evidence",
               "tier for the module's DEFINITION, then probe count.",
               "Bold = survives.")
  head <- .rp_head(c(4.2, 17, 4.0, 6.5), ttl, sub)
  graphics::plot(NA, xlim = c(0.002, 1.3), ylim = c(.4, n + .6), log = "x",
                 axes = FALSE, xlab = "", ylab = "")
  at <- c(.005, .01, .02, .05, .1, .2, .5, 1)
  graphics::axis(1, at = at, labels = format(at), cex.axis = .7,
                 col = "grey55", col.axis = "grey30")
  graphics::abline(v = .05, lty = 2, col = "grey55")
  ## module labels wrap: they are sentences, not symbols
  lab <- vapply(as.character(r$unit), function(z)
    paste(strwrap(z, width = 34), collapse = "\n"), character(1))
  graphics::mtext(lab, side = 2, at = yy, las = 1, line = .4, cex = .52,
                  font = ifelse(r$selected, 2, 1))
  ## The family boundary has to be visible, because it is the claim the panel
  ## makes: each module was corrected against its OWN system's modules, not
  ## against all ten rows on the page. Name each block and rule between them.
  sysv <- if (!is.null(r$system)) as.character(r$system) else
    as.character(r$system_id)
  grp <- split(seq_len(n), factor(sysv, levels = unique(sysv)))
  for (k in seq_along(grp)) {
    idx <- grp[[k]]
    graphics::mtext(paste(strwrap(names(grp)[k], width = 18), collapse = "\n"),
                    side = 2, at = mean(yy[idx]), las = 1, line = 11.6,
                    cex = .55, font = 2, col = "grey30")
    if (k < length(grp))
      graphics::abline(h = min(yy[idx]) - .5, col = "grey85", lwd = .8)
  }
  graphics::points(r$p_coherence_adj, yy + .18, pch = 19, cex = .7, col = COH)
  graphics::points(r$p_composite_adj, yy, pch = 19, cex = .7, col = CMP)
  graphics::points(r$p_diffuse_adj, yy - .18, pch = 19, cex = .7, col = DIF)
  ## right margin: the audit tier, and how much the module rests on
  graphics::mtext(sprintf("%s%d %s",
                          ifelse(is.na(tier), "", paste0(tier, ", ")),
                          r$n_probes, ifelse(r$n_probes == 1, "probe", "probes")),
                  side = 4, at = yy, las = 1, line = .4, cex = .5,
                  col = ifelse(is.na(tier) | tier == "High", "grey40",
                               "firebrick"))
  head()
  graphics::mtext("family-adjusted p", side = 1, line = 2.4, cex = .7)
  graphics::legend("bottomleft", bty = "n", cex = .6, pch = 19,
                   col = c(COH, CMP, DIF),
                   legend = c("coherence", "composite", "diffuse"))
  invisible(fp)
}

## ---- the overview figure per outcome --------------------------------------
## Spec 29-38 (2026-08-29): the overview obeys FIGURE RULES instead of
## growing without bound. One panel PER SYSTEM, each in that system's own
## accent colour; at most 40 units per panel, further units paginated into
## _p2, _p3, ... files; units whose calibrated q (re-ruled 2026-10-15) falls below
## twice the FDR level get a light background band so near-threshold rows catch the eye
## without claiming significance. Returns every file written.
.rp_fig_overview <- function(frame, res, oc, file) {
  r0 <- res[res$outcome == oc & res$n_probes > 0, , drop = FALSE]
  if (!nrow(r0)) return(invisible(NULL))
  pal <- .rp_pal(frame); COH <- pal[1]; CMP <- pal[3]; DIF <- pal[4]
  sids <- unique(r0$system_id)
  syscol <- grDevices::hcl.colors(max(3L, length(sids)), "Dark 3")
  files <- character(0)
  .slug <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
    substr(gsub("^_|_$", "", x), 1, 28)
  }
  ## everything in this panel is keyed to ONE statistic - the best lens's
  ## family-adjusted p, which is the naming rule (the dots ARE the per-lens
  ## adjusted p's, so a dot left of the alpha line and a bold label and a
  ## shaded row now all say the same thing; PI, 2026-08-29). Re-ruled
  ## 2026-10-15: the naming statistic is the dependence-calibrated q, drawn
  ## as a black diamond against the FDR line; order, bold and shading follow
  ## it. The three coloured dots stay as the pre-registered any-lens view.
  r0$.lens_adj <- pmin(r0$p_coherence_adj, r0$p_composite_adj,
                       r0$p_diffuse_adj, na.rm = TRUE)
  .fdr <- frame$fdr %||% .10
  if (is.null(r0$q_unit_dep)) r0$q_unit_dep <- NA_real_
  for (si in seq_along(sids)) {
    rs <- r0[r0$system_id == sids[si], , drop = FALSE]
    rs <- rs[order(rs$q_unit_dep, rs$p_unit, rs$.lens_adj, rs$p_omnibus),
             , drop = FALSE]
    pages <- split(seq_len(nrow(rs)), ceiling(seq_len(nrow(rs)) / 40))
    for (pg in seq_along(pages)) {
      r <- rs[pages[[pg]], , drop = FALSE]
      sysnm <- r$system[1]
      fbase <- paste0(file,
                      if (length(sids) > 1L) paste0("_", .slug(sysnm)) else "",
                      if (length(pages) > 1L) paste0("_p", pg) else "")
      fp <- .rp_dev(frame, fbase, width = 7.2,
                    height = max(3.2, 0.24 * nrow(r) + 1.9))
      n <- nrow(r); yy <- rev(seq_len(n))
      ttl <- sprintf("%s - %s%s", .lab(frame, oc), sysnm,
                     if (length(pages) > 1L)
                       sprintf(" (%d of %d)", pg, length(pages)) else "")
      ## no explanatory subtitle (PI, 2026-08-29): the visual codes are
      ## defined once in summary.md; the figure carries only its title
      head <- .rp_head(c(4.2, 11, 3.4, 5.5), ttl, "", sub_cex = .6)
      graphics::plot(NA, xlim = c(0.002, 1.3), ylim = c(.4, n + .6),
                     log = "x", yaxt = "n", xaxt = "n", xlab = "", ylab = "",
                     bty = "n")
      ## near-threshold shading BEFORE the points so it sits underneath -
      ## keyed to the SAME statistic as bold and the dots (the naming rule)
      near <- which(is.finite(r$q_unit_dep) & r$q_unit_dep < 2 * .fdr)
      for (k in near)
        graphics::rect(0.0015, yy[k] - .42, 1.45, yy[k] + .42,
                       col = grDevices::adjustcolor(syscol[si], .10),
                       border = NA)
      graphics::abline(v = frame$alpha, col = "grey55", lty = 2)
      graphics::abline(v = .fdr, col = "black", lty = 3)
      at <- c(.005, .01, .02, .05, .1, .2, .5, 1)
      graphics::axis(1, at = at, labels = at, cex.axis = .62, col = "grey60",
                     col.axis = "grey30")
      .ec <- if (!is.null(r$exact_confirmed)) r$exact_confirmed %in% TRUE else
        rep(FALSE, nrow(r))
      .oc2 <- if (!is.null(r$omnibus_confirmed))
        r$omnibus_confirmed %in% TRUE else rep(FALSE, nrow(r))
      graphics::mtext(paste0(as.character(r$unit), ifelse(.ec, " \u2020", ""),
                             ifelse(.oc2, " \u2021", "")),
                      side = 2, at = yy, las = 1,
                      line = .4, cex = .6, font = ifelse(r$selected, 2, 1),
                      col = syscol[si])
      graphics::points(r$p_coherence_adj, yy + .16, pch = 19, cex = .7,
                       col = COH)
      graphics::points(r$p_composite_adj, yy, pch = 19, cex = .7, col = CMP)
      graphics::points(r$p_diffuse_adj, yy - .16, pch = 19, cex = .7,
                       col = DIF)
      graphics::points(r$q_unit_dep, yy, pch = 18, cex = .95, col = "black")
      head()
      graphics::mtext("calibrated q / family-adjusted p", side = 1,
                      line = 2.4, cex = .7)
      graphics::mtext(sysnm, side = 3, line = -0.1, adj = 1, cex = .6,
                      col = syscol[si], font = 2)
      graphics::legend("bottomleft", bty = "n", cex = .6, pch = c(18, 19, 19, 19),
                       col = c("black", COH, CMP, DIF),
                       legend = c("calibrated q", "coherence",
                                  "composite", "diffuse"))
      grDevices::dev.off()
      files <- c(files, fp)
    }
  }
  invisible(files)
}

## locus panel per surviving gene
.rp_fig_locus <- function(frame, oc, sid, gene, dir_ok, file, fits = NULL) {
  s <- frame$sets[[as.character(sid)]]
  gmap <- s$map[s$map$gene == gene, , drop = FALSE]
  if (!nrow(gmap)) return(invisible(NULL))
  ## spec 52: ONE estimator for the figure and the test. The gene test pooled
  ## per-probe fits of the WINSORISED, STANDARDISED, covariate-adjusted (and,
  ## when a random intercept was declared, RI-transformed) probe on the term;
  ## those b/se are what the panel draws when the engine handed them over
  ## (attr "probe_fits" of the gene-level triangulation). The nominal p that
  ## greys a probe comes from that same fit. The report-side OLS on raw
  ## M-values is now only the fallback for a panel drawn without a test.
  .ix <- if (!is.null(fits)) match(gmap$column, fits$column) else NA_integer_
  if (!is.null(fits) && !anyNA(.ix)) {
    .dfr <- attr(fits, "dfr") %||% (nrow(frame$data) - 2L)
    b <- fits$b[.ix]; se <- fits$se[.ix]
    pr <- data.frame(probe = gmap$probe, b = b, se = se,
                     p = 2 * stats::pt(-abs(b / se), df = .dfr),
                     d = gmap$best_direction, stringsAsFactors = FALSE)
    .est <- "engine"
  } else {
    f <- .rp_probe_fits(frame, oc, gmap$column)
    pr <- data.frame(probe = gmap$probe, b = f$b, se = f$se, p = f$p,
                     d = gmap$best_direction, stringsAsFactors = FALSE)
    .est <- "ols_fallback"
  }
  .rows <- data.frame(outcome = oc, system_id = s$system_id %||% as.character(sid),
                      gene = gene,
                      probe = pr$probe, column = gmap$column, d = pr$d,
                      b = pr$b, se = pr$se, z = pr$b / pr$se, p_nominal = pr$p,
                      estimator = .est,
                      scale = if (.est == "engine")
                        "winsorised + standardised M per SD of the term (the test's own per-probe fit)" else
                        "raw M (report OLS; no test ran for this panel)",
                      stringsAsFactors = FALSE)
  ## Gene-level significance CAN come with no individually significant CpG -
  ## the lenses pool small direction-consistent shifts, so evidence
  ## accumulates at the gene even when every CpG sits below the noise line.
  ## That is counter-intuitive on a figure whose probes are all grey, so when
  ## it happens the panel says so under its own title instead of leaving the
  ## reader to reconcile "significant gene" with "no significant probe".
  .none <- !any(is.finite(pr$p) & pr$p < frame$alpha)
  .ctx <- paste0(.lab(frame, oc),
                 if (.est == "engine")
                   " - bars: the gene test's own per-probe fit (winsorised, standardised M)"
                 else " - bars: report-side OLS on raw M (no test ran here)")
  if (.none)
    .ctx <- sprintf(paste0("%s - no single CpG clears p < %.2g on its own ",
                           "(grey): the gene-level result pools small, ",
                           "direction-consistent shifts across all %d CpGs"),
                    .ctx, frame$alpha, nrow(pr))
  ## Coordinates come from the CASCADE, looked up by probe id, because it is
  ## the same table the gene model is aligned to (hg38) and it is always
  ## present. dmsa_probe_coords() is the fallback only: it reads whatever
  ## manifest the companion package ships, which for EPIC is hg19, and mixing
  ## an hg19 probe with an hg38 exon puts a CpG ~19 kb from its own gene
  ## without either layer being wrong on its own.
  pos <- NULL
  cas <- .rp_try("locus: cascade coordinates", dmsa_sets(if (is.null(frame$sets_source)) "alpha"
                            else frame$sets_source)$cascade,
                 outcome = oc, level = "gene", unit = gene,
                 fallback = "manifest coordinates, else evenly spaced probes")
  if (!is.null(cas) && all(c("cpg", "pos_hg38") %in% names(cas))) {
    i <- match(gmap$probe, cas$cpg)
    if (sum(!is.na(i)) == length(i)) {
      pr$chr <- cas$chr[i]
      pr$pos <- suppressWarnings(as.numeric(cas$pos_hg38[i]))
      pr$genome <- "hg38"
      pos <- "cascade"
    }
  }
  if (is.null(pos)) {
    co <- .rp_try("locus: manifest coordinates", dmsa_probe_coords(gmap$probe),
                  outcome = oc, level = "gene", unit = gene,
                  fallback = "probes evenly spaced")
    if (!is.null(co) && all(c("probe", "pos") %in% names(co)))
      pr <- merge(pr, co[, intersect(c("probe", "chrom", "pos"), names(co))],
                  by = "probe", all.x = TRUE)
  }
  ## `gene_models` on the frame decides whether the report goes to the network
  ## for exon structure. A report must not silently make HTTP calls, so this is
  ## opt-in, and a failed fetch degrades to the coordinate axis rather than to
  ## an invented one.
  ## under "auto", test runs and R CMD check never touch the network
  ## (Bioconductor policy; a fake test gene would fail the lookup anyway) -
  ## an explicit gene_models = TRUE still fetches even there
  .no_net <- identical(frame$gene_models, "auto") &&
    (identical(Sys.getenv("TESTTHAT"), "true") ||
     nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")))
  gm <- if (!.no_net && (isTRUE(frame$gene_models) ||
            identical(frame$gene_models, "auto"))) {
    ## "auto" (the default): the exon/TSS model is fetched for NAMED genes
    ## only - this function is only ever called for them - and a failed
    ## fetch degrades to the bare coordinate axis with a message instead of
    ## quietly dropping the exon view (PI, 2026-08-29)
    .gerr <- NULL
    .g <- tryCatch(dmsa_gene_model(gene, quiet = TRUE),
                   error = function(e) { .gerr <<- conditionMessage(e); NULL })
    if (is.null(.g) || (is.data.frame(.g) && !nrow(.g))) {
      ## spec 53: a missing exon model is a recorded fallback, not a hole
      .rp_qc_note("locus: gene model (Ensembl)",
                  .gerr %||% "fetch returned no model", outcome = oc,
                  level = "gene", unit = gene, severity = "fell back",
                  fallback = "coordinate axis without exons")
      message("gene model for ", gene, " could not be fetched from Ensembl ",
              "(offline?); the locus panel keeps true coordinates but ",
              "cannot draw the exons. dmsa_gene_model() built once and ",
              "passed as gene_models = <table> works offline.")
    }
    .g
  } else if (.no_net || isFALSE(frame$gene_models) ||
             identical(frame$gene_models, "auto")) NULL else
    frame$gene_models          # NULL, or a table the caller built once
  if (is.data.frame(gm) && "gene" %in% names(gm))
    gm <- gm[gm$gene == gene, , drop = FALSE]
  ## A NAMED gene's locus panel MUST be printed (PI, 2026-08-29): the full
  ## panel (gene model, genomic coordinates) is attempted first; if any of
  ## that layered context fails to draw, the panel falls back to the bare,
  ## always-drawable form - probes evenly spaced, no model - rather than
  ## silently disappearing from the report. Only if even the bare panel
  ## errors does this return NULL, and then the caller records the failure
  ## in summary.md instead of leaving a hole where a figure should be.
  ok <- tryCatch({
    dmsa_plot_locus(pr, gene = gene, gene_model = gm,
                    signal_p = frame$alpha, context = .ctx,
                    file = paste0(file, ".", frame$plot_type))
    TRUE
  }, error = function(e) {
    message("full locus panel for ", gene, " failed (",
            conditionMessage(e), "); drawing the bare panel instead")
    .rp_qc_note("figure: locus panel (full)", conditionMessage(e), outcome = oc,
                level = "gene", unit = gene, severity = "fell back",
                fallback = "bare panel drawn (probes evenly spaced, no model)")
    FALSE
  })
  if (!ok)
    ok <- tryCatch({
      suppressMessages(dmsa_plot_locus(
        pr[, intersect(c("probe", "b", "se", "p", "d"), names(pr)),
           drop = FALSE],
        gene = gene, gene_model = NULL, signal_p = frame$alpha,
        context = .ctx, file = paste0(file, ".", frame$plot_type)))
      TRUE
    }, error = function(e) {
      message("locus panel for ", gene, " could not be drawn at all: ",
              conditionMessage(e))
      .rp_qc_note("figure: locus panel", conditionMessage(e), outcome = oc,
                  level = "gene", unit = gene, fallback = "no locus figure")
      FALSE
    })
  if (ok) invisible(structure(file, probe_rows = .rows)) else invisible(NULL)
}

## ---- progress and completion signal ---------------------------------------
## A DMSA report is minutes of permutation with nothing on screen: the parent
## cohort at B = 1999 across 19 systems runs 35 minutes. A bar is not a nicety
## at that length, it is the difference between "working" and "hung".
##
## Both are OFF outside an interactive session by default. A progress bar
## writes carriage returns, which turn a piped log into unreadable soup, and a
## package that makes NOISE during R CMD check or a testthat run is a package
## people learn to mute. interactive() is TRUE at the console and FALSE under
## Rscript, check and testthat, which is exactly the line wanted.
.rp_bar <- function(total, on = TRUE) {
  noop <- list(tick = function(n = 1L) invisible(NULL),
               msg  = function(...) invisible(NULL),
               done = function() invisible(NULL))
  if (!isTRUE(on) || !is.finite(total) || total < 1) return(noop)
  pb <- utils::txtProgressBar(min = 0, max = total, style = 3, width = 42)
  i <- 0L; shut <- FALSE
  list(
    tick = function(n = 1L) {
      i <<- min(i + n, total); utils::setTxtProgressBar(pb, i); invisible(NULL)
    },
    ## print a line without leaving the bar half-drawn above it
    msg = function(...) {
      cat("\n"); message(...); utils::setTxtProgressBar(pb, i); invisible(NULL)
    },
    ## idempotent: the run closes the bar before its final message, and the
    ## on.exit handler closes it again on any exit path. Closing twice would
    ## print a stray newline into the middle of the summary.
    done = function() {
      if (shut) return(invisible(NULL))
      shut <<- TRUE
      utils::setTxtProgressBar(pb, total); close(pb); invisible(NULL)
    })
}

## beepr is a Suggests: absent, this is silently a no-op, which is the correct
## behaviour for a convenience that must never be load-bearing.
.rp_beep <- function(beep) {
  if (is.null(beep) || isFALSE(beep)) return(invisible(NULL))
  snd <- if (isTRUE(beep)) 8 else beep
  if (requireNamespace("beepr", quietly = TRUE))
    try(beepr::beep(sound = snd), silent = TRUE)
  invisible(NULL)
}

## how many triangulation calls the declared design will make - the unit of
## work that actually costs time, so the bar is proportional rather than a
## decorative sweep
.rp_work <- function(frame, levels_on) {
  n_sys <- length(frame$sets)
  n_mod <- if (!is.null(frame$modules))
    length(intersect(unique(frame$modules$system_id),
                     as.numeric(names(frame$sets)))) else 0L
  per_outcome <- ("system" %in% levels_on) +
    ("gene" %in% levels_on) * n_sys + ("module" %in% levels_on) * n_mod
  max(1L, as.integer(length(frame$outcome) * per_outcome))
}

#' @rdname dmsa_report
#' @details For a \code{dmsa_frame} this is call 2 of the two-call interface:
#'   it runs each requested level x outcome with the locked statistics,
#'   applies correction ONLY within level-local families, runs the always-on
#'   brain-bridge post-hoc on surviving probes, and writes everything under
#'   \code{frame$outdir}.
#' @param progress Show a 0-100% progress bar. \code{NULL} (default) takes the
#'   value declared on the frame, which itself defaults to \code{interactive()}.
#' @param beep Completion sound. \code{NULL} (default) takes the frame's value.
#'   \code{TRUE} is \pkg{beepr} sound 8; a number picks another; \code{FALSE}
#'   is silent. Without \pkg{beepr} installed this is a no-op.
#' @param overwrite Replace an existing report in `outdir`? Default FALSE:
#'   a directory that already holds DMSA output is refused rather than
#'   silently overwritten.
#' @examples
#' set.seed(1)
#' map <- data.frame(gene = rep(c("NR3C1", "FKBP5"), each = 2), system_id = 1L,
#'                   system = "HPA axis", probe = paste0("cg0", 1:4),
#'                   column = paste0("cg0", 1:4), best_direction = c(-1, 1, -1, 1),
#'                   p_plus = c(.1, .9, .1, .9))
#' dat <- data.frame(anx = rnorm(40), cov1 = rnorm(40), cID = rep(1:20, each = 2))
#' ## methylation tracks anxiety in each probe's own expression direction
#' for (i in 1:4)
#'   dat[[map$column[i]]] <- plogis(rnorm(40) + .6 * dat$anx * map$best_direction[i])
#'
#' ## everything lands under outdir; figures and tables are switched off here
#' ## only to keep the example quick
#' fr <- dmsa_frame(dat, map = map, outcome = "anx", covariates = "cov1",
#'                  B = 49, seed = 1, outdir = tempfile(),
#'                  plots = FALSE, tables = FALSE)
#' r <- dmsa_report(fr)
#' r$results[, c("level", "unit", "n_probes", "p_omnibus")]
#' @rdname dmsa_report
#' @export
dmsa_report.dmsa_frame <- function(frame, progress = NULL, beep = NULL,
                                   overwrite = FALSE, ...) {
  .rp_env$gt_said <- FALSE
  stopifnot(inherits(frame, "dmsa_frame"))
  t_start <- Sys.time()
  ## report argument wins, then the frame's declared value, then the default
  if (is.null(progress)) progress <- frame$progress %||% interactive()
  if (is.null(beep))     beep     <- frame$beep     %||% interactive()
  outdir <- frame$outdir
  ## spec 54: never silently replace an earlier analysis. A DMSA report is a
  ## scientific artefact; overwriting one without saying so destroys the record
  ## a result was read from. Presence is judged by DMSA's own output shape, so
  ## an unrelated directory the user happens to point at is not misread.
  if (!isTRUE(overwrite)) {
    .prev <- c(list.files(file.path(outdir, "figures"), pattern = "\\.(png|pdf)$"),
               list.files(file.path(outdir, "tables"), pattern = "\\.(csv|html|docx|rtf)$"),
               ## a plots=FALSE, tables=FALSE run leaves ONLY summary.md;
               ## it is no less an analysis record than the figures are
               if (file.exists(file.path(outdir, "summary.md")))
                 "summary.md")
    if (length(.prev))
      stop("`outdir` already contains a DMSA report: ",
           normalizePath(outdir, winslash = "/", mustWork = FALSE), "\n",
           length(.prev), " existing output file(s), e.g. ",
           paste(utils::head(basename(.prev), 3), collapse = ", "), "\n",
           "Refusing to replace an earlier analysis. Point `outdir` at a new ",
           "directory, or call dmsa_report(frame, overwrite = TRUE) if you ",
           "meant to discard it.\nNo report was written.", call. = FALSE)
  }
  dir.create(file.path(outdir, "figures"), recursive = TRUE,
             showWarnings = FALSE)
  dir.create(file.path(outdir, "tables"), recursive = TRUE,
             showWarnings = FALSE)
  levels_on <- names(frame$levels)[frame$levels]
  ## stale-frame detection: same package version, DIFFERENT installation
  .cur_built <- tryCatch(
    utils::packageDescription("dmsa")[["Built"]] %||% NA_character_,
    error = function(e) NA_character_)
  if (!is.null(frame$built_stamp) && !is.na(frame$built_stamp) &&
      !is.na(.cur_built) && !identical(frame$built_stamp, .cur_built))
    message("NOTE: this frame was built by a DIFFERENT dmsa installation ",
            "than the one now loaded. Stored fields (outcome families, ",
            "labels, maps) may predate recent fixes - rerun dmsa_frame() ",
            "to rebuild the frame before trusting this report.")
  ## moderation requires gaussian families - checked for EVERY outcome before
  ## any permutation work is spent, not mid-loop after outcome 1's battery ran
  if (frame$moderation && identical(frame$frame_role, "predictor")) {
    ## the gaussian requirement applies ONLY when the OUTCOME is the
    ## response of the moderated model (frame_role = "predictor": oc ~ S x
    ## mod). Under frame_role = "outcome" the model is S ~ oc x mod - the
    ## TONE SCORE is the response, and a two-level oc on the right-hand
    ## side is an ordinary group contrast, so refusing it was wrong (PI's
    ## pills x sex battery, 2026-08-29).
    .bad <- Filter(function(oc) .rp_otype(frame, oc) != "gaussian",
                   frame$outcome)
    if (length(.bad))
      stop("moderation with frame_role = \"predictor\" requires ",
           "outcome_type = 'gaussian' (the moderated composite model puts ",
           "the outcome on the left and is linear); outcome(s) ",
           paste(.bad, collapse = ", "), " declare another family.\n",
           "No report was written.", call. = FALSE)
  }
  all_res <- list(); mod_res <- list(); shape_rows <- list()
  ## spec 52: the per-probe fits each gene test was built from (one
  ## data.frame per outcome x system), so the locus panel draws the
  ## test's own numbers; and the rows every locus panel actually drew
  probe_fits <- list(); locus_rows <- list()
  .rp_qc_reset()   # spec 53: fresh ledger per run
  sys_skip <- FALSE   # set when the system level is skipped for missing polarity
  bar <- .rp_bar(.rp_work(frame, levels_on) + 2L, on = progress)
  on.exit({ try(bar$done(), silent = TRUE); .rp_beep(beep) }, add = TRUE)

  for (oc in frame$outcome) {
    .oty <- .rp_otype(frame, oc)   # spec 50: this outcome's declared family
    if (frame$moderation) {
      ## composite-only, per spec: interactions ride the composite lens
      ## (families validated for every outcome before the loop).
      ## MODERATION NO LONGER REPLACES THE MAIN ANALYSIS (PI, 2026-08-29:
      ## a moderation run produced no system/module/gene results or
      ## figures at all): the moderated battery runs IN ADDITION to the
      ## main triangulation, so one report carries both. The price is
      ## runtime - both batteries permute.
      for (lv in intersect(levels_on, c("system", "gene", "module", "probe"))) {
        units_list <- .rp_units(frame, lv)
        for (fam in names(units_list)) {
          mo <- .report_moderated(frame, oc, units_list[[fam]], fam)
          mo$level <- lv; mo$outcome <- oc
          mod_res[[paste(oc, lv, fam)]] <- mo
        }
      }
    }
    ## The shape scan runs at every declared level, not only at the system.
    ## An aligned tone score exists for a module and for a gene too, so there
    ## was never a reason to ask about functional form at one level only.
    ## The shape scan puts the OUTCOME on the left (oc ~ S + covariates), unlike
    ## the main-effect path which puts methylation there. So a non-gaussian
    ## outcome is a real problem here and not in the main path: the arms would
    ## be a linear probability model, and the Sasabuchi/Fieller machinery on top
    ## of it would be meaningless. Moderation already refuses this; the scan
    ## must too rather than quietly producing shape statistics for a 0/1.
    if (frame$type != "linear" && .oty != "gaussian") {
      ## The Sasabuchi/Fieller arms genuinely cannot run here - they assume a
      ## least-squares fit, and a 0/1 response fitted that way is a linear
      ## probability model. But a binary outcome IS non-linear, and the response
      ## curve is well defined on the logit scale, so a figure is drawn and the
      ## curvature is tested there rather than nothing being reported at all.
      message("outcome_type = '", .oty, "' for outcome '", oc,
              "': the Sasabuchi and ",
              "Fieller arms need a least-squares fit and are not run. For a ",
              "surviving unit of a two-level outcome, a logistic response ",
              "curve is drawn instead, with curvature tested ",
              "on the logit scale (LRT against the linear logistic model). The ",
              "unit-level DMSA results are unaffected - there the outcome is a ",
              "predictor of methylation, not a response.")
    } else if (frame$type != "linear") {
      for (lv in intersect(levels_on, c("system", "module", "gene", "probe"))) {
        ul <- .rp_try("shape scan: units", .rp_units(frame, lv), outcome = oc,
                      level = lv, fallback = "level absent from the shape scan")
        if (is.null(ul)) next
        for (fam in names(ul)) {
          sr <- .rp_try("shape scan: family", .report_shape_family(frame, oc, ul[[fam]], fam, lv),
                        outcome = oc, level = lv, unit = fam,
                        fallback = "family absent from the shape scan")
          if (!is.null(sr)) shape_rows[[paste(oc, lv, fam)]] <- sr
        }
      }
    }
    ## ---- system level: one family = the named systems --------------------
    if ("system" %in% levels_on && !isTRUE(sys_skip)) {
      sy <- .rp_system_frame(frame)
      if (inherits(sy, "dmsa_no_polarity")) {
        message("system level SKIPPED: ", attr(sy, "why"), ".\n",
                "  A system score needs each gene's sign against the system's ",
                "activation tone.\n  Weighting every gene +1 would report a ",
                "purely activating system whatever the biology says, so it is ",
                "not done.\n  Gene and probe levels are unaffected. Supply ",
                "polarity with your reference to enable the system level.")
        ## flag, don't mutate levels_on: the shape scan and moderation for
        ## LATER outcomes read levels_on before this block, so removing
        ## "system" here made outcome 1 and outcome 2 scan different levels
        sys_skip <- TRUE
        bar$tick()
      } else {
      ## spec 44: the declared chip reaches EVERY level. `.report_gene_level()`
      ## has always passed it, so gene, module and probe were fitted with the
      ## chip random intercept while the system level - the headline - was
      ## fitted without it, from the same frame and the same declaration. The
      ## argument was simply not forwarded here.
      r <- dmsa_triangulate(sy$M, frame$data, .rp_rhs(frame, oc), oc,
                            sy$units, sy$alignment, block = frame$block,
                            B = frame$B, correction = frame$correction,
                            weighting = frame$weighting %||% "combined",
                            w_floor = frame$w_floor %||% 1.5, seed = frame$seed,
                            ri_group = if (.nz(frame$chip_random))
                                         frame$data[[frame$chip_random]]
                                       else NULL)
      r$system_id <- NA; r$system <- r$unit; r$outcome <- oc; r$level <- "system"
      .mb <- attr(r, "union_null_min")
      r$fwer_realized <- if (length(.mb) && any(is.finite(.mb)))
        mean(.mb < frame$alpha, na.rm = TRUE) else NA_real_
      .qm <- attr(r, "unit_null_qmin")
      r$null_alarm_bh <- if (length(.qm) && any(is.finite(.qm)))
        mean(.qm <= (frame$fdr %||% .10), na.rm = TRUE) else NA_real_
      all_res[[paste(oc, "system")]] <- r
      bar$tick()
      }
    }
    ## ---- gene level: family = the system's genes --------------------------
    if ("gene" %in% levels_on) {
      for (sid in names(frame$sets)) {
        r <- .report_gene_level(frame, oc, sid)
        probe_fits[[paste(oc, sid)]] <- attr(r, "probe_fits")
        r$level <- "gene"
        all_res[[paste(oc, "gene", sid)]] <- r
        bar$tick()
      }
    }
    ## ---- module level: family = the system's modules ----------------------
    if ("module" %in% levels_on && !is.null(frame$modules)) {
      for (sid in intersect(unique(frame$modules$system_id),
                            as.numeric(names(frame$sets)))) {
        s <- frame$sets[[as.character(sid)]]
        mm <- frame$modules[frame$modules$system_id == sid, ]
        gm <- s$map; gm$module <- mm$module[match(gm$gene, mm$gene)]
        keep <- !is.na(gm$module)
        if (sum(keep) < 2) next
        al <- dmsa_align(data.frame(cpg = gm$probe[keep],
                                    d = gm$best_direction[keep],
                                    p_plus = gm$p_plus[keep]),
                         genes = gm$gene[keep], level = "gene")
        r <- .report_gene_level(frame, oc, sid,
                                units = gm$module[keep], alignment = al,
                                Ms = frame$M[, gm$column[keep], drop = FALSE])
        r$level <- "module"
        all_res[[paste(oc, "module", sid)]] <- r
        bar$tick()
      }
    }
  }

  RES <- if (length(all_res)) do.call(rbind, all_res) else NULL
  if (!is.null(RES)) {
    ## spec 27: EACH LENS GETS ITS OWN SURVIVOR FLAG. Each per-lens adjusted p
    ## is maxT/minP-controlled within its own lens family, so each flag is an
    ## honest within-lens claim; any_lens_hit and n_lenses_hit summarise them.
    .hit <- function(z) is.finite(z) & z < frame$alpha & RES$n_probes > 0
    RES$selected_coherence <- .hit(RES$p_coherence_adj)
    RES$selected_composite <- .hit(RES$p_composite_adj)
    RES$selected_diffuse   <- .hit(RES$p_diffuse_adj)
    RES$n_lenses_hit <- RES$selected_coherence + RES$selected_composite +
                        RES$selected_diffuse
    RES$any_lens_hit <- RES$n_lenses_hit > 0L
    ## E2 RE-RULED (2026-08-29, second PI ruling, after simulation): a unit is
    ## NAMED by the ANY-LENS rule - some lens's family-adjusted p below alpha
    ## - exactly as the MS states and pre-registers ("survives its family
    ## correction on the lens built for its structure"). The realized
    ## family-wise error of this union rule is disclosed in the MS (.04-.12;
    ## ~.09 under this data's lens dependence, sims in dmsa_patch/
    ## undermine_naming*.R). Every EXACT-alpha alternative (the joint max
    ## statistic, continuous ACAT, Stouffer, hybrids) costs ~10 power points
    ## in every simulated regime - the extra power IS the disclosed
    ## inflation, and the PI ruled for power with disclosure. p_unit_adj (the
    ## exact joint union test) stays computed and printed beside every named
    ## unit as the honesty line, but no longer gates.
    ## RE-RULED (PI, 2026-10-15; ERRCTL_DECISION_v2.md). The any-lens rule
    ## above stays computed and is reported beside the naming rule as the
    ## PRE-REGISTERED rule (`prereg_any_lens`), but it no longer names. On the
    ## pre-registered single engine its coherence lens corrects RAW statistics,
    ## so the family's alpha goes to the units with the widest null (97% to
    ## AVP in the 12-gene oxytocin family), and the union of three lens
    ## families fires in ~13% of null families. A unit is now NAMED when its
    ## dependence-calibrated Benjamini-Hochberg q-value on the unit p (best of
    ## the three lenses, priced by the permutation null of that maximum) is at
    ## or below `fdr` within its level-local family: every unit is judged on
    ## its own calibrated scale, the cross-lens and cross-unit dependence are
    ## measured rather than assumed, and the stated error rate is the true
    ## one (null families: fdr; FDR <= fdr in simulation).
    RES$prereg_any_lens <- RES$any_lens_hit
    .fdr <- frame$fdr %||% .10
    RES$selected <- is.finite(RES$q_unit_dep) & RES$q_unit_dep <= .fdr &
      RES$n_probes > 0
    ## the EXACT-CONFIRMED badge (re-ruled 2026-10-15): a named unit whose
    ## joint union test - the max over the three standardised lenses,
    ## family-corrected on the same stream (p_unit_adj) - also clears alpha
    ## carries an exact, balanced family-wise claim on top of its FDR naming.
    RES$exact_confirmed <- RES$selected &
      is.finite(RES$p_unit_adj) & RES$p_unit_adj < frame$alpha
    ## the OMNIBUS-CONFIRMED badge (PI-approved 2026-08-29): the
    ## family-corrected ACAT omnibus defends against BOTH multiplicities
    ## (cross-lens by the ACAT combination, cross-gene by permutation minP
    ## on the same stream) and rewards cross-lens agreement - the signal
    ## class the union test is weakest for. Additive, like the union badge.
    RES$omnibus_confirmed <- RES$selected &
      is.finite(RES$p_omnibus_adj) & RES$p_omnibus_adj < frame$alpha
    lens <- c("coherence", "composite", "diffuse")
    RES$best_lens <- lens[apply(cbind(RES$p_coherence_adj, RES$p_composite_adj,
                                      RES$p_diffuse_adj), 1, function(z)
                                        if (all(!is.finite(z))) NA else
                                          which.min(z))]
  }
  MOD <- if (length(mod_res)) do.call(rbind, mod_res) else NULL
  SHAPE <- if (length(shape_rows)) do.call(rbind, shape_rows) else NULL
  ## A module whose only direction-called gene is X, and that gene whose only
  ## probe is p, are THE SAME DATA at three level names. They produce identical
  ## statistics and were being listed as three separate survivors, which
  ## silently multiplies the apparent number of findings. Mark the duplicates.
  if (!is.null(MOD) && nrow(MOD)) {
    key <- paste(MOD$outcome, round(MOD$t, 8), round(MOD$b, 8))
    rank <- match(MOD$level, c("system", "module", "gene", "probe"))
    MOD$same_as <- NA_character_
    for (k in unique(key)) {
      i <- which(key == k)
      if (length(i) < 2L) next
      i <- i[order(rank[i])]                       # coarsest level first
      MOD$same_as[i[-1]] <- MOD$unit[i[1]]
      ## A duplicate must be corrected as the unit it duplicates. A gene with a
      ## single direction-called probe creates a PROBE family of one, where maxT
      ## has nothing to correct against, so the identical measurement came out at
      ## the raw p (e.g. .04) while the gene-level twin, corrected across its
      ## 6-gene family, failed at .12. The finer unit then "survived" purely
      ## because its family had one member - and, being a duplicate, it was also
      ## excluded from the figures, so the run reported a survivor it never drew.
      ## Inherit the coarsest twin's adjusted p; keep each unit's own value in a
      ## separate column so the inheritance is auditable.
      for (cc in c("p_composite_adj", "p_curv_adj")) {
        if (!cc %in% names(MOD)) next
        own <- paste0(cc, "_own")
        if (!own %in% names(MOD)) MOD[[own]] <- MOD[[cc]]
        MOD[[cc]][i[-1]] <- MOD[[cc]][i[1]]
      }
    }
    MOD$distinct <- is.na(MOD$same_as)
  }
  if (!is.null(MOD)) {
    MOD$selected <- is.finite(MOD$p_composite_adj) &
                    MOD$p_composite_adj < frame$alpha
    ## The moderated-curvature arm (S^2 x mod) carries its OWN adjusted p.
    ## Keying `selected` on the linear product alone meant a run whose only
    ## survivors were moderated curvature named nine of them in prose and drew
    ## no figure at all - and any user filtering moderation$selected saw none
    ## of them.
    MOD$curv_selected <- if (is.null(MOD$p_curv_adj)) rep(FALSE, nrow(MOD)) else
      is.finite(MOD$p_curv_adj) & MOD$p_curv_adj < frame$alpha
    MOD$any_selected <- MOD$selected | MOD$curv_selected
  }

  ## ---- probe level inside surviving genes -------------------------------
  probe_res <- NULL
  if (!is.null(RES) && "probe" %in% levels_on) {
    hits <- RES[RES$level == "gene" & RES$selected, , drop = FALSE]
    pr_l <- list()
    for (i in seq_len(nrow(hits))) {
      sid <- hits$system_id[i]; g <- hits$unit[i]; oc <- hits$outcome[i]
      s <- frame$sets[[as.character(sid)]]
      k <- s$map$gene == g
      if (sum(k) < 1) next
      al <- dmsa_align(data.frame(cpg = s$map$probe[k],
                                  d = s$map$best_direction[k],
                                  p_plus = s$map$p_plus[k]),
                       genes = s$map$probe[k], level = "gene")
      r <- .report_gene_level(frame, oc, sid, units = s$map$probe[k],
                              alignment = al,
                              Ms = frame$M[, s$map$column[k], drop = FALSE])
      r$level <- "probe"; r$gene <- g
      pr_l[[paste(oc, g)]] <- r
    }
    if (length(pr_l)) probe_res <- do.call(rbind, pr_l)
  }

  ## ---- brain bridge on surviving probes ---------------------------------
  surv_probes <- character(0)
  if (!is.null(RES)) for (i in which(RES$selected & RES$level == "gene")) {
    s <- frame$sets[[as.character(RES$system_id[i])]]
    surv_probes <- c(surv_probes, s$map$probe[s$map$gene == RES$unit[i]])
  }
  bar$tick()                                   # probe level done
  bridge <- .report_bridge(surv_probes)

  ## ---- figures -----------------------------------------------------------
  figs <- character(0); locus_fail <- character(0)
  if (frame$plots && !is.null(RES)) {
    for (oc in frame$outcome) {
      gl <- RES[RES$level == "gene", , drop = FALSE]
      if (nrow(gl[gl$outcome == oc, ])) {
        fp <- file.path(outdir, "figures", paste0("overview_", oc))
        ## spec 29-38: one panel per system, paginated at 40 units - the
        ## function returns every file it actually wrote
        .ovf <- .rp_fig_overview(frame, gl, oc, fp)
        if (length(.ovf)) figs <- c(figs, .ovf)
      }
      ## module level gets its own panel: it is a declared level of the
      ## hierarchy and was the only one with no figure
      if (any(RES$level == "module" & RES$outcome == oc)) {
        mp <- file.path(outdir, "figures", paste0("modules_", oc))
        if (!is.null(.rp_fig_modules(frame, RES, oc, mp)))
          figs <- c(figs, paste0(mp, ".", frame$plot_type))
      }
    }
    ## EVERY named gene gets its locus panel (PI, 2026-08-29). .rp_fig_locus
    ## falls back to the bare panel on any failure of the full one; a gene
    ## that still could not be drawn is collected here and NAMED in
    ## summary.md, so a missing figure is a stated fact, never a hole.
    hits <- RES[RES$level == "gene" & RES$selected, , drop = FALSE]
    for (i in seq_len(nrow(hits))) {
      fp <- file.path(outdir, "figures",
                      paste0("locus_", hits$unit[i], "_", hits$outcome[i]))
      .lf <- .rp_fig_locus(frame, hits$outcome[i], hits$system_id[i],
                           hits$unit[i], TRUE, fp,
                           fits = probe_fits[[paste(hits$outcome[i],
                                                    hits$system_id[i])]])
      if (!is.null(.lf)) {
        figs <- c(figs, paste0(fp, ".", frame$plot_type))
        if (!is.null(attr(.lf, "probe_rows")))
          locus_rows[[paste(hits$outcome[i], hits$unit[i])]] <-
            attr(.lf, "probe_rows")
      } else locus_fail <- c(locus_fail,
                           sprintf("%s (%s)", hits$unit[i],
                                   .lab(frame, hits$outcome[i])))
      ## A binary outcome is non-linear by construction, so every surviving unit
      ## gets its logistic response curve - the effect itself, drawn.
      kk <- (frame$outcome_kind %||% list())[[hits$outcome[i]]]
      if (!is.null(kk) && kk$kind %in% c("binary", "wave")) {
        fb <- file.path(outdir, "figures",
                        paste0("response_",
                               .rp_unit_file(frame, hits$unit[i], "gene"),
                               "_", hits$outcome[i]))
        r <- .rp_try("figure: response curve",
                     .rp_fig_binary(frame, hits$outcome[i], "gene",
                                    hits$family[i] %||%
                                      as.character(hits$system[i]),
                                    hits$unit[i], fb),
                     outcome = hits$outcome[i], level = "gene",
                     unit = hits$unit[i], fallback = "no response-curve figure")
        if (!is.null(r)) figs <- c(figs, paste0(fb, ".", frame$plot_type))
      }
    }
  }
  ## A surviving interaction is the result that most needs a picture and had
  ## none: one panel per survivor, simple slopes over their own partial
  ## residuals plus Johnson-Neyman.
  if (frame$plots && !is.null(MOD)) {
    ## draw linear-product AND curvature survivors; `distinct` keeps a unit that
    ## repeats the same data at another level from getting a duplicate panel
    mh <- MOD[MOD$any_selected &
              (if (is.null(MOD$distinct)) TRUE else MOD$distinct), ,
              drop = FALSE]
    for (i in seq_len(nrow(mh))) {
      fp <- file.path(outdir, "figures",
                      paste0("moderation_", gsub("[^A-Za-z0-9]+", "_",
                               .rp_unit_file(frame, mh$unit[i], mh$level[i])),
                             "_", mh$outcome[i]))
      if (!is.null(.rp_fig_moderation(frame, mh[i, ], fp)))
        figs <- c(figs, paste0(fp, ".", frame$plot_type))
      ## companion surface for a moderated non-linear model
      fps <- sub("^moderation_", "surface_", basename(fp))
      fps <- file.path(dirname(fp), fps)
      if (!is.null(.rp_try("figure: moderation surface",
                           .rp_fig_surface(frame, mh[i, ], fps),
                           outcome = mh$outcome[i], level = mh$level[i],
                           unit = mh$unit[i], fallback = "no surface figure")))
        figs <- c(figs, paste0(fps, ".", frame$plot_type))
    }
  }
  ## type != "linear" asserts something about functional form. Draw it.
  if (frame$plots && !is.null(SHAPE)) {
    ## Draw only what survived its family: a figure per scanned unit would be
    ## dozens of panels, and the ones that matter are the ones a reader would
    ## otherwise have to take on trust.
    ## Gate on ANY surviving departure arm. Keying this on the quadratic alone
    ## meant an exponential run whose own arm survived drew no figure at all.
    adjm <- cbind(SHAPE$p_quad_adj, SHAPE$p_thr_adj, SHAPE$p_expo_adj)
    keep <- which(apply(adjm, 1, function(z)
      any(is.finite(z) & z < frame$alpha)))
    for (i in keep) {
      fp <- file.path(outdir, "figures",
                      paste0("shape_", SHAPE$level[i], "_",
                             gsub("[^A-Za-z0-9]+", "_",
                               .rp_unit_file(frame, SHAPE$unit[i], SHAPE$level[i])),
                             "_", SHAPE$outcome[i]))
      if (!is.null(.rp_fig_shape(frame, SHAPE$outcome[i], SHAPE$level[i],
                                 SHAPE$family[i], SHAPE$unit[i], SHAPE[i, ], fp)))
        figs <- c(figs, paste0(fp, ".", frame$plot_type))
    }
    if (!length(keep))
      message("shape scan: no unit's departure term survived its level-local ",
              "family correction, so no shape figure was drawn. ",
              "tables/shape.csv has every unit's raw and adjusted p.")
  }

  bar$tick()                                   # figures done
  ## ---- tables ------------------------------------------------------------
  tabs <- character(0)
  if (frame$tables) {
    if (!is.null(RES)) {
      keep <- intersect(c("level", "outcome", "system", "unit", "n_probes",
                          "concordance", "direction", "p_coherence",
                          "p_coherence_adj", "p_composite", "p_composite_adj",
                          "p_diffuse", "p_diffuse_adj",
                          "p_unit", "p_unit_adj", "p_omnibus",
                          "best_lens", "selected_coherence",
                          "selected_composite", "selected_diffuse",
                          "n_lenses_hit", "any_lens_hit", "prereg_any_lens",
                          "q_unit", "q_unit_dep", "null_alarm_bh", "selected",
                          "p_union_exact", "exact_confirmed",
                          "p_omnibus_adj", "omnibus_confirmed",
                          "fwer_realized"),
                        names(RES))
      .rp_write_table(.rp_add_unit_label(frame, RES[, keep]), file.path(outdir, "tables", "units"),
                      frame, "DMSA - every unit, three lenses")
      hits <- RES[RES$selected, keep, drop = FALSE]
      if (!is.null(bridge$table) && isTRUE(bridge$any_hit) && nrow(hits))
        hits$brain_bridge <- vapply(seq_len(nrow(hits)), function(i) {
          if (hits$level[i] != "gene") return(NA_real_)
          sid <- RES$system_id[RES$selected][i]
          if (is.na(sid)) return(NA_real_)
          s <- frame$sets[[as.character(sid)]]
          pb <- s$map$probe[s$map$gene == hits$unit[i]]
          sum(bridge$table$bridge_hit[bridge$table$cpg_id %in% pb])
        }, numeric(1))
      .rp_write_table(.rp_add_unit_label(frame, hits), file.path(outdir, "tables", "hits"), frame,
                      "DMSA - units surviving their level-local family")
      tabs <- c(tabs, file.path(outdir, "tables", c("units", "hits")))
    }
    if (!is.null(MOD)) {
      .rp_write_table(.rp_add_unit_label(frame, MOD), file.path(outdir, "tables", "moderation"), frame,
                      sprintf("DMSA - frame x %s%s (composite lens)",
                              frame$mod,
                              if (nzchar(frame$mod2))
                                paste0(" x ", frame$mod2) else ""))
      tabs <- c(tabs, file.path(outdir, "tables", "moderation"))
    }
    if (!is.null(SHAPE)) {
      .rp_write_table(.rp_add_unit_label(frame, SHAPE), file.path(outdir, "tables", "shape"), frame,
                      sprintf("DMSA - shape scan (%s), system-level aligned score",
                              frame$type))
      tabs <- c(tabs, file.path(outdir, "tables", "shape"))
    }
    ## spec 52: the numbers every locus panel drew, so figure and test can
    ## be audited against each other row by row
    if (length(locus_rows)) {
      .lr <- do.call(rbind, locus_rows); rownames(.lr) <- NULL
      utils::write.csv(.lr, file.path(outdir, "tables", "locus_probes.csv"),
                       row.names = FALSE)
      tabs <- c(tabs, file.path(outdir, "tables", "locus_probes"))
    }
    ## spec 40: per-gene polarity status behind every system score
    .pa_t <- .rp_polarity_audit(frame)
    if (!is.null(.pa_t) && isTRUE(attr(.pa_t, "has_polarity"))) {
      utils::write.csv(.pa_t, file.path(outdir, "tables", "polarity_audit.csv"),
                       row.names = FALSE)
      tabs <- c(tabs, file.path(outdir, "tables", "polarity_audit"))
    }
    ## spec 53: the QC ledger is ALWAYS written when tables are on - a
    ## header-only file states "nothing failed" as plainly as rows state
    ## what did. Written here, after every figure and table above it.
    utils::write.csv(.rp_qc_rows(),
                     file.path(outdir, "tables", "qc_report.csv"),
                     row.names = FALSE)
    tabs <- c(tabs, file.path(outdir, "tables", "qc_report"))
    if (!is.null(probe_res)) {
      .rp_write_table(probe_res, file.path(outdir, "tables", "probes"),
                      frame, "DMSA - probe level within surviving genes")
      ## probes.csv was written but never added to `tabs`, so the console line
      ## under-reported the table count by one in every run that produced it
      tabs <- c(tabs, file.path(outdir, "tables", "probes"))
    }
    ## spec 17: on the pair path, every report carries the FULL pair ledger -
    ## every discovered CpG x gene pair, used or not, with its reason. This is
    ## the table that makes a silent drop impossible: a gene that lost all of
    ## its evidence is visible here even when nothing else mentions it.
    if (!is.null(frame$pair_ledger)) {
      utils::write.csv(frame$pair_ledger,
                       file.path(outdir, "tables", "cpg_gene_pair_ledger.csv"),
                       row.names = FALSE)
      tabs <- c(tabs, file.path(outdir, "tables", "cpg_gene_pair_ledger"))
    }
    ## the ANALYSIS SET: the probes actually tested, annotated, one row per
    ## CpG x gene pair - the supplement table a paper needs and the file a
    ## user validates the selection against (PI, 2026-08-29). The ledger
    ## above additionally holds every pair NOT used; this one is only what
    ## entered the analysis. Both paths (pair and bundled) get it.
    tryCatch({
      dmsa_save_analysis_set(frame,
                             file.path(outdir, "tables", "analysis_set.csv"))
      tabs <- c(tabs, file.path(outdir, "tables", "analysis_set"))
    }, error = function(e)
      message("analysis_set.csv skipped: ", conditionMessage(e)))
  }

  ## ---- summary.md --------------------------------------------------------
  smy <- file.path(outdir, "summary.md")
  if (frame$summary) {
    n_hit <- if (!is.null(RES)) sum(RES$selected & RES$level == "gene") else 0
    md <- c(sprintf("# DMSA report - %s",
                    paste(vapply(frame$outcome, function(.o) .lab(frame, .o),
                                 character(1)), collapse = ", ")),
      "",
      sprintf("**n = %d. %d system(s), %d mapped probes (cpg_map = '%s'). %s within level-local families, B = %d, seed %d.%s**",
              nrow(frame$data), nrow(frame$systems), nrow(frame$map),
              frame$cpg_map, frame$correction, frame$B, frame$seed,
              if (frame$frame_role == "outcome")
                " frame_role = outcome: the named columns PREDICT the methylation frame; direction language is expression tone." else ""),
      "", "## Summary", "")
    if (frame$moderation && !is.null(MOD)) {
      mh <- MOD[MOD$selected, , drop = FALSE]
      tp <- MOD[which.min(MOD$p_composite), ]
      md <- c(md, sprintf(
        "Moderation (composite lens - the one lens that carries to products): %d unit(s) tested across %d level-local famil%s for frame x %s%s, B = %d each. %s. Strongest interaction: %s (%s level, b = %+.3f, t = %+.2f, raw p %.4f, family-adjusted %.4f).%s",
        nrow(MOD), length(unique(MOD$family)),
        if (length(unique(MOD$family)) == 1) "y" else "ies",
        frame$mod,
        if (nzchar(frame$mod2)) paste0(" x ", frame$mod2) else "", frame$B,
        if (nrow(mh)) paste0(nrow(mh), " survive",
          { nd <- sum(mh$distinct, na.rm = TRUE)
            if (nd < nrow(mh)) sprintf(
              " (%d distinct; the rest are the SAME data at another level - a module whose only direction-called gene is that gene, or a gene whose only probe is that probe)",
              nd) else "" }, ": ",
          paste(sprintf("%s (%s, adj %.4f%s%s)", .rp_unit_label(frame, mh$unit, mh$level),
                        .lab(frame, mh$outcome), mh$p_composite_adj,
                        ifelse(!is.na(mh$n_family) & mh$n_family == 1L,
                               " - FAMILY OF 1, so uncorrected", ""),
                        ifelse(!is.na(mh$same_as),
                               paste0(" = ", mh$same_as), "")),
                collapse = "; "))
        else "None survives its family correction - a calibrated null, not a failed run",
        .rp_unit_label(frame, tp$unit, tp$level), tp$level, tp$b, tp$t, tp$p_composite, tp$p_composite_adj,
        if (!nrow(mh)) " The full per-unit table is tables/moderation.csv."
        else ""))
      if (frame$type != "linear") {
        ch <- MOD[is.finite(MOD$p_curv_adj) & MOD$p_curv_adj < frame$alpha, ,
                  drop = FALSE]
        tq <- if (any(is.finite(MOD$p_curv)))
          MOD[which.min(MOD$p_curv), ] else NULL
        md <- c(md, "", sprintf(paste0(
          "Moderated non-linearity (`type = \"%s\"`): the same %d unit(s) were ",
          "also tested for whether the CURVATURE of the tone score depends on ",
          "%s - tested term %s x %s, with S, S^2, %s and S x %s all in the ",
          "model. %s%s"),
          frame$type, nrow(MOD), .lab(frame, frame$mod), "S^2",
          .lab(frame, frame$mod), .lab(frame, frame$mod), .lab(frame, frame$mod),
          if (nrow(ch)) paste0(nrow(ch), " survive",
            { nd <- sum(ch$distinct, na.rm = TRUE)
              if (nd < nrow(ch)) sprintf(
                " (%d distinct; the rest repeat the same data at another level)",
                nd) else "" }, ": ",
            paste(sprintf("%s (%s, adj %.4f%s)", .rp_unit_label(frame, ch$unit, ch$level),
                          .lab(frame, ch$outcome),
                          ch$p_curv_adj,
                          ifelse(!is.na(ch$same_as), paste0(" = ", ch$same_as), "")),
                  collapse = "; "))
          else "None survives its family correction",
          if (!is.null(tq) && is.finite(tq$p_curv))
            sprintf(". Strongest: %s (%s level, b = %+.3f, t = %+.2f, raw p %.4f, adjusted %.4f).",
                    .rp_unit_label(frame, tq$unit, tq$level), tq$level, tq$b_curv, tq$t_curv, tq$p_curv,
                    tq$p_curv_adj) else "."))
      }
    } else if (!is.null(RES)) {
      hits <- RES[RES$selected & RES$level == "gene", , drop = FALSE]
      .fdr <- frame$fdr %||% .10
      ## the naming rule (PI ruling 2026-10-15): dependence-calibrated BH on
      ## the unit p within each level-local family. The carrying lens is still
      ## stated, the exact balanced union test (p_unit_adj) rides along as the
      ## family-wise badge, and the pre-registered any-lens rule is reported
      ## right after, so a reader sees both decisions on the same evidence.
      md <- c(md, if (nrow(hits)) sprintf(
        "%d gene(s) named at FDR %.2g (Benjamini-Hochberg on each gene's best-of-three-lens permutation p within its level-local family, level calibrated to the family's own dependence; q below is that calibrated q): %s.",
        nrow(hits), .fdr,
        paste(sprintf("%s (%s, q = %.3f, unit p %.4f, carried by the %s lens; joint union adj %.4f; family-corrected omnibus %.4f%s%s)",
                      .rp_unit_label(frame, hits$unit, hits$level),
                      vapply(hits$outcome, function(.o) .lab(frame, .o), character(1)),
                      hits$q_unit_dep, hits$p_unit, hits$best_lens,
                      hits$p_unit_adj, hits$p_omnibus_adj,
                      ifelse(hits$exact_confirmed %in% TRUE,
                             " - EXACT-CONFIRMED", ""),
                      ifelse(hits$omnibus_confirmed %in% TRUE,
                             " - OMNIBUS-CONFIRMED", "")),
              collapse = "; "))
        else sprintf("No gene is named at FDR %.2g: no gene's calibrated q reaches the bar in its family.", .fdr))
      ## the pre-registered any-lens rule, on the same evidence, beside it
      .pr <- RES[RES$level == "gene" & RES$prereg_any_lens %in% TRUE, , drop = FALSE]
      md <- c(md, "", if (nrow(.pr)) sprintf(
        "Pre-registered any-lens rule (some lens's family-adjusted p < %.2g, %s within each lens's own family), reported beside the naming rule: %s.",
        frame$alpha, frame$correction,
        paste(sprintf("%s (%s, %s lens, lens adj %.4f%s)",
                      .rp_unit_label(frame, .pr$unit, .pr$level),
                      vapply(.pr$outcome, function(.o) .lab(frame, .o), character(1)),
                      .pr$best_lens,
                      pmin(.pr$p_coherence_adj, .pr$p_composite_adj,
                           .pr$p_diffuse_adj, na.rm = TRUE),
                      ifelse(.pr$selected %in% TRUE, "; also named", "; not named at the FDR level")),
              collapse = "; "))
        else sprintf("Pre-registered any-lens rule (some lens's family-adjusted p < %.2g), reported beside the naming rule: no gene passes.", frame$alpha))
      ## the run's own measured error rates, per gene family
      .fw <- unique(RES[RES$level == "gene" & (is.finite(RES$fwer_realized) |
                          is.finite(RES$null_alarm_bh %||% NA_real_)),
                        intersect(c("system", "outcome", "fwer_realized",
                                    "null_alarm_bh"), names(RES)), drop = FALSE])
      if (nrow(.fw))
        md <- c(md, "", sprintf(
          paste0("Error rates measured from this run's own permutation ",
                 "stream, per gene family (system/outcome: realized ",
                 "family-wise error of the pre-registered any-lens rule; ",
                 "null false-alarm of ordinary BH at %.2g, before the ",
                 "dependence calibration that sets it to %.2g): %s. Two ",
                 "badges ride on a named gene, each an exact family-wise ",
                 "%.2g claim: EXACT-CONFIRMED (dagger) = the joint union ",
                 "test (max over the three standardised lenses, ",
                 "family-corrected) clears alpha; OMNIBUS-CONFIRMED (double ",
                 "dagger) = the family-corrected ACAT omnibus clears alpha ",
                 "(rewards cross-lens agreement)."),
          .fdr, .fdr,
          paste(sprintf("%s/%s %.3f, %s", .fw$system,
                        vapply(.fw$outcome, function(.o) .lab(frame, .o),
                               character(1)),
                        .fw$fwer_realized,
                        if (is.null(.fw$null_alarm_bh)) "-" else
                          sprintf("%.3f", .fw$null_alarm_bh)),
                collapse = "; "),
          frame$alpha))
      ## every named gene's panel is pointed to by name - and if one could
      ## not be drawn even by the bare fallback, that is said here, not
      ## discovered as a missing file
      if (nrow(hits) && isTRUE(frame$plots)) {
        .drawn <- !(sprintf("%s (%s)", hits$unit,
                            vapply(hits$outcome, function(.o) .lab(frame, .o),
                                   character(1))) %in% locus_fail)
        if (any(.drawn))
          md <- c(md, "", sprintf(
            "Each named gene has its locus panel: %s. The bars in panel C are the gene test's own per-probe estimates (winsorised, standardised methylation per SD of the tested term - the same preprocessing, adjustment and random-intercept transform as the test); the numbers drawn are in `tables/locus_probes.csv`, one row per probe, so figure and test can be checked against each other.",
            paste(sprintf("`figures/locus_%s_%s.%s`", hits$unit[.drawn],
                          hits$outcome[.drawn], frame$plot_type),
                  collapse = ", ")))
        if (length(locus_fail))
          md <- c(md, "", sprintf(
            paste0("NOTE: the locus panel(s) for %s could not be drawn even ",
                   "in bare form (the reason was printed to the console at ",
                   "run time); every other output for these named gene(s) ",
                   "is present."),
            paste(locus_fail, collapse = "; ")))
      }

      ## GENE RESULTS are stated even when nothing survives (PI, 2026-08-29):
      ## a summary that jumps from "no survivors" to the module table leaves
      ## the reader with no idea what the genes DID. Top of each system's
      ## family by the naming statistic, with the full table's location.
      gl_all <- RES[RES$level == "gene" & RES$n_probes > 0, , drop = FALSE]
      if (nrow(gl_all)) {
        md <- c(md, "", "### Gene results", "", paste0(
          "A gene's significance comes from its CpGs jointly, not from any ",
          "one of them: each lens pools small, direction-consistent shifts ",
          "across the gene's CpGs, so a gene can survive while no single ",
          "CpG is significant on its own. When that happens the gene's ",
          "locus figure draws every CpG in grey and says so under its ",
          "title - grey means \"below the noise line individually\", not ",
          "\"excluded from the result\".",
          if (isTRUE(frame$tables))
            paste0(" The exact CpGs tested per gene, with their annotation, ",
                   "are in `tables/analysis_set.csv`.") else ""), "")
        ## ONE bar everywhere (PI, 2026-08-29; statistic re-ruled 2026-10-15):
        ## the table is ordered by the SAME statistic that names a finding -
        ## the dependence-calibrated q. Bold = named. The pre-registered
        ## any-lens column (best lens's family-adjusted p) sits beside it, and
        ## the raw omnibus says (raw) in its own header so an uncorrected
        ## number can never read as a verdict.
        gl_all$.lens_adj <- pmin(gl_all$p_coherence_adj,
                                 gl_all$p_composite_adj,
                                 gl_all$p_diffuse_adj, na.rm = TRUE)
        for (.oc2 in unique(gl_all$outcome)) {
          g2 <- gl_all[gl_all$outcome == .oc2, , drop = FALSE]
          for (.sy in unique(g2$system)) {
            gs2 <- g2[g2$system == .sy, , drop = FALSE]
            gs2 <- gs2[order(gs2$q_unit_dep, gs2$p_unit, gs2$p_omnibus),
                       , drop = FALSE]
            topn <- utils::head(gs2, 10L)
            md <- c(md, sprintf(
              "**%s** - %s: %d gene(s) tested. Bold = named (calibrated q <= %.2g); \u2020 = exact-confirmed (joint union adj < %.2g); \u2021 = omnibus-confirmed (family-corrected ACAT); both are exact family-wise %.2g claims. The pre-registered any-lens column is the best lens's family-adjusted p (%s within the lens).%s",
              .sy, .lab(frame, .oc2), nrow(gs2), frame$fdr %||% .10,
              frame$alpha, frame$alpha, frame$correction,
              if (nrow(gs2) > nrow(topn))
                sprintf(" Top %d by the naming statistic; all %d are in `tables/units.csv`.",
                        nrow(topn), nrow(gs2)) else ""), "",
              "| gene | CpGs | concord. | dir | best lens | q (calibrated) | unit p | joint union adj | pre-registered: lens adj | omnibus (raw) | omnibus adj |",
              "|---|---|---|---|---|---|---|---|---|---|---|",
              sprintf("| %s%s%s%s | %d | %.2f | %s | %s | %s | %s | %s | %s | %.4f | %s |",
                      ifelse(topn$selected, "**", ""),
                      paste0(topn$unit, ifelse(topn$selected, "**", "")),
                      ifelse(topn$exact_confirmed %in% TRUE, " \u2020", ""),
                      ifelse(topn$omnibus_confirmed %in% TRUE, " \u2021", ""),
                      topn$n_probes, topn$concordance,
                      ifelse(is.na(topn$direction), "-",
                             ifelse(topn$direction > 0, "+1", "-1")),
                      topn$best_lens,
                      sprintf("%.3f", topn$q_unit_dep),
                      sprintf("%.4f", topn$p_unit),
                      sprintf("%.4f", topn$p_unit_adj),
                      sprintf("%.4f%s", topn$.lens_adj,
                              ifelse(topn$prereg_any_lens %in% TRUE, " *", "")),
                      topn$p_omnibus,
                      sprintf("%.4f", topn$p_omnibus_adj)), "")
          }
        }
      }

      ## MODULE is a declared level of the hierarchy, and its survivors were
      ## written to hits.csv but named nowhere in the prose. In one battery run
      ## four modules survived - the strongest at adj 0.002, better than two of
      ## the three genes that WERE reported - and a reader of summary.md alone
      ## would never have known. Say them.
      modr <- RES[RES$selected & RES$level == "module", , drop = FALSE]
      if (nrow(modr)) {
        md <- c(md, "", sprintf(
          "%d module(s) named at FDR %.2g within their own system's family: %s.",
          nrow(modr), frame$fdr %||% .10,
          paste(sprintf("%s (%s, q = %.3f, carried by the %s lens; pre-registered lens adj %.4f; joint union adj %.4f)",
                        modr$unit,
                        vapply(modr$outcome, function(.o) .lab(frame, .o), character(1)),
                        modr$q_unit_dep, modr$best_lens,
                        pmin(modr$p_coherence_adj, modr$p_composite_adj,
                             modr$p_diffuse_adj, na.rm = TRUE),
                        modr$p_unit_adj), collapse = "; ")))
      }
      sysr <- RES[RES$level == "system", , drop = FALSE]
      ## the system line leads with the SAME naming statistic as every other
      ## level (re-ruled 2026-10-15: the calibrated q), with the
      ## pre-registered best-lens adjusted p beside it; the star marks
      ## naming, and the raw omnibus is labelled raw so it cannot read as
      ## the verdict
      if (nrow(sysr)) {
        .sl <- pmin(sysr$p_coherence_adj, sysr$p_composite_adj,
                    sysr$p_diffuse_adj, na.rm = TRUE)
        md <- c(md, "",
        sprintf("System level (family = the %d named system(s); * = named at FDR %.2g): %s.",
                nrow(frame$systems), frame$fdr %||% .10,
                paste(sprintf("%s/%s q %.3f%s (unit p %.4f, carried by %s; pre-registered lens adj %.4f; omnibus raw %.4f; joint union adj %.4f)",
                              sysr$unit,
                              vapply(sysr$outcome, function(.o) .lab(frame, .o),
                                     character(1)),
                              sysr$q_unit_dep, ifelse(sysr$selected, " *", ""),
                              sysr$p_unit, sysr$best_lens, .sl,
                              sysr$p_omnibus, sysr$p_unit_adj),
                      collapse = "; ")))
      }
    }
    ## spec 40: say what the system tone was built from - a system score
    ## from 30 of 41 testable genes must not read like one from 41
    .pa <- .rp_polarity_audit(frame)
    if (!is.null(.pa) && isTRUE(attr(.pa, "has_polarity")) &&
        "system" %in% levels_on && !sys_skip) {
      md <- c(md, "", sprintf(paste0(
        "System-level composition (spec 40): %s. A gene with unresolved ",
        "polarity has no curated activation sign and is weighted 0 in the ",
        "system score - it is NOT assumed +1; it simply does not contribute. ",
        "Full per-gene status: tables/polarity_audit.csv."),
        paste(sprintf("%s - tone built from %d of %d testable gene(s)%s%s",
                      .pa$system, .pa$n_polarity_signed, .pa$n_genes_testable,
                      ifelse(.pa$n_polarity_off_axis > 0,
                             sprintf(", %d explicitly off-axis (w_g = 0)",
                                     .pa$n_polarity_off_axis), ""),
                      ifelse(.pa$n_polarity_unresolved > 0,
                             sprintf(", %d with UNRESOLVED polarity (%s)",
                                     .pa$n_polarity_unresolved,
                                     gsub(";", ", ", .pa$genes_unresolved)),
                             "")),
              collapse = "; ")))
    }
    md <- c(md, "", sprintf("Brain bridge: %s.", bridge$status), "")
    ## The shape scan, said out loud. It used to be computed and discarded, so a
    ## run with type = "non-linear" returned numbers identical to the linear run
    ## and never mentioned it - the user had no way to know.
    if (!is.null(SHAPE)) {
      md <- c(md, "", "## Shape scan", "",
        sprintf(paste0("`type = \"%s\"` scans functional form at every declared ",
          "level (%s). Each departure term is fitted **with the linear term in ",
          "the model**, so its p is the 1-df incremental test of curvature over ",
          "and above the linear component, and the quadratic term is corrected ",
          "by Westfall-Young maxT **inside its level-local family**, on the same ",
          "shared permutation stream used everywhere else - as is every other ",
          "departure arm (B = %d per arm)."),
          frame$type,
          paste(unique(SHAPE$level), collapse = " > "),
          max(499L, frame$B %/% 4L)),
        "",
        paste0("**This does not change any unit-level DMSA result.** The ",
               "p-values reported above are the linear tests. The scan is a ",
               "shape diagnostic on each unit's aligned score."),
        "",
        if (any(is.finite(SHAPE$p_expo)))
          "| outcome | level | unit | linear | +quad (adj) | +thr (adj) | +exp (adj) | ACAT | shape test | extremum [95% Fieller] |"
        else
          "| outcome | level | unit | linear | +quad (adj) | +thr (adj) | ACAT | shape test | extremum [95% Fieller] |",
        if (any(is.finite(SHAPE$p_expo))) "|---|---|---|---|---|---|---|---|---|---|"
        else "|---|---|---|---|---|---|---|---|---|")
      fmtp <- function(x) if (!is.finite(x)) "n/a" else
        sprintf("%.4f%s", x, if (x < frame$alpha) " *" else "")
      ## Rank by the best adjusted departure p across ARMS. Ranking on the
      ## quadratic alone pushed units that survived on the threshold or
      ## exponential arm off the bottom of the table while they still got a
      ## figure - the table and the figures then disagreed.
      best <- apply(cbind(SHAPE$p_quad_adj, SHAPE$p_thr_adj, SHAPE$p_expo_adj),
                    1, function(z) if (all(!is.finite(z))) NA_real_ else
                      min(z, na.rm = TRUE))
      ord <- order(best, SHAPE$p_nonlin)
      surv <- which(is.finite(best) & best < frame$alpha)
      shown <- union(surv, utils::head(ord, 12L))
      shown <- shown[order(match(shown, ord))]
      for (i in shown) {
        shp <- if (!is.finite(SHAPE$p_ushape[i])) "n/a" else
          if (SHAPE$p_ushape[i] < frame$alpha && isTRUE(SHAPE$turn_inside[i]))
            sprintf("%s p = %.4f *%s", SHAPE$shape[i], SHAPE$p_ushape[i],
                    if (is.finite(SHAPE$p_ushape_trim[i]) &&
                        SHAPE$p_ushape_trim[i] >= frame$alpha)
                      sprintf(" (p = %.3f on central 90%%)",
                              SHAPE$p_ushape_trim[i]) else "")
          else sprintf("no %s (p = %.3f)", SHAPE$shape[i], SHAPE$p_ushape[i])
        ci <- if (!is.finite(SHAPE$turn[i])) "n/a" else
          if (!is.finite(SHAPE$turn_lo[i]) || !is.finite(SHAPE$turn_hi[i]))
            sprintf("%.2f [unbounded]", SHAPE$turn[i])
          else sprintf("%.2f [%.2f, %.2f]%s", SHAPE$turn[i], SHAPE$turn_lo[i],
                       SHAPE$turn_hi[i],
                       if (isTRUE(SHAPE$turn_ci_inside[i])) "" else " (not contained)")
        pa <- function(r, a) sprintf("%s (%s)", fmtp(r), fmtp(a))
        md <- c(md, if (any(is.finite(SHAPE$p_expo)))
          sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
            SHAPE$outcome[i], SHAPE$level[i], SHAPE$unit[i],
            fmtp(SHAPE$p_lin[i]), pa(SHAPE$p_quad[i], SHAPE$p_quad_adj[i]),
            pa(SHAPE$p_thr[i], SHAPE$p_thr_adj[i]),
            pa(SHAPE$p_expo[i], SHAPE$p_expo_adj[i]),
            fmtp(SHAPE$p_nonlin[i]), shp, ci)
          else
          sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |",
            SHAPE$outcome[i], SHAPE$level[i], SHAPE$unit[i],
            fmtp(SHAPE$p_lin[i]), pa(SHAPE$p_quad[i], SHAPE$p_quad_adj[i]),
            pa(SHAPE$p_thr[i], SHAPE$p_thr_adj[i]),
            fmtp(SHAPE$p_nonlin[i]), shp, ci))
        if (isTRUE(SHAPE$n_family[i] == 1L))
          md[length(md)] <- sub("\\|$", "| <- family of 1: no correction possible",
                                md[length(md)])
      }
      if (length(ord) > length(shown))
        md <- c(md, sprintf(
          if (any(is.finite(SHAPE$p_expo)))
            "| ... | | %d further unit(s) in `tables/shape.csv` | | | | | | | |"
          else
            "| ... | | %d further unit(s) in `tables/shape.csv` | | | | | | |",
          length(ord) - length(shown)))
      md <- c(md, "",
        paste0("**Reading the shape columns.** A significant quadratic term is ",
               "*not* evidence of a U or an inverted U: a relationship that ",
               "rises throughout but with increasing steepness produces one ",
               "too, and the quadratic approximation then invents a turning ",
               "point (Lind & Mehlum 2010; Simonsohn 2018). The shape column ",
               "is the Sasabuchi intersection-union test - the fitted slope ",
               "must be **negative at the low end and positive at the high end ",
               "of the observed range**, with the p-value taken as the larger ",
               "of the two one-sided p-values so both ends must hold. The ",
               "extremum carries a **Fieller** 95% interval; Lind & Mehlum ",
               "treat its containment in the observed range as the equivalent ",
               "test and warn that the delta method is severely biased in ",
               "finite samples. An unbounded interval means the turning point ",
               "is simply not determined. Claim a U only when the adjusted ",
               "quadratic p, the Sasabuchi test and the Fieller containment ",
               "all agree; otherwise report a curvilinear component."),
        "",
        if (any(is.finite(SHAPE$p_thr_adj) & SHAPE$p_thr_adj < frame$alpha))
          paste0("**On the threshold arm, read with caution.** Its break point ",
                 "is fixed at the score's own mean (score > 0); it is NOT ",
                 "estimated from the data. An arbitrary break tested against ",
                 "the same data that defined it is exactly the design ",
                 "Simonsohn (2018) shows inflates false positives, and his ",
                 "two-lines procedure exists to choose the break properly. ",
                 "Treat a threshold-only survivor as a flag to look at the ",
                 "figure, not as a result - and note that the Sasabuchi and ",
                 "Fieller columns describe the QUADRATIC fit, so they say ",
                 "nothing about a step.") else NULL,
        if (any(is.finite(SHAPE$p_thr_adj) & SHAPE$p_thr_adj < frame$alpha))
          "" else NULL,
        if (any(is.finite(SHAPE$p_expo)))
          paste0("On the exponential arm: exp(score) is a monotone transform ",
                 "and correlates strongly with the linear score, so the ",
                 "incremental test has low power and a null there is weak ",
                 "evidence. It answers whether the association accelerates, ",
                 "not whether it is non-monotone - only the shape test speaks ",
                 "to that.") else NULL,
        if (any(is.finite(SHAPE$p_expo))) "" else NULL,
        paste0("Two cautions. The Sasabuchi and Fieller results use parametric ",
               "standard errors, not the block permutation, so under clustering ",
               "(here: ", if (length(frame$block_cols))
                 paste(frame$block_cols, collapse = ", ") else "none",
               ") treat them as anticonservative - only the quadratic p and its ",
               "maxT adjustment are permutation-based. And curvature in a ",
               "multivariable model can be induced by an omitted product term ",
               "(Ganzach 1997), so re-check a curvilinear result with the ",
               "relevant interaction in the model before interpreting it."),
        "",
        paste0("Full table: `tables/shape.csv`. A figure is drawn for any unit ",
               "whose adjusted p clears alpha on ANY departure arm - quadratic, ",
               "threshold or exponential - not on the quadratic alone."), "")
    }
    ## Module-level evidence. A module result is only as good as the module's
    ## definition, so the strength of that definition and the references behind
    ## it are reported next to the result rather than left in a data file.
    md <- c(md, .report_evidence_md(frame, RES))
    md <- c(md, "## Design notes", "")
    if (nrow(frame$corrections))
      md <- c(md, "Test-drive corrections:",
              sprintf("- %s: %s -> %s", frame$corrections$field,
                      frame$corrections$issue, frame$corrections$action), "")
    if (nrow(frame$map_conflicts))
      md <- c(md,
        sprintf("Direction-map disagreements (confidence vs full) - analysed under '%s', flip elsewhere: %s.",
                frame$cpg_map,
                paste(sprintf("%s %s (%+d vs %+d)",
                              frame$map_conflicts$level,
                              frame$map_conflicts$unit,
                              frame$map_conflicts$sign_confidence,
                              frame$map_conflicts$sign_full), collapse = "; ")),
        "")
    md <- c(md, sprintf(
      "Coverage: %d gene(s) in the chosen systems carry no usable direction-called probe under cpg_map = '%s' and are untestable%s.",
      nrow(frame$untestable), frame$cpg_map,
      if (nrow(frame$untestable))
        paste0(" (", paste(utils::head(frame$untestable$gene, 12),
                           collapse = ", "),
               if (nrow(frame$untestable) > 12) ", ..." else "", ")")
      else ""),
      { mm <- as.numeric(Sys.time() - t_start, units = "mins")
        sprintf("Runtime: %s (pilot ETA was ~%.0f min).",
                if (mm < 1) sprintf("%d s", round(mm * 60)) else
                  sprintf("%.1f min", mm),
                if (!is.null(frame$pilot)) frame$pilot$eta_minutes else NA) })
    ## spec 53: report QC - every component that failed or fell back, named,
    ## so an absent section is never mistaken for "not applicable"
    .qc <- .rp_qc_rows()
    md <- c(md, "", "## Report QC", "",
            if (!nrow(.qc))
              "No report component failed or fell back in this run."
            else c(sprintf(paste0("%d report component(s) failed or fell back ",
                                  "(full ledger: `tables/qc_report.csv`). The ",
                                  "analysis itself is unaffected; each line says ",
                                  "what is missing and what took its place."),
                           nrow(.qc)), "",
                   sprintf("- **%s** (%s%s%s): %s. Fallback: %s. Message: %s",
                           .qc$component, .qc$severity,
                           ifelse(is.na(.qc$unit), "", paste0("; ", .qc$unit)),
                           ifelse(is.na(.qc$outcome), "",
                                  paste0("; ", vapply(.qc$outcome, function(.o)
                                    if (is.na(.o)) "" else .lab(frame, .o),
                                    character(1)))),
                           ifelse(.qc$severity == "failed",
                                  "not produced", "degraded"),
                           ifelse(is.na(.qc$fallback), "none", .qc$fallback),
                           .qc$message)))
    ## the reporting block goes LAST, so it reads as the take-away rather than
    ## as another section of diagnostics
    md <- c(md, if (frame$moderation && !is.null(MOD))
                  .rp_mod_block(frame, MOD, nrow(frame$data))
                else .rp_report_block(frame, RES, nrow(frame$data), SHAPE))
    ## summary.md is UTF-8 whatever the session locale: the badge daggers
    ## (\u2020 / \u2021) were written as "<U+2020>" under a C/POSIX locale,
    ## because writeLines() translates to the native encoding first - even
    ## through a UTF-8 connection. useBytes writes the UTF-8 bytes as they are.
    writeLines(enc2utf8(md), smy, useBytes = TRUE)
  }

  ## Say where the output went, every time, whether or not the result is
  ## printed. outdir is relative by default - which is the right default,
  ## because it puts results in the working directory where every other R tool
  ## puts them - but "dmsa_output" is only findable if you already know what it
  ## is relative TO. Resolving it removes the one question a user is otherwise
  ## left holding after a 35-minute run.
  bar$done()
  ## A moderation run draws a panel per SURVIVOR. With no survivor there is
  ## nothing to draw - say that, rather than leaving an empty figures/ folder
  ## looking like a failure.
  if (frame$moderation && isTRUE(frame$plots) && !length(figs))
    message("no figures: a moderation run draws one panel per unit that ",
            "survives its family correction ",
            if (!identical(frame$type, "linear"))
              "on EITHER the linear product or the moderated-curvature arm, "
            else "on the tested product, ",
            "and none did. The unit is collapsed to one subject-level tone ",
            "score before the interaction is fitted, so there is no per-probe ",
            "panel either. Run the same frame with moderation = FALSE for the ",
            "main-effect figures.")
  message(sprintf("DMSA report written to:\n  %s\n  %d figure(s), %d table(s)%s",
                  normalizePath(outdir, winslash = "/", mustWork = FALSE),
                  length(figs), length(tabs),
                  if (frame$summary) ", summary.md" else ""))
  .qc_final <- .rp_qc_rows()
  if (nrow(.qc_final))
    message(sprintf(paste0("  %d report component(s) failed or fell back - ",
                           "see tables/qc_report.csv%s (the analysis itself ",
                           "is unaffected)"),
                    nrow(.qc_final),
                    if (frame$summary) " and the 'Report QC' section of summary.md"
                    else ""))

  structure(list(results = RES, moderation = MOD, shape = SHAPE,
                 probes = probe_res, qc = .rp_qc_rows(),
                 bridge = bridge, figures = figs, tables = tabs,
                 summary = if (frame$summary) smy else NULL,
                 outdir = outdir, frame_args = frame[c(
                   "outcome", "correction", "B", "alpha", "fdr", "cpg_map",
                   "moderation", "mod", "mod2", "frame_role", "type")],
                 minutes = as.numeric(Sys.time() - t_start, units = "mins")),
            class = "dmsa_report")
}

## units per level for the moderation path: named list of families, each a
## named list of units with $columns and $alignment
.rp_units <- function(frame, level) {
  out <- list()
  if (level == "system") {
    fam <- list()
    for (sid in names(frame$sets)) {
      s <- frame$sets[[sid]]
      fam[[s$system]] <- list(columns = s$columns, alignment = s$alignment)
    }
    out[["systems"]] <- fam
  } else if (level == "gene") {
    for (sid in names(frame$sets)) {
      s <- frame$sets[[sid]]; fam <- list()
      for (g in unique(s$map$gene)) {
        k <- s$map$gene == g
        if (!any(is.finite(s$map$best_direction[k]))) next
        al <- dmsa_align(data.frame(cpg = s$map$probe[k],
                                    d = s$map$best_direction[k],
                                    p_plus = s$map$p_plus[k]),
                         genes = s$map$gene[k], level = "gene")
        fam[[g]] <- list(columns = s$map$column[k], alignment = al)
      }
      if (length(fam)) out[[s$system]] <- fam
    }
  } else if (level == "probe") {
    ## Family = the probes of one gene, which is the same level-local family
    ## the main-effect path corrects within. Without this branch the shape scan
    ## and the moderation path silently stopped at the gene level while the
    ## declared hierarchy said system > module > gene > probe.
    for (sid in names(frame$sets)) {
      s <- frame$sets[[sid]]
      for (g in unique(s$map$gene)) {
        k <- which(s$map$gene == g & is.finite(s$map$best_direction))
        if (!length(k)) next
        fam <- list()
        for (i in k) {
          al <- dmsa_align(data.frame(cpg = s$map$probe[i],
                                      d = s$map$best_direction[i],
                                      p_plus = s$map$p_plus[i]),
                           genes = s$map$probe[i], level = "gene")
          fam[[s$map$probe[i]]] <- list(columns = s$map$column[i],
                                        alignment = al)
        }
        if (length(fam)) out[[paste0(s$system, " > ", g)]] <- fam
      }
    }
  } else if (level == "module" && !is.null(frame$modules)) {
    for (sid in intersect(unique(frame$modules$system_id),
                          as.numeric(names(frame$sets)))) {
      s <- frame$sets[[as.character(sid)]]
      mm <- frame$modules[frame$modules$system_id == sid, ]
      md <- mm$module[match(s$map$gene, mm$gene)]
      fam <- list()
      for (m in unique(stats::na.omit(md))) {
        k <- which(md == m)
        al <- dmsa_align(data.frame(cpg = s$map$probe[k],
                                    d = s$map$best_direction[k],
                                    p_plus = s$map$p_plus[k]),
                         genes = s$map$gene[k], level = "gene")
        fam[[m]] <- list(columns = s$map$column[k], alignment = al)
      }
      if (length(fam)) out[[s$system]] <- fam
    }
  }
  out
}

## the system-level frame: all chosen systems' probes, units = system,
## spec 40: what the system-level tone is actually built from. Per analysed
## system: the testable genes, how many carry a signed polarity, how many
## are explicitly off-axis (w_g = 0) and how many are UNRESOLVED (no entry,
## or NA) - those last are weighted 0 by .rp_system_frame() and therefore
## drop out of the system score, which must be said, not implied.
.rp_polarity_audit <- function(frame) {
  mp <- frame$map
  if (is.null(mp) || !nrow(mp)) return(NULL)
  pol <- frame$polarity_table
  if (inherits(pol, "dmsa_polarity")) pol <- pol$polarity
  pol <- if (is.null(pol)) NULL else as.data.frame(pol)
  ok <- !is.null(pol) && all(c("system_id", "gene", "w_g") %in% names(pol))
  g <- unique(mp[, c("system_id", "system", "gene")])
  key <- paste(g$system_id, g$gene)
  w <- if (ok) pol$w_g[match(key, paste(pol$system_id, pol$gene))] else
    rep(NA_real_, nrow(g))
  g$status <- ifelse(is.na(w), "unresolved", ifelse(w == 0, "off_axis", "signed"))
  sp <- split(g, as.character(g$system_id))
  out <- do.call(rbind, lapply(sp, function(d) data.frame(
    system_id = as.character(d$system_id[1]), system = d$system[1],
    n_genes_testable = nrow(d),
    n_polarity_signed = sum(d$status == "signed"),
    n_polarity_off_axis = sum(d$status == "off_axis"),
    n_polarity_unresolved = sum(d$status == "unresolved"),
    genes_unresolved = paste(sort(d$gene[d$status == "unresolved"]),
                             collapse = ";"),
    stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  attr(out, "has_polarity") <- ok
  out
}

## alignment weighted by curated polarity where available
.rp_system_frame <- function(frame) {
  mp <- frame$map
  ## spec 39/40: POLARITY COMES FROM THE ACTIVE REFERENCE, never unconditionally
  ## from Alpha. Calling alpha_polarity() here meant a user analysing their own
  ## reference silently had Project Alpha's gene-to-system signs applied to
  ## THEIR genes wherever the symbols happened to collide - a wrong adjustment
  ## nothing announced. And when nothing matched, every gene was weighted +1,
  ## which is the "silent all +1" substitution the spec forbids: it reports a
  ## purely activating system with no brake, whatever the biology says.
  pol <- frame$polarity_table
  if (inherits(pol, "dmsa_polarity")) pol <- pol$polarity
  pol <- if (is.null(pol)) NULL else as.data.frame(pol)
  ok <- !is.null(pol) && all(c("system_id", "gene", "w_g") %in% names(pol))
  w <- rep(0, nrow(mp))
  if (ok) {
    key <- paste(mp$system_id, mp$gene)
    hit <- match(key, paste(pol$system_id, pol$gene))
    w <- ifelse(is.na(hit), 0, pol$w_g[hit])
  }
  if (!ok || all(w == 0)) {
    ## Signed system analysis is UNAVAILABLE - but the gene and probe levels
    ## are unaffected (they align on d alone and need no w_g), so the report
    ## continues without the system level rather than losing valid analysis.
    attr(mp, "why") <- if (!ok)
      "the active reference carries no polarity table"
      else "no polarity entry in the active reference matches the chosen systems"
    return(structure(list(), class = "dmsa_no_polarity", why = attr(mp, "why")))
  }
  d <- mp$best_direction; pp <- pmin(pmax(mp$p_plus, 0), 1)
  al <- data.frame(probe = mp$probe, gene = mp$gene, d = d, p_plus = pp,
                   w_g = w, s = d * w,
                   p_s_plus = pp * (w > 0) + (1 - pp) * (w < 0),
                   usable = is.finite(d) & w != 0,
                   reason = "", stringsAsFactors = FALSE)
  list(M = frame$M[, mp$column, drop = FALSE], units = mp$system,
       alignment = al)
}

#' @export
print.dmsa_report <- function(x, ...) {
  dur <- if (x$minutes < 1) sprintf("%ds", round(x$minutes * 60)) else
    sprintf("%.1f min", x$minutes)
  cat("DMSA report -", dur, "->",
      normalizePath(x$outdir, winslash = "/", mustWork = FALSE), "\n")
  if (!is.null(x$results)) {
    hits <- x$results[x$results$selected & x$results$level == "gene", ,
                      drop = FALSE]
    cat(sprintf("  gene level: %d unit(s) named (calibrated BH, FDR %.2g)\n",
                nrow(hits), x$frame_args$fdr %||% .10))
    for (i in seq_len(nrow(hits)))
      cat(sprintf("   - %-10s %-24s %s  q %.3f (%s lens)  joint union adj %.4f  omnibus adj %.4f%s%s\n",
                  hits$unit[i], substr(hits$system[i], 1, 24),
                  hits$outcome[i],
                  hits$q_unit_dep[i] %||% NA_real_,
                  hits$best_lens[i],
                  hits$p_unit_adj[i] %||% NA_real_,
                  hits$p_omnibus_adj[i] %||% NA_real_,
                  if (hits$exact_confirmed[i] %in% TRUE)
                    "  EXACT-CONFIRMED" else "",
                  if (hits$omnibus_confirmed[i] %in% TRUE)
                    "  OMNIBUS-CONFIRMED" else ""))
  }
  if (!is.null(x$moderation)) {
    M <- x$moderation
    mh <- M[M$selected, , drop = FALSE]
    tp <- M[which.min(M$p_composite), ]
    cat(sprintf("  moderation (composite): %d unit(s) tested in %d famil%s -> %d survive\n",
                nrow(M), length(unique(M$family)),
                if (length(unique(M$family)) == 1) "y" else "ies", nrow(mh)))
    cat(sprintf("   strongest: %-12s %s level  b %+.3f  t %+.2f  p %.4f  adj %.4f\n",
                tp$unit, tp$level, tp$b, tp$t, tp$p_composite,
                tp$p_composite_adj))
    ## one line per DISTINCT survivor: a probe that also survives as its
    ## gene's or module's only member used to print three near-identical
    ## lines with a cryptic "[= cgX]" tail (PI, 2026-08-29)
    .dk <- paste(mh$unit, mh$outcome)
    for (i in which(!duplicated(.dk))) {
      .nlv <- sum(.dk == .dk[i])
      cat(sprintf("   - %s (%s) adj %.4f%s%s\n", mh$unit[i], mh$outcome[i],
                  mh$p_composite_adj[i],
                  if (isTRUE(mh$n_family[i] == 1L)) "  [family of 1: uncorrected]" else "",
                  if (.nlv > 1L)
                    sprintf("  [same data at %d levels: %s]", .nlv,
                            paste(mh$level[.dk == .dk[i]], collapse = "/"))
                  else ""))
    }
    if (!is.null(M$p_curv) && any(is.finite(M$p_curv))) {
      cq <- M[which.min(M$p_curv), ]
      ch <- M[is.finite(M$p_curv_adj) &
              M$p_curv_adj < (x$frame_args$alpha %||% 0.05), , drop = FALSE]
      cat(sprintf("  moderated non-linearity (S^2 x %s): %d survive%s; strongest %s b %+.3f t %+.2f p %.4f adj %.4f\n",
                  x$frame_args$mod, nrow(ch),
                  { nd <- sum(ch$distinct, na.rm = TRUE)
                    if (nrow(ch) && nd < nrow(ch)) sprintf(" (%d distinct)", nd) else "" },
                  cq$unit, cq$b_curv, cq$t_curv, cq$p_curv, cq$p_curv_adj))
    }
  }
  if (!is.null(x$shape)) {
    S <- x$shape
    cat(sprintf("  shape scan (%s, does NOT change the unit results above):\n",
                x$frame_args$type %||% "non-linear"))
    cat(sprintf("   %d unit(s) scanned across %s; departure terms fitted WITH the linear term,\n",
                nrow(S), paste(unique(S$level), collapse = " > ")))
    cat("   every departure term maxT-corrected inside its level-local family\n")
    .al <- (x$frame_args$alpha %||% 0.05)
    am <- cbind(S$p_quad_adj, S$p_thr_adj, S$p_expo_adj)
    sig <- S[apply(am, 1, function(z) any(is.finite(z) & z < .al)), ,
             drop = FALSE]
    if (!nrow(sig)) cat("   no unit's departure term survives its family\n")
    for (i in seq_len(nrow(sig))) {
      ## the display arm is the FINITE minimum over all three adjusted arms -
      ## NA-safe (a unit can survive on the threshold arm with the quadratic
      ## arm NA, and `if (x < NA)` is an error, not a choice)
      .aa <- c(quad = sig$p_quad_adj[i], thr = sig$p_thr_adj[i],
               expo = sig$p_expo_adj[i])
      .aa[!is.finite(.aa)] <- Inf
      .arm <- names(which.min(.aa))
      .raw <- switch(.arm, quad = sig$p_quad[i], thr = sig$p_thr[i],
                     expo = sig$p_expo[i])
      cat(sprintf("   %-22s %-6s %-13s lin %.4f  +%s %.4f  adj %.4f  %s\n",
                  substr(sig$unit[i], 1, 22), sig$level[i],
                  substr(sig$outcome[i], 1, 13),
                  sig$p_lin[i], .arm,
                  .raw,
                  min(.aa),
                  if (!is.finite(sig$p_ushape[i])) "" else
                  if (sig$p_ushape[i] < 0.05 && isTRUE(sig$turn_inside[i]))
                    sprintf("%s p %.3f", sig$shape[i], sig$p_ushape[i])
                  else sprintf("no %s (p %.3f)", sig$shape[i], sig$p_ushape[i])))
    }
  }
  cat(" ", x$bridge$status, "\n")
  if (length(x$figures)) cat("  figures:", length(x$figures), "\n")
  if (!is.null(x$summary)) cat("  summary:", x$summary, "\n")
  invisible(x)
}

## ---- "how to report this" -------------------------------------------------
## A DMSA result carries three things a reader has never seen together before:
## WHICH LENS carried it, WHICH FAMILY it was corrected inside, and a direction
## that is a claim about EXPRESSION TONE rather than about methylation level.
## Reported wrongly - "p = .012 for AVP methylation" - the result is unreadable
## and, worse, sounds like an ordinary EWAS hit. So the report writes the
## sentences, with this run's own numbers in them, rather than leaving a user to
## reverse-engineer the phrasing from a table.
.rp_tone <- function(d) if (is.na(d)) "an undetermined expression tone" else
  if (d > 0) "a HIGHER expression tone" else "a LOWER expression tone"


## The direction sentence, correct for what the outcome actually is.
## `outcome_label` may be given as TWO strings for a two-level outcome, naming
## the levels (e.g. c("T1", "T4") or c("female", "male")); otherwise the levels
## are named by their coded values so the sentence is still unambiguous.
.rp_level_names <- function(frame, oc) {
  k <- (frame$outcome_kind %||% list())[[oc]]
  ol <- frame$outcome_levels
  ## the NAMED LIST form labels each outcome's levels explicitly and works
  ## for any number of outcomes; the bare length-2 vector is declared once
  ## for the whole frame, so with several outcomes it can only describe one
  ## of them and is honoured only for a single-outcome frame - otherwise the
  ## wave outcome's labels ("T1"/"T4") would title every other outcome's
  ## figure. The result carries attr "custom" so the direction sentence
  ## knows whether to append the coded values in brackets.
  if (is.list(ol) && oc %in% names(ol))
    return(structure(as.character(ol[[oc]])[seq_len(2)], custom = TRUE))
  if (!is.null(ol) && !is.list(ol) && length(ol) >= 2L &&
      length(frame$outcome) == 1L)
    return(structure(as.character(ol)[seq_len(2)], custom = TRUE))
  nm <- .lab(frame, oc)
  if (is.null(k) || is.na(k$lo)) return(c(nm, nm))
  c(paste0(nm, " = ", k$lo), paste0(nm, " = ", k$hi))
}

.rp_dir_sentence <- function(frame, oc, direction, what) {
  k <- (frame$outcome_kind %||% list())[[oc]]
  kind <- if (is.null(k)) "continuous" else k$kind
  tone <- .rp_tone(direction)
  if (identical(kind, "wave")) {
    lv <- .rp_level_names(frame, oc)
    sprintf(paste0("Within person, from %s to %s, the methylation of %s ",
                   "shifted toward %s. This is a change contrast: the ",
                   "permutation is restricted to exchangeable blocks, so both ",
                   "waves of the same participant stay together."),
            lv[1], lv[2], what, tone)
  } else if (identical(kind, "binary")) {
    lv <- .rp_level_names(frame, oc)
    ## with declared labels the coded values ride along in brackets so the
    ## sentence is self-verifying against the data ("taking pills [1] ...
    ## than not taking pills [0]" - PI, 2026-08-29); without labels the
    ## fallback names already carry the values, so no brackets are added
    if (isTRUE(attr(lv, "custom")) && !is.null(k) && !is.na(k$lo)) {
      lv <- c(sprintf("%s [%s]", lv[1], k$lo),
              sprintf("%s [%s]", lv[2], k$hi))
    }
    sprintf(paste0("%s was associated with methylation consistent with %s ",
                   "of %s, relative to %s. This is a group contrast, not a ",
                   "dose-response relationship."),
            sub("^(.)", "\\U\\1", lv[2], perl = TRUE), tone, what, lv[1])
  } else {
    sprintf("Higher %s was associated with methylation consistent with %s of %s.",
            .lab(frame, oc), tone, what)
  }
}

.rp_lens_gloss <- function(l) switch(
  as.character(l),
  coherence = "cross-probe sign agreement - the gene's probes move together",
  composite = "a block shift detected by the averaging statistic",
  diffuse   = "thin signal spread across many probes",
  "an unnamed lens")

## Moderation figure: simple slopes carrying their own partial residuals, and
## Johnson-Neyman. A surviving interaction is the one result a reader cannot
## picture from a coefficient, and a slopes plot drawn without the residuals it
## was fitted on asks them to take the lines on trust. The J-N panel answers the
## question the slopes panel raises - over what range of the moderator is the
## association actually distinguishable from zero - and carries a rug of the
## observed moderator, because a J-N boundary outside the data is not a finding.


## short form for cramped 3D axis labels
.rp_short <- function(x, n = 18L) {
  x <- as.character(x)[1]
  if (is.na(x) || nchar(x) <= n) return(x)
  paste0(substr(x, 1L, n - 3L), "...")
}


## ---- binary outcome: the response curve ----------------------------------
## A binary outcome IS non-linear - probability is bounded at 0 and 1 - so the
## relationship between an aligned tone score and a binary outcome is a logistic
## curve, and a reader is entitled to see it. Through 1.19.0 the shape scan
## simply refused: it puts the OUTCOME on the left, and a 0/1 response fitted by
## least squares is a linear probability model on which the Sasabuchi and
## Fieller machinery is meaningless. That reasoning was right about the LINEAR
## fit and wrong to conclude nothing could be drawn. On the LOGIT scale the
## quadratic term is perfectly well defined, so the curvature is tested with a
## likelihood-ratio test against the linear logistic model, and the response
## curve is drawn on the probability scale with the observed data binned onto it.
.rp_fig_binary <- function(frame, oc, level, family, unit, file) {
  k <- (frame$outcome_kind %||% list())[[oc]]
  if (is.null(k) || !identical(k$kind, "binary") && !identical(k$kind, "wave"))
    return(invisible(NULL))
  ul <- try(.rp_units(frame, level), silent = TRUE)
  if (inherits(ul, "try-error") || is.null(ul[[family]])) return(invisible(NULL))
  u <- ul[[family]][[unit]]
  if (is.null(u)) return(invisible(NULL))
  Sr <- dmsa_scores(frame$M[, u$columns, drop = FALSE], u$alignment,
                    weighting = frame$weighting %||% "combined",
                    w_floor = frame$w_floor %||% 1.5)$aligned
  d <- frame$data
  d$S <- as.numeric(scale(Sr))
  y <- d[[oc]]
  yb <- as.integer(y == k$hi)
  ok <- is.finite(d$S) & is.finite(yb)
  if (sum(ok) < 40L || length(unique(yb[ok])) < 2L) return(invisible(NULL))
  covs <- Filter(nzchar, c(setdiff(frame$outcome, oc), frame$covariates,
                           frame$chip))
  covs <- intersect(covs, names(d))
  dd <- d[ok, , drop = FALSE]; dd$.y <- yb[ok]
  f1 <- stats::as.formula(paste(".y ~", paste(c("S", covs), collapse = " + ")))
  f2 <- stats::as.formula(paste(".y ~", paste(c("S", "I(S^2)", covs), collapse = " + ")))
  m1 <- try(stats::glm(f1, dd, family = stats::binomial()), silent = TRUE)
  m2 <- try(stats::glm(f2, dd, family = stats::binomial()), silent = TRUE)
  if (inherits(m1, "try-error") || inherits(m2, "try-error"))
    return(invisible(NULL))
  lrt <- .rp_try("response curve: likelihood-ratio test",
                 stats::anova(m1, m2, test = "LRT"),
                 fallback = "curve drawn without the LRT line")
  p_curv <- if (is.null(lrt)) NA_real_ else lrt[["Pr(>Chi)"]][2]

  ## predicted probability across the score, covariates held at their means
  nd0 <- dd[rep(1L, 200L), , drop = FALSE]
  for (cv in covs) if (is.numeric(dd[[cv]])) nd0[[cv]] <- mean(dd[[cv]], na.rm = TRUE)
  gs <- seq(min(dd$S), max(dd$S), length.out = 200L)
  nd0$S <- gs
  pr <- function(m) stats::predict(m, nd0, type = "response")
  p1 <- tryCatch(pr(m1), error = function(e) rep(NA_real_, 200L))
  p2 <- tryCatch(pr(m2), error = function(e) rep(NA_real_, 200L))
  qlim <- stats::quantile(dd$S, c(.05, .95), na.rm = TRUE)
  inside <- gs >= qlim[1] & gs <= qlim[2]

  ## observed proportions in equal-count bins, with Wilson intervals
  nb <- max(4L, min(10L, floor(nrow(dd) / 25)))
  br <- stats::quantile(dd$S, seq(0, 1, length.out = nb + 1L), na.rm = TRUE)
  br[1] <- br[1] - 1e-9; br[length(br)] <- br[length(br)] + 1e-9
  cutb <- cut(dd$S, br, include.lowest = TRUE)
  bx <- tapply(dd$S, cutb, mean); bn <- tapply(dd$.y, cutb, length)
  bk <- tapply(dd$.y, cutb, sum); bp <- bk / bn
  z <- 1.959964
  den <- 1 + z^2 / bn
  ctr <- (bp + z^2 / (2 * bn)) / den
  hw  <- z * sqrt(bp * (1 - bp) / bn + z^2 / (4 * bn^2)) / den

  fp <- .rp_dev(frame, file, width = 7.6, height = 5.2)
  on.exit(grDevices::dev.off(), add = TRUE)
  lv <- .rp_level_names(frame, oc)
  ttl <- sprintf("%s (%s) -> P(%s)", .rp_unit_label(frame, unit, level), level,
                 lv[2])
  sub <- sprintf(paste0("logistic response; curvature tested on the LOGIT scale ",
                        "(LRT vs linear, p = %s). solid = central 90%% of the ",
                        "score, dashed = extrapolation"),
                 if (is.finite(p_curv)) sprintf("%.4f", p_curv) else "NA")
  tf <- .rp_fit(ttl, cex = .8); sf <- .rp_fit(sub, cex = .58)
  graphics::par(mar = c(4.4, 4.6, 2.2 + 1.1 * (length(tf$lines) +
                        length(sf$lines)), 1.4))
  graphics::plot(gs, p2, type = "n", ylim = c(0, 1),
                 xlab = sprintf("%s aligned tone score (z)",
                                .rp_unit_label(frame, unit, level)),
                 ylab = sprintf("P(%s)", lv[2]),
                 cex.lab = .82, cex.axis = .75, col.axis = "grey30", bty = "n")
  graphics::abline(h = mean(dd$.y), col = "grey85", lty = 3)
  ## observed, with intervals
  graphics::segments(bx, pmax(0, ctr - hw), bx, pmin(1, ctr + hw),
                     col = "grey55", lwd = 1.6)
  graphics::points(bx, bp, pch = 19, cex = .75, col = "grey30")
  ## the two fits
  for (v in list(list(p = p1, col = "#4C3B73", lty = 1),
                 list(p = p2, col = "#1F6F8B", lty = 1))) {
    graphics::lines(gs[inside], v$p[inside], col = v$col, lwd = 2.4, lty = v$lty)
    graphics::lines(gs[!inside & gs < qlim[1]], v$p[!inside & gs < qlim[1]],
                    col = v$col, lwd = 2.0, lty = 2)
    graphics::lines(gs[!inside & gs > qlim[2]], v$p[!inside & gs > qlim[2]],
                    col = v$col, lwd = 2.0, lty = 2)
  }
  graphics::rug(dd$S, col = grDevices::adjustcolor("grey30", alpha.f = .35))
  graphics::legend("topleft", bty = "o", box.col = NA,
                   bg = grDevices::adjustcolor("white", alpha.f = .82),
                   cex = .6, lwd = 2.2, col = c("#4C3B73", "#1F6F8B", "grey30"),
                   lty = c(1, 1, NA), pch = c(NA, NA, 19),
                   legend = c("logistic, linear in the score",
                              "logistic, + quadratic",
                              sprintf("observed (%d bins, 95%% Wilson)", nb)))
  top <- 0.9 + 0.95 * (length(tf$lines) - 1L) + 0.75 * length(sf$lines)
  for (i in seq_along(tf$lines))
    graphics::mtext(tf$lines[i], side = 3, line = top - 0.95 * (i - 1L),
                    adj = 0, font = 2, cex = tf$cex)
  for (i in seq_along(sf$lines))
    graphics::mtext(sf$lines[i], side = 3,
                    line = top - 0.95 * length(tf$lines) - 0.72 * (i - 1L),
                    adj = 0, cex = sf$cex, col = "grey35")
  invisible(list(file = fp, p_curv = p_curv))
}

## ---- moderated surface ---------------------------------------------------
## A surface needs TWO predictors. The plain shape scan has only the tone score,
## so there is no second axis to put on it and none is drawn. A MODERATED
## non-linear model does have one - tone score x moderator - and the surface
## shows the curvature changing continuously instead of at three sliced levels.
##
## The corners of a quadratic-by-moderator surface extrapolate violently, and a
## 3D surface makes extrapolation look like data in a way a scatterplot never
## does. Grid cells with no observation nearby are drawn grey, and the fraction
## of the surface that is unsupported is printed on the figure.
.rp_fig_surface <- function(frame, row, file) {
  ## drawn only for a moderated NON-LINEAR model: a plain shape scan has one
  ## predictor and therefore no second axis, and a three-way has two moderators
  ## (the sliced-panel figure is the right object there).
  if (identical(frame$type, "linear") || nzchar(frame$mod2 %||% ""))
    return(invisible(NULL))
  ul <- try(.rp_units(frame, row$level), silent = TRUE)
  if (inherits(ul, "try-error") || is.null(ul[[row$family]]))
    return(invisible(NULL))
  u <- ul[[row$family]][[row$unit]]
  if (is.null(u)) return(invisible(NULL))
  oc <- row$outcome
  Sr <- dmsa_scores(frame$M[, u$columns, drop = FALSE], u$alignment,
                    weighting = frame$weighting %||% "combined",
                    w_floor = frame$w_floor %||% 1.5)$aligned
  d <- frame$data
  d$S <- as.numeric(scale(Sr))
  mu <- mean(d[[frame$mod]], na.rm = TRUE)
  sg <- stats::sd(d[[frame$mod]], na.rm = TRUE)
  if (!is.finite(sg) || sg <= 0) return(invisible(NULL))
  d$.mc <- (d[[frame$mod]] - mu) / sg
  d$.sq <- d$S^2 - mean(d$S^2, na.rm = TRUE)
  covs <- Filter(nzchar, c(setdiff(frame$outcome, oc), frame$covariates,
                           frame$chip))
  ff <- stats::as.formula(paste(oc, "~",
          paste(c("S * .mc + .sq * .mc", covs), collapse = " + ")))
  fit <- try(stats::lm(ff, d), silent = TRUE)
  if (inherits(fit, "try-error")) return(invisible(NULL))
  cf <- stats::coef(fit)
  if (!"S" %in% names(cf)) return(invisible(NULL))
  S <- d$S; M <- d$.mc
  ok <- is.finite(S) & is.finite(M)
  S <- S[ok]; M <- M[ok]
  if (length(S) < 20L) return(invisible(NULL))
  ng <- 40L
  gs <- seq(min(S), max(S), length.out = ng)
  gm <- seq(min(M), max(M), length.out = ng)
  msq <- mean(d$S^2, na.rm = TRUE)
  g <- function(n) if (n %in% names(cf)) unname(cf[n]) else 0
  ## the quadratic-by-moderator column is written `.mc:.sq` by R, not `.sq:.mc`
  qn <- intersect(c(".sq:.mc", ".mc:.sq"), names(cf))
  bq <- if (length(qn)) unname(cf[qn[1]]) else 0
  Z <- outer(gs, gm, function(x, m)
    g("(Intercept)") + g("S") * x + g(".sq") * (x^2 - msq) + g(".mc") * m +
    g("S:.mc") * x * m + bq * (x^2 - msq) * m)
  ## support: is there an observation within one bandwidth of this cell?
  hx <- 0.12 * diff(range(S)); hy <- 0.12 * diff(range(M))
  sup <- outer(gs, gm, function(x, m)
    vapply(seq_along(x), function(k)
      any(abs(S - x[k]) <= hx & abs(M - m[k]) <= hy), logical(1)))
  pct_extrap <- 100 * mean(!sup)

  fp <- .rp_dev(frame, file, width = 10.6, height = 5.2)
  on.exit(grDevices::dev.off(), add = TRUE)
  ttl <- sprintf("%s x %s -> %s   (%s level, moderated %s surface)",
                 .rp_unit_label(frame, row$unit, row$level),
                 .lab(frame, frame$mod), .lab(frame, oc), row$level,
                 if (identical(frame$type, "exponential")) "exponential" else
                   "quadratic")
  sub <- sprintf("grey / blank = unsupported: %.0f%% of this surface is extrapolation",
                 pct_extrap)
  tf <- .rp_fit(ttl, cex = .8); sf <- .rp_fit(sub, cex = .6)
  top <- 0.9 + 0.95 * (length(tf$lines) - 1L) + 0.75 * length(sf$lines)
  ## the contour panel needs room for its axis labels; 2.2 lines clipped them
  graphics::par(mfrow = c(1, 2), mar = c(3.8, 3.8, 2.4, 1.2),
                oma = c(0, 0, max(2.8, top + 1.1), 0))

  ## ---- panel 1: 3D perspective ----
  pal <- grDevices::colorRampPalette(c("#3B6EA5", "#F2F2F2", "#B2182B"))(64)
  zf <- (Z[-1, -1] + Z[-1, -ng] + Z[-ng, -1] + Z[-ng, -ng]) / 4
  sf4 <- sup[-1, -1] & sup[-1, -ng] & sup[-ng, -1] & sup[-ng, -ng]
  idx <- cut(zf, breaks = 64, labels = FALSE)
  fcol <- pal[idx]; fcol[!sf4] <- "grey88"
  graphics::persp(gs, gm, Z, col = fcol, theta = -35, phi = 22, expand = .62,
                  border = grDevices::adjustcolor("grey40", alpha.f = .25),
                  ticktype = "detailed", nticks = 4, cex.axis = .6,
                  ## persp draws all three labels along the box edges, where
                  ## two long ones overlap each other. Keep them short here -
                  ## the title carries the full names.
                  xlab = "tone score (z)",
                  ylab = .rp_short(.lab(frame, frame$mod), 18L),
                  zlab = .rp_short(.lab(frame, oc), 18L))
  graphics::mtext("perspective", side = 3, line = .4, cex = .62, col = "grey35")

  ## ---- panel 2: filled contour with the data on it ----
  Zm <- Z; Zm[!sup] <- NA
  graphics::image(gs, gm, Zm, col = pal, useRaster = TRUE,
                  xlab = "", ylab = "", cex.axis = .7, col.axis = "grey30")
  graphics::contour(gs, gm, Zm, add = TRUE, drawlabels = TRUE, labcex = .5,
                    col = grDevices::adjustcolor("grey20", alpha.f = .55))
  graphics::points(S, M, pch = 19, cex = .32,
                   col = grDevices::adjustcolor("grey15", alpha.f = .35))
  graphics::mtext("aligned tone score (z)", side = 1, line = 2.1, cex = .7)
  graphics::mtext(as.character(.lab(frame, frame$mod)), side = 2, line = 2.1,
                  cex = .7)
  graphics::mtext("contour, over the observed data (blank = unsupported)",
                  side = 3, line = .4, cex = .62, col = "grey35")

  for (i in seq_along(tf$lines))
    graphics::mtext(tf$lines[i], side = 3, outer = TRUE,
                    line = top - 0.95 * (i - 1L), font = 2, cex = tf$cex)
  for (i in seq_along(sf$lines))
    graphics::mtext(sf$lines[i], side = 3, outer = TRUE,
                    line = top - 0.95 * length(tf$lines) - 0.75 * (i - 1L),
                    cex = sf$cex, col = "grey35")
  invisible(fp)
}

.rp_fig_moderation <- function(frame, row, file) {
  ## panel captions: wrap to the panel and never spill into the neighbour
  ## Wrapping by a fixed CHARACTER count truncated long labels to "... " even
  ## when they would have fitted; measure against the panel instead, and only
  ## shrink if three lines still will not fit.
  .cap <- function(txt, max_lines = 3L) {
    f <- .rp_fit(txt, cex = .58, width = graphics::par("pin")[1],
                 min_cex = .40, max_lines = max_lines)
    for (k in seq_along(f$lines))
      graphics::mtext(f$lines[k], side = 3,
                      line = 1.05 + (length(f$lines) - k) * .72 -
                             (length(f$lines) - 1) * .72,
                      adj = 0, cex = f$cex, col = "grey30")
  }
  ul <- try(.rp_units(frame, row$level), silent = TRUE)
  if (inherits(ul, "try-error") || is.null(ul[[row$family]])) return(invisible(NULL))
  u <- ul[[row$family]][[row$unit]]
  if (is.null(u)) return(invisible(NULL))
  oc <- row$outcome
  S <- dmsa_scores(frame$M[, u$columns, drop = FALSE], u$alignment,
                   weighting = frame$weighting %||% "combined",
                   w_floor = frame$w_floor %||% 1.5)$aligned
  d <- frame$data
  d$S <- as.numeric(scale(S))
  mu <- mean(d[[frame$mod]], na.rm = TRUE); sg <- stats::sd(d[[frame$mod]], na.rm = TRUE)
  if (!is.finite(sg) || sg <= 0) return(invisible(NULL))
  d$.mc <- (d[[frame$mod]] - mu) / sg
  has2 <- nzchar(frame$mod2)
  if (has2) {
    mu2 <- mean(d[[frame$mod2]], na.rm = TRUE)
    sg2 <- stats::sd(d[[frame$mod2]], na.rm = TRUE)
    if (!is.finite(sg2) || sg2 <= 0) return(invisible(NULL))
    d$.m2c <- (d[[frame$mod2]] - mu2) / sg2
  }
  covs <- Filter(nzchar, c(setdiff(frame$outcome, oc), frame$covariates, frame$chip))
  ## When the run asked for a non-linear moderation, the panel must draw the
  ## CURVED model - otherwise the table reports moderated curvature while the
  ## figure shows straight simple slopes, and the picture contradicts the
  ## analysis it is supposed to illustrate.
  curved <- frame$type != "linear" && !has2
  d$.sq <- d$S^2 - mean(d$S^2, na.rm = TRUE)
  inter <- if (has2) "S * .mc * .m2c" else
    if (curved) "S * .mc + .sq * .mc" else "S * .mc"
  .k_oc <- (frame$outcome_kind %||% list())[[oc]]
  .bin <- !is.null(.k_oc) && identical(.k_oc$kind, "binary") &&
          !is.na(.k_oc$hi)
  ## belt: a frame object built by an OLDER dmsa installation can carry a
  ## stale outcome_kind (the PI's pills battery drew the linear display for
  ## a 0/1 outcome exactly this way) - so the two-level check is repeated
  ## on the data itself
  if (!.bin) {
    .uoc <- if (is.numeric(d[[oc]]))
      sort(unique(d[[oc]][is.finite(d[[oc]])])) else
      sort(unique(as.character(d[[oc]][!is.na(d[[oc]])])))
    if (length(.uoc) == 2L) {
      .bin <- TRUE
      .k_oc <- list(kind = "binary", lo = .uoc[1], hi = .uoc[2])
    }
  }
  if (.bin) {
    ## a two-level outcome gets a LOGISTIC display (PI, 2026-08-29: "the
    ## slope figure should predict likelihood"): same formula, binomial
    ## family, so panel 1 draws probability curves and panel 2 log-odds
    ## slopes. The TESTED statistic is unchanged - the permutation-
    ## calibrated linear product; this fit is the figure's display model.
    d$.oc01 <- as.integer(as.character(d[[oc]]) == as.character(.k_oc$hi))
    ff <- stats::as.formula(paste(".oc01 ~", paste(c(inter, covs),
                                                   collapse = " + ")))
    fit <- try(stats::glm(ff, d, family = stats::binomial()), silent = TRUE)
  } else {
    ff <- stats::as.formula(paste(oc, "~", paste(c(inter, covs),
                                                 collapse = " + ")))
    fit <- try(stats::lm(ff, d), silent = TRUE)
  }
  if (inherits(fit, "try-error")) return(invisible(NULL))
  V <- stats::vcov(fit); cf <- stats::coef(fit)[rownames(V)]
  if (!"S" %in% names(cf)) return(invisible(NULL))

  ## d(outcome)/d(tone score) at a given moderator position, as a contrast over
  ## the coefficient vector - one expression that covers two- and three-way.
  ## name-safe lookup: R writes the quadratic-by-moderator column as `.mc:.sq`
  qn <- { nm <- names(cf); pr <- strsplit(nm, ":", fixed = TRUE)
          hit <- vapply(pr, function(z) length(z) == 2L &&
                          setequal(z, c(".sq", ".mc")), logical(1))
          if (any(hit)) nm[which(hit)[1L]] else NA_character_ }
  ctr <- function(m, m2 = 0, at_s = 0) {
    cc <- stats::setNames(numeric(length(cf)), names(cf))
    cc["S"] <- 1
    if (curved) {
      if (".sq" %in% names(cc)) cc[".sq"] <- 2 * at_s
      if (!is.na(qn)) cc[qn] <- 2 * at_s * m
    }
    if ("S:.mc" %in% names(cc)) cc["S:.mc"] <- m
    if (has2) {
      if ("S:.m2c" %in% names(cc)) cc["S:.m2c"] <- m2
      if ("S:.mc:.m2c" %in% names(cc)) cc["S:.mc:.m2c"] <- m * m2
    }
    cc
  }
  slope <- function(m, m2 = 0, at_s = 0) sum(ctr(m, m2, at_s) * cf)
  sse   <- function(m, m2 = 0, at_s = 0) { cc <- ctr(m, m2, at_s)
                                           sqrt(drop(cc %*% V %*% cc)) }
  ## curvature = d2(y)/dS2 = 2 * (b_sq + b_sq:mc * m); its own contrast
  curvc <- function(m) { cc <- stats::setNames(numeric(length(cf)), names(cf))
    if (".sq" %in% names(cc)) cc[".sq"] <- 2
    if (!is.na(qn)) cc[qn] <- 2 * m
    cc }
  dfr <- stats::df.residual(fit); crit <- stats::qt(.975, dfr)

  mfd <- fit$model
  Sv <- mfd$S; mcv <- mfd$.mc
  m2v <- if (has2) mfd$.m2c else rep(0, length(Sv))
  at2 <- if (has2) c(-1, 1) else 0

  ## Where to draw the simple slopes. The +-1 SD convention puts a line at an
  ## IMPOSSIBLE moderator value whenever the moderator is bounded and skewed -
  ## an adversity count has a floor at zero, and mean - 1 SD is below it. Use
  ## observed quantiles instead, and fall back down a ladder when ties (a
  ## zero-inflated count) collapse two of them onto the same value.
  qlad <- list(c(.10, .50, .90), c(.05, .50, .95), c(.25, .50, .95), c(0, .5, 1))
  qs <- NULL
  for (q in qlad) {
    v <- unname(stats::quantile(mcv, q, na.rm = TRUE))
    if (length(unique(round(v, 8))) == 3L) { qs <- list(q = q, v = v); break }
  }
  if (is.null(qs)) {
    rr1 <- range(mcv, na.rm = TRUE)
    qs <- list(q = c(0, .5, 1),
               v = c(rr1[1], mean(rr1), rr1[2]))
  }
  slev <- qs$v
  slab <- sprintf("%s = %.2f  (%dth pct)", .lab(frame, frame$mod), mu + slev * sg,
                  as.integer(round(qs$q * 100)))
  ## a TWO-LEVEL moderator (sex coded 1/2, a 0/1 group) gets two lines at
  ## its actual levels - the quantile ladder would draw a third line at an
  ## impossible middle value (sex = 1.5) - labelled by mod_levels when
  ## declared (PI, 2026-08-29)
  .u2 <- sort(unique(mcv[is.finite(mcv)]))
  if (length(.u2) == 2L) {
    slev <- .u2
    slab <- if (!is.null(frame$mod_levels))
      sprintf("%s [%.4g]", frame$mod_levels, mu + .u2 * sg)
    else sprintf("%s = %.4g", .lab(frame, frame$mod), mu + .u2 * sg)
  }

  pal <- .rp_pal(frame, 6); LO <- pal[1]; MI <- pal[3]; HI <- pal[5]
  .lcols <- if (length(slev) == 2L) c(LO, HI) else c(LO, MI, HI)
  ## a TWO-LEVEL moderator draws ONE panel (PI ruling below), so the device
  ## is narrowed accordingly instead of stretching one panel to full width
  .mod2lvl <- length(unique(mcv[is.finite(mcv)])) == 2L && !curved
  fp <- .rp_dev(frame, file, width = if (.mod2lvl) 6.8 else 10.5,
                height = 4.6 * length(at2))
  on.exit(grDevices::dev.off(), add = TRUE)
  pq  <- if (is.null(row$p_curv_adj)) NA_real_ else row$p_curv_adj
  lin_ok  <- is.finite(row$p_composite_adj) && row$p_composite_adj < frame$alpha
  curv_ok <- curved && is.finite(pq) && pq < frame$alpha
  parm <- if (curv_ok && !lin_ok)
      sprintf("moderated curvature adj p = %.4f   |   linear product adj p = %.4f",
              pq, row$p_composite_adj)
    else if (curved && is.finite(pq))
      sprintf("family-adjusted p = %.4f   |   moderated curvature adj p = %.4f",
              row$p_composite_adj, pq)
    else sprintf("family-adjusted p = %.4f", row$p_composite_adj)
  ## Fit the title BEFORE par(mfrow) is set. Two reasons: par("cex") is still 1
  ## here, which is the scale outer text is drawn at (inside a 2x2 layout it
  ## drops to 0.83 and the measurement comes out 1.2x too small); and the outer
  ## margin has to be sized for the number of lines the title actually needs -
  ## a fixed 3.2 lines silently clipped the first line of a 3-line title off the
  ## top of the device.
  .ttl0 <- sprintf("%s x %s%s -> %s   (%s level)",
                   .rp_unit_label(frame, row$unit, row$level),
                   .lab(frame, frame$mod),
                   if (has2) paste0(" x ", .lab(frame, frame$mod2)) else "",
                   .lab(frame, oc), row$level)
  .sub0 <- parm
  if (has2)
    .sub0 <- paste0(.sub0, "   |   rows: ", .lab(frame, frame$mod2),
                    " held at -1 SD (top) and +1 SD (bottom)")
  .tf <- .rp_fit(.ttl0, cex = .8)
  .sf <- .rp_fit(.sub0, cex = .62)
  .top <- 0.9 + 0.95 * (length(.tf$lines) - 1L) + 0.75 * length(.sf$lines)
  ## a TWO-LEVEL moderator gets NO second panel at all (PI, 2026-08-29): a
  ## significant interaction means, by definition, that the two slopes
  ## differ - each level's slope and its significance are stated in the
  ## legend of the one panel instead
  graphics::par(mfrow = c(length(at2), if (.mod2lvl) 1L else 2L),
                mar = c(4.4, 4.6, 4.2, 1.4),
                oma = c(0, 0, max(2.6, .top + 1.1), 0))

  for (a2 in at2) {
    ## ---- panel 1: simple slopes over the partial residuals ----------------
    if (.bin) pres <- NULL else if (curved) {
      bS_i <- unname(cf["S"]) + (if ("S:.mc" %in% names(cf))
        unname(cf["S:.mc"]) * mcv else 0)
      bQ_i <- (if (".sq" %in% names(cf)) unname(cf[".sq"]) else 0) +
              (if (!is.na(qn)) unname(cf[qn]) * mcv else 0)
      pres <- as.numeric(stats::residuals(fit)) + bS_i * Sv +
              bQ_i * (Sv^2 - mean(Sv^2, na.rm = TRUE))
    } else {
      sl_i <- vapply(seq_along(Sv), function(i)
        slope(mcv[i], if (has2) m2v[i] else 0), numeric(1))
      pres <- as.numeric(stats::residuals(fit)) + sl_i * Sv
    }
    ## Colour the points on a continuous ramp rather than by tertile: a
    ## moderator like an adversity count is heavily tied at zero, and its
    ## lower quantiles collapse onto the same break.
    rr0 <- range(mcv, na.rm = TRUE)
    tt <- if (diff(rr0) > 0) (mcv - rr0[1]) / diff(rr0) else rep(.5, length(mcv))
    ramp <- grDevices::colorRamp(c(LO, MI, HI))
    cols <- grDevices::rgb(ramp(tt), maxColorValue = 255)
    ## Reserve headroom for the legend. Drawn at "topleft" over a full-height
    ## y range it sat on top of the scatter and the slope lines, and the moderator
    ## values in the legend text were unreadable where they crossed a point.
    if (.bin) {
      ## probability display: observed 0/1 jittered, fitted P(high level)
      ## per moderator level, covariates absorbed into an average intercept
      lp <- as.numeric(stats::predict(fit))          # linear predictor
      msq1 <- mean(Sv^2, na.rm = TRUE)
      g1 <- function(n) if (n %in% names(cf)) unname(cf[n]) else 0
      .contrib <- function(x, m, m2) {
        e <- g1("S") * x + g1(".mc") * m + g1("S:.mc") * x * m
        if (has2) e <- e + g1(".m2c") * m2 + g1("S:.m2c") * x * m2 +
          g1(".mc:.m2c") * m * m2 + g1("S:.mc:.m2c") * x * m * m2
        if (curved) e <- e + g1(".sq") * (x^2 - msq1) +
          (if (!is.na(qn)) unname(cf[qn]) else 0) * (x^2 - msq1) * m
        e
      }
      eta0 <- mean(lp - .contrib(Sv, mcv, if (has2) m2v else 0),
                   na.rm = TRUE)
      y01 <- mfd$.oc01
      .hilab <- .rp_level_names(frame, oc)[2]
      graphics::plot(Sv, jitter(y01, amount = .035), pch = 19, cex = .5,
                     col = grDevices::adjustcolor(cols, alpha.f = .55),
                     ylim = c(-.08, 1.28),
                     xlab = sprintf("%s aligned tone score (z)",
                        .rp_unit_label(frame, row$unit, row$level)),
                     ylab = sprintf("P(%s)", .hilab),
                     cex.lab = .8, cex.axis = .75, col.axis = "grey30",
                     bty = "n")
      graphics::abline(h = c(0, 1), col = "grey85")
      gs <- seq(min(Sv), max(Sv), length.out = 200)
      for (k in seq_along(slev))
        graphics::lines(gs, stats::plogis(eta0 + .contrib(gs, slev[k], a2)),
                        col = .lcols[k], lwd = 2.2)
    } else {
    .yr <- range(pres, na.rm = TRUE)
    .ylim <- c(.yr[1], .yr[2] + 0.26 * diff(.yr))
    graphics::plot(Sv, pres, pch = 19, cex = .5,
                   col = grDevices::adjustcolor(cols, alpha.f = .55),
                   ylim = .ylim,
                   xlab = sprintf("%s aligned tone score (z)",
                      .rp_unit_label(frame, row$unit, row$level)),
                   ylab = sprintf("%s | partial residual", .lab(frame, oc)),
                   cex.lab = .8, cex.axis = .75, col.axis = "grey30",
                   bty = "n")
    graphics::abline(h = 0, col = "grey85")
    if (curved) {
      ## partial component for S AND S^2, so the drawn curve is the fitted one
      gs <- seq(min(Sv), max(Sv), length.out = 200)
      for (k in seq_along(slev)) {
        m <- slev[k]
        bS <- unname(cf["S"]) + (if ("S:.mc" %in% names(cf))
          unname(cf["S:.mc"]) * m else 0)
        bQ <- (if (".sq" %in% names(cf)) unname(cf[".sq"]) else 0) +
              (if (!is.na(qn)) unname(cf[qn]) * m else 0)
        graphics::lines(gs, bS * gs + bQ * (gs^2 - mean(Sv^2, na.rm = TRUE)),
                        col = .lcols[k], lwd = 2.2)
      }
    } else
    for (k in seq_along(slev))
      graphics::abline(a = 0, b = slope(slev[k], a2), col = .lcols[k],
                       lwd = 2.2)
    }
    ## opaque backing so the legend stays readable even if a point strays under it
    .leg <- slab
    if (.mod2lvl && length(slev) == 2L) {
      .sl2 <- vapply(slev, slope, numeric(1), m2 = a2)
      .se2 <- vapply(slev, sse, numeric(1), m2 = a2)
      .pp2 <- 2 * stats::pt(-abs(.sl2 / .se2), df = dfr)
      .leg <- sprintf("%s:  b %+.3f, p %s%s", slab, .sl2,
                      ifelse(.pp2 < .001, "< .001",
                             sub("^0", "", sprintf("%.3f", .pp2))),
                      if (.bin) " (log-odds)" else "")
    }
    graphics::legend("topleft", bty = "o", box.col = NA,
                     bg = grDevices::adjustcolor("white", alpha.f = .82),
                     cex = .62, lwd = 2.2, col = .lcols, legend = .leg)
    ## In a three-way panel the row label lived only in the outer note, so a
    ## reader had to remember it and count rows. Name the level in the panel.
    rowlab <- if (has2)
      sprintf("%s = %.2f (%s SD).  ", .lab(frame, frame$mod2), mu2 + a2 * sg2,
              if (a2 < 0) "-1" else "+1") else ""
    ## Short, and WRAPPED to the panel. The full model statement belongs in the
    ## summary; four unwrapped lines here collided with the legend, and one long
    ## line ran into the next panel's caption in the three-way layout.
    .cap(paste0(rowlab,
                if (.bin) sprintf(paste0("Fitted probability of %s per level ",
                    "of %s (logistic display; the tested statistic is the ",
                    "permutation-calibrated linear product)"),
                    .rp_level_names(frame, oc)[2], .lab(frame, frame$mod))
                else paste0(if (curved) "Simple curves" else "Simple slopes",
                            " at ", length(slev), " levels of ",
                            .lab(frame, frame$mod),
                            ", over partial residuals"),
                if (.mod2lvl) paste0(". The tested interaction IS the ",
                    "difference between the two slopes; each slope\'s own ",
                    "b and p are in the legend.") else ""))

    ## ---- panel 2: Johnson-Neyman -----------------------------------------
    ## no second panel for a TWO-LEVEL moderator (PI ruling): the per-level
    ## slopes and their p-values live in the first panel\'s legend
    if (.mod2lvl) next
    rng <- range(mcv, na.rm = TRUE)
    gz <- seq(rng[1], rng[2], length.out = 400)
    if (curved) {
      ## J-N is about a slope, and with curvature the slope depends on the score
      ## as well as the moderator - so the honest analogue is the CURVATURE as a
      ## function of the moderator: where is the relationship bent at all?
      sl <- vapply(gz, function(z) sum(curvc(z) * cf), numeric(1))
      se <- vapply(gz, function(z) { cc <- curvc(z)
                                     sqrt(drop(cc %*% V %*% cc)) }, numeric(1))
    } else {
    sl <- vapply(gz, slope, numeric(1), m2 = a2)
    se <- vapply(gz, sse,   numeric(1), m2 = a2)
    }
    lo <- sl - crit * se; hi <- sl + crit * se
    graw <- mu + gz * sg
    graphics::plot(graw, sl, type = "n", ylim = range(lo, hi, 0),
                   xlab = sprintf("%s (observed range)", .lab(frame, frame$mod)),
                   ylab = if (curved) "curvature of the tone-score effect" else "slope of the tone score", cex.lab = .8,
                   cex.axis = .75, col.axis = "grey30", bty = "n")
    sig <- abs(sl / se) >= crit
    ## shade only where the slope is distinguishable from zero AND inside data
    rr <- rle(sig)
    e <- cumsum(rr$lengths); s <- e - rr$lengths + 1L
    for (i in which(rr$values))
      graphics::rect(graw[s[i]], graphics::par("usr")[3], graw[e[i]],
                     graphics::par("usr")[4],
                     col = grDevices::adjustcolor(HI, alpha.f = .10), border = NA)
    graphics::polygon(c(graw, rev(graw)), c(lo, rev(hi)),
                      col = grDevices::adjustcolor(MI, alpha.f = .18), border = NA)
    graphics::lines(graw, sl, col = MI, lwd = 2.2)
    graphics::abline(h = 0, lty = 2, col = "grey45")
    ## Data-density strip, drawn explicitly rather than with rug(): rug() puts
    ## its ticks on the axis line, where they are invisible against it.
    usr <- graphics::par("usr")
    graphics::segments(mu + mcv * sg, usr[3],
                       mu + mcv * sg, usr[3] + .045 * diff(usr[3:4]),
                       col = grDevices::adjustcolor("grey15", .40), lwd = 1.1)
    bnd <- graw[which(diff(sign(abs(sl / se) - crit)) != 0)]
    for (b in bnd) graphics::abline(v = b, lty = 3, col = "firebrick", lwd = 1.6)
    ## How many people actually sit in the significant region. A J-N boundary
    ## near the top of a skewed moderator can be true and still describe almost
    ## nobody, and that is not visible from the boundary value alone.
    mraw <- mu + mcv * sg
    nsig <- sum(vapply(mraw, function(x) {
      i <- which.min(abs(graw - x)); isTRUE(sig[i]) }, logical(1)))
    jn <- if (!length(bnd))
        sprintf("%s %s across the whole observed range",
                if (curved) "Curvature" else "Slope (Johnson-Neyman)",
                if (all(sig)) "differs from zero" else
                  "is not distinguishable from zero")
      else sprintf(paste0("Shaded: ", if (curved) "curvature" else "slope",
                          " differs from zero (boundary %s); %d/%d (%.0f%%)"),
                   paste(sprintf("%.2f", bnd), collapse = ", "),
                   nsig, length(mraw), 100 * nsig / length(mraw))
    ## no rowlab here: the outer subtitle and the left panel both name the row
    .cap(jn)
  }
  ## Lead with the arm that actually carried this unit. A panel drawn because
  ## the CURVATURE survived used to be headed "family-adjusted p = 1.0000" (the
  ## linear product), which reads as a null result for a figure of a finding.
  ## The p-values used to be appended to the title, which is how a module-level
  ## 3-way title reached 202 characters. The title identifies the model; the
  ## statistics belong on their own line.
  for (i in seq_along(.tf$lines))
    graphics::mtext(.tf$lines[i], side = 3, outer = TRUE,
                    line = .top - 0.95 * (i - 1L), font = 2, cex = .tf$cex)
  for (i in seq_along(.sf$lines))
    graphics::mtext(.sf$lines[i], side = 3, outer = TRUE,
                    line = .top - 0.95 * length(.tf$lines) - 0.75 * (i - 1L),
                    cex = .sf$cex, col = "grey35")
  invisible(fp)
}

## What was actually in the model, in words. A reader cannot audit a moderation
## from "all lower-order terms present": they need the terms enumerated, the
## covariates named, and - the one that is easy to get wrong here - to be told
## that cID enters as a PERMUTATION BLOCK, not as a random intercept, and that
## with several outcomes each is adjusted for the others.
.rp_mod_model_words <- function(frame, oc) {
  has2 <- nzchar(frame$mod2)
  multi <- length(frame$outcome) > 1L
  ## With several outcomes the model is refitted per outcome, so name the
  ## outcome generically rather than pinning the sentence to whichever one
  ## happened to carry the first surviving unit.
  ocw <- if (multi) "each outcome" else oc
  lhs <- if (frame$frame_role == "outcome") "the tone score" else ocw
  sc <- if (frame$frame_role == "outcome") ocw else "the tone score"
  m1 <- frame$mod; m2 <- frame$mod2
  terms <- if (has2)
    c(sc, m1, m2, paste(sc, "x", m1), paste(sc, "x", m2), paste(m1, "x", m2),
      paste(sc, "x", m1, "x", m2))
  else c(sc, m1, paste(sc, "x", m1))
  other <- setdiff(frame$outcome, oc)
  covs <- Filter(nzchar, c(frame$covariates, frame$chip))
  list(
    model = sprintf(
      paste0("The fitted model was %s ~ %s, plus covariates: %s are present, ",
             "so the tested %s is the only term not already accounted for."),
      lhs, paste(terms, collapse = " + "),
      if (has2)
        "all three main effects and all three two-way interactions"
      else "both main effects",
      if (has2) "three-way product" else "two-way product"),
    covs = sprintf(
      paste0("Covariates: %s.%s Cluster structure (%s) was handled by ",
             "restricting the permutation to exchangeable blocks rather than ",
             "by a random intercept - no mixed model was fitted, so there is ",
             "no variance component to report."),
      if (length(covs)) paste(covs, collapse = ", ") else "none",
      if (multi) sprintf(
        paste0(" Because %d outcomes were named (%s), each was additionally ",
               "adjusted for the other%s - the model is refitted per outcome."),
        length(frame$outcome),
        paste(vapply(frame$outcome, function(.o) .lab(frame, .o), character(1)),
              collapse = ", "),
        if (length(frame$outcome) > 2L) "s" else "") else "",
      if (length(frame$block_cols)) paste(frame$block_cols, collapse = ", ")
      else "none"))
}

## The moderated branch. A reader meeting mDMSA has to be told three things
## the main-effect wording never has to say: that the unit was collapsed into
## ONE subject-level tone score before the interaction was fitted, that only
## the composite lens is reported because coherence and diffuse do not carry
## through a product, and what a slope-on-a-slope actually means in words.
.rp_mod_block <- function(frame, MOD, n_used) {
  hits <- MOD[MOD$selected, , drop = FALSE]
  mw <- .rp_mod_model_words(frame, MOD$outcome[1])
  term <- paste0("tone score x ", frame$mod,
                 if (nzchar(frame$mod2)) paste0(" x ", frame$mod2) else "")
  out <- c("", "## How to report this", "",
    paste("Copy and adapt. A moderated DMSA result is not an interaction",
          "between a CpG and a moderator. The set is first collapsed into one",
          "**subject-level aligned tone score** - positive means the set reads",
          "as higher activation tone - and the interaction is fitted on that",
          "score. Say so, or a reader will assume probe-by-moderator testing",
          "and ask why there is no probe-level correction."), "")
  if (!nrow(hits)) {
    tp <- MOD[which.min(MOD$p_composite), ]
    out <- c(out, "### Results (nothing survived)", "",
      sprintf(paste0("> The %s interaction was tested at %d unit(s) across %d ",
        "level-local famil%s (composite lens, %s, B = %s). No unit survived ",
        "its family correction. The strongest interaction was %s at the %s ",
        "level (b = %+.3f, t = %+.2f, raw p = %.4f, family-adjusted p = %.4f). ",
        "This is a calibrated null rather than a failed run: %d participants ",
        "and %d direction-called probe%s were available to the strongest unit."),
        term, nrow(MOD), length(unique(MOD$family)),
        if (length(unique(MOD$family)) == 1) "y" else "ies",
        frame$correction, format(frame$B, big.mark = ","),
        tp$unit, tp$level, tp$b, tp$t, tp$p_composite, tp$p_composite_adj,
        n_used, tp$n_probes, if (isTRUE(tp$n_probes == 1L)) "" else "s"), "")
  } else {
    out <- c(out, "### Results", "")
    for (i in seq_len(nrow(hits))) {
      h <- hits[i, ]
      fam <- sum(MOD$family == h$family & MOD$level == h$level &
                   MOD$outcome == h$outcome, na.rm = TRUE)
      dir <- if (h$b > 0) "more positive" else "more negative"
      pair <- if (frame$frame_role == "outcome")
        sprintf("the association between %s and the %s tone score",
                .lab(frame, h$outcome), h$unit)
      else
        sprintf("the association between the %s tone score and %s",
                h$unit, .lab(frame, h$outcome))
      ## With two moderators the coefficient is the THREE-way term: it is the
      ## change in the mod-moderation per SD of mod2, not the change in the
      ## simple slope. Saying the latter would misreport what was estimated.
      sent <- if (nzchar(frame$mod2))
        sprintf(paste0("each standard deviation of %s made the %s moderation ",
                       "of %s %.3f %s - that is, how much %s moderates the ",
                       "association itself depends on %s"),
                frame$mod2, frame$mod, pair, abs(h$b),
                if (h$b > 0) "larger" else "smaller", frame$mod, frame$mod2)
      else
        sprintf(paste0("each standard deviation of %s changed that slope by ",
                       "%+.3f, so %s became %s as %s rose"),
                frame$mod, h$b, pair, dir, frame$mod)
      out <- c(out, sprintf(
        paste0("> At the %s level, the %s interaction survived correction ",
               "within its %d-unit level-local family (composite lens, %s, ",
               "B = %s, %s engine): **%s**, b = %+.3f, t = %+.2f, ",
               "family-adjusted p = %.4f, on %d direction-called CpG%s.%s ",
               "In words: %s."),
        h$level, term, fam, frame$correction, format(frame$B, big.mark = ","),
        frame$weighting %||% "combined", h$unit, h$b, h$t, h$p_composite_adj,
        h$n_probes, if (isTRUE(h$n_probes == 1L)) "" else "s",
        if (isTRUE(h$n_family == 1L))
          " **That family has only one member, so no correction was possible and this p is the raw permutation p.**"
        else "", sent), "")
    }
  }
  if (!is.null(MOD$p_curv) && any(is.finite(MOD$p_curv)))
    out <- c(out, sprintf(paste0(
      "> Whether the moderation is itself non-linear was tested separately: ",
      "the squared tone score by %s product was entered with S, S^2, %s and ",
      "S x %s all present, and corrected by maxT within the same level-local ",
      "family. %s"),
      .lab(frame, frame$mod), .lab(frame, frame$mod), .lab(frame, frame$mod),
      { cq <- MOD[is.finite(MOD$p_curv_adj) & MOD$p_curv_adj < frame$alpha, ,
                  drop = FALSE]
        if (nrow(cq)) paste0(nrow(cq), " unit(s) survived",
          { nd <- sum(cq$distinct, na.rm = TRUE)
            if (nd < nrow(cq)) sprintf(
              " (%d distinct; the remainder are the same data reported at another level)",
              nd) else "" }, ": ",
          paste(sprintf("%s (%s, b = %+.3f, adjusted p = %.4f%s)", cq$unit,
                        .lab(frame, cq$outcome), cq$b_curv, cq$p_curv_adj,
                        ifelse(!is.na(cq$same_as), paste0(", = ", cq$same_as), "")),
                collapse = "; "),
          ". A positive coefficient means the tone-score curve becomes more ",
          "convex (U-like) as the moderator rises; a negative one means more ",
          "concave (inverted-U-like).")
        else "No unit survived, so the moderation is reported as linear in the tone score." }), "")
  out <- c(out, "### Methods", "",
    sprintf(paste0("> Moderated Directional Methylation Set Analysis (mDMSA; ",
      "dmsa R package v%s) was used. Probe effects were aligned to their ",
      "predicted consequence for gene expression through m_j = d_j x w_g x ",
      "(2p+ - 1), and each unit's aligned probes were collapsed into a single ",
      "standardised subject-level tone score under the %s weighting engine. ",
      "%s was standardised before entering the model. %s %s ",
      "The tested term was the highest-order product (%s).%s ",
      "Only the composite lens is reported: coherence and ",
      "diffuse summarise agreement and spread across probes and do not carry ",
      "through a product, so an interaction has no counterpart under them. ",
      "Significance was obtained by block permutation (B = %s, blocks = %s) ",
      "and family-wise error controlled by Westfall-Young step-down %s within ",
      "level-local families on a shared permutation stream. Analyses used ",
      "n = %d and the %s direction map."),
      utils::packageVersion("dmsa"), frame$weighting %||% "combined",
      paste(c(frame$mod, if (nzchar(frame$mod2)) frame$mod2),
            collapse = " and "),
      mw$model, mw$covs, term,
      if (!is.null(MOD$p_curv) && any(is.finite(MOD$p_curv)))
        sprintf(paste0(" Because type = \"%s\" was requested, a second model ",
          "was fitted per unit adding a centred quadratic tone score and its ",
          "product with %s (S, S^2, %s, S x %s, S^2 x %s), and the ",
          "highest-order term S^2 x %s was tested and corrected the same way; ",
          "that is the moderated non-linearity reported above."),
          frame$type, .lab(frame, frame$mod), .lab(frame, frame$mod),
          .lab(frame, frame$mod), .lab(frame, frame$mod), .lab(frame, frame$mod))
      else "",
      format(frame$B, big.mark = ","),
      if (length(frame$block_cols)) paste(frame$block_cols, collapse = ", ")
      else "none",
      frame$correction, n_used, frame$cpg_map), "",
    paste("The tone score is a claim about expression, not about methylation",
          "level: a positive score means the set's methylation pattern is",
          "consistent with higher expression of its member genes, whichever",
          "way the individual probes moved."), "")
  c(out, .rp_column_glossary(frame))
}

## The shape scan's own Results and Methods prose. Without this, a run with
## type = "non-linear" hands the user a Methods paragraph describing a purely
## linear analysis - which is what they would then publish.
.rp_shape_report <- function(frame, SHAPE, n_used) {
  if (is.null(SHAPE) || !nrow(SHAPE)) return(character(0))
  am <- cbind(SHAPE$p_quad_adj, SHAPE$p_thr_adj, SHAPE$p_expo_adj)
  hit <- which(apply(am, 1, function(z) any(is.finite(z) & z < frame$alpha)))
  Bper <- max(499L, frame$B %/% 4L)
  out <- c("", "### Results - shape", "")
  if (!length(hit)) {
    out <- c(out, sprintf(paste0(
      "> Departure from linearity was screened at %d unit(s) across %s. No ",
      "unit's departure term survived correction within its level-local ",
      "family, so the associations are reported as linear."),
      nrow(SHAPE), paste(unique(SHAPE$level), collapse = " > ")), "")
  } else {
    for (i in hit) {
      r <- SHAPE[i, ]
      lab <- if (is.finite(r$p_expo_adj) && r$p_expo_adj < frame$alpha &&
                 (!is.finite(r$p_quad_adj) || r$p_expo_adj < r$p_quad_adj))
        list(nm = "exponential", p = r$p_expo, adj = r$p_expo_adj)
      else if (is.finite(r$p_thr_adj) && r$p_thr_adj < frame$alpha &&
               (!is.finite(r$p_quad_adj) || r$p_thr_adj < r$p_quad_adj))
        list(nm = "threshold", p = r$p_thr, adj = r$p_thr_adj)
      else list(nm = "quadratic", p = r$p_quad, adj = r$p_quad_adj)
      ## The shape claim, stated at exactly the strength the tests support.
      shp <- if (!is.finite(r$p_ushape)) "" else
        if (r$p_ushape < frame$alpha && isTRUE(r$turn_inside) &&
            (!is.finite(r$p_ushape_trim) || r$p_ushape_trim < frame$alpha) &&
            isTRUE(r$turn_ci_inside))
          sprintf(paste0(" The slope took opposite signs at the two ends of ",
            "the observed range (Sasabuchi p = %.4f), and the extremum at ",
            "%.2f had a Fieller 95%% interval [%.2f, %.2f] inside the data, ",
            "so the relationship is reported as %s-shaped."),
            r$p_ushape, r$turn, r$turn_lo, r$turn_hi, r$shape)
        else if (r$p_ushape < frame$alpha && isTRUE(r$turn_inside))
          sprintf(paste0(" The slope-sign test was significant over the full ",
            "observed range (Sasabuchi p = %.4f)%s%s Because %s, this is ",
            "reported as a curvilinear component rather than an established ",
            "%s shape."),
            r$p_ushape,
            if (is.finite(r$p_ushape_trim))
              sprintf(", but p = %.4f on the central 90%% of the score.",
                      r$p_ushape_trim) else ".",
            if (is.finite(r$turn))
              sprintf(" The extremum was %.2f (Fieller 95%% %s).", r$turn,
                      if (is.finite(r$turn_lo) && is.finite(r$turn_hi))
                        sprintf("[%.2f, %.2f]", r$turn_lo, r$turn_hi)
                      else "unbounded") else "",
            paste(c(if (is.finite(r$p_ushape_trim) &&
                        r$p_ushape_trim >= frame$alpha)
                      "the shape does not hold on the central 90% of the score",
                    if (!isTRUE(r$turn_ci_inside))
                      "the Fieller interval for the extremum is not contained in the observed range"),
                  collapse = " and "),
            r$shape)
        else sprintf(paste0(" The slope did not take opposite signs at the ",
            "two ends of the observed range (Sasabuchi p = %.4f), so no U or ",
            "inverted-U is claimed; the departure is reported as a ",
            "curvilinear component only."), r$p_ushape)
      out <- c(out, sprintf(paste0(
        "> At the %s level, the aligned tone score of **%s** showed a ",
        "departure from linearity in its association with %s that the linear ",
        "term does not capture: the linear term gave p = %.4f, while the %s ",
        "term entered alongside it gave p = %.4f (Westfall-Young maxT within ",
        "the %s-level family of %d unit(s), adjusted p = %.4f; B = %s ",
        "permutations per arm, n = %d).%s"),
        r$level, r$unit, .lab(frame, r$outcome), r$p_lin, lab$nm, lab$p, r$level,
        sum(SHAPE$family == r$family & SHAPE$level == r$level &
              SHAPE$outcome == r$outcome, na.rm = TRUE),
        lab$adj, format(Bper, big.mark = ","), n_used, shp), "")
    }
  }
  c(out, "### Methods - shape", "",
    sprintf(paste0(
      "> Departure from linearity was assessed with a shape scan on each ",
      "unit's aligned tone score, at every declared level (%s). A linear, a ",
      "centred quadratic and a threshold term%s were each entered **alongside ",
      "the linear term**, so every reported p is the 1-df incremental test of ",
      "that departure over and above the linear component; entering a squared ",
      "term without its linear counterpart would impose a turning point at the ",
      "centring value and is not evidence of curvature (marginality principle; ",
      "Morris et al., 2023). Each departure term was tested by block ",
      "permutation (B = %s, blocks = %s) and corrected by Westfall-Young ",
      "step-down maxT within its level-local family. Because a monotone but ",
      "convex relationship also produces a significant quadratic term ",
      "(Simonsohn, 2018), a U or inverted-U shape was claimed only when the ",
      "fitted slope took opposite signs at the two ends of the observed score ",
      "range by the Sasabuchi intersection-union test (Sasabuchi, 1980) as ",
      "applied by Lind and Mehlum (2010), and when the Fieller (1943) interval ",
      "for the extremum fell inside that range. The slope-sign and extremum ",
      "intervals use parametric standard errors rather than the block ",
      "permutation and are therefore anticonservative under clustering. ",
      "Curvilinear terms can also be induced by an omitted product term ",
      "(Ganzach, 1997)."),
      paste(unique(SHAPE$level), collapse = " > "),
      if (any(is.finite(SHAPE$p_expo))) ", and an exponential term" else "",
      format(Bper, big.mark = ","),
      if (length(frame$block_cols))
        paste(frame$block_cols, collapse = ", ") else "none"), "",
    "### References for the shape tests", "",
    paste0("- Fieller, E. C. (1943). Fundamental formula in the statistics of ",
           "biological assay, and some applications. *Quarterly Journal of ",
           "Pharmacy and Pharmacology*, 17, 117-123."),
    paste0("- Ganzach, Y. (1997). Misleading interaction and curvilinear ",
           "terms. *Psychological Methods*, 2(3), 235-247. ",
           "doi:10.1037/1082-989X.2.3.235"),
    paste0("- Lind, J. T., & Mehlum, H. (2010). With or without U? The ",
           "appropriate test for a U-shaped relationship. *Oxford Bulletin of ",
           "Economics and Statistics*, 72(1), 109-118. ",
           "doi:10.1111/j.1468-0084.2009.00569.x"),
    paste0("- Morris, T. P., et al. (2023). The marginality principle ",
           "revisited: should \"higher-order\" terms always be accompanied by ",
           "\"lower-order\" terms in regression analyses? *Biometrical ",
           "Journal*. doi:10.1002/bimj.202300069"),
    paste0("- Sasabuchi, S. (1980). A test of a multivariate normal mean with ",
           "composite hypotheses determined by linear inequalities. ",
           "*Biometrika*, 67(2), 429-439. doi:10.1093/biomet/67.2.429"),
    paste0("- Simonsohn, U. (2018). Two lines: a valid alternative to the ",
           "invalid testing of U-shaped relationships with quadratic ",
           "regressions. *Advances in Methods and Practices in Psychological ",
           "Science*, 1(4), 538-555. doi:10.1177/2515245918805755"), "")
}

.rp_report_block <- function(frame, RES, n_used, SHAPE = NULL) {
  if (is.null(RES)) return(character(0))
  hits <- RES[RES$selected & RES$level == "gene", , drop = FALSE]
  ## Surviving MODULES were written to hits.csv and then never mentioned. They
  ## are a declared level of the hierarchy, so they belong in the Results the
  ## same way genes do - and "nothing survived" must not be printed when a
  ## module did.
  modhits <- RES[RES$selected & RES$level == "module", , drop = FALSE]
  out <- c("", "## How to report this", "",
    paste("Copy and adapt. Every DMSA result needs three things a reader will",
          "not supply for themselves: the **carrying lens**, the **family** the",
          "unit was corrected inside, and the direction read as **expression",
          "tone**, not as methylation level. Reporting only an adjusted p makes",
          "a DMSA finding look like an ordinary EWAS hit, which is the one",
          "thing it is not."), "")
  if (!nrow(hits) && !nrow(modhits)) {
    out <- c(out,
      sprintf("**Results (nothing named).** No unit reached the naming bar (calibrated q <= %.2g)", frame$fdr %||% .10),
      "within its level-local family, so nothing is named.",
      "This is a calibrated null, not a failed run: report the families",
      "tested, the number of units in each, and the engine, so a reader can",
      "see what the design could have found.", "")
  } else {
    out <- c(out, "### Results", "")
    for (i in seq_len(nrow(hits))) {
      h <- hits[i, ]
      fam <- sum(RES$level == "gene" & RES$system == h$system &
                   RES$outcome == h$outcome & RES$n_probes > 0, na.rm = TRUE)
      out <- c(out, sprintf(
        paste0("> In the %s system, **%s** was named at FDR %.2g within its ",
               "%d-gene level-local family (dependence-calibrated ",
               "Benjamini-Hochberg q = %.3f on its best-of-three-lens ",
               "permutation p = %.4f; B = %s, %s engine), carried by the ",
               "**%s** lens. Its %d direction-called CpGs showed sign ",
               "concordance %.2f; the joint union test (max over the three ",
               "standardised lenses, family-corrected) gave p = %.4f%s, and ",
               "the ACAT omnibus across the three lenses was p = %.4f raw ",
               "and p = %.4f family-corrected by permutation minP on the ",
               "same stream%s. Under the pre-registered any-lens rule ",
               "(best lens's family-adjusted p = %.4f, %s within the lens) ",
               "it was %s. %s"),
        h$system, h$unit, frame$fdr %||% .10, fam, h$q_unit_dep, h$p_unit,
        format(frame$B, big.mark = ","), frame$weighting %||% "combined",
        h$best_lens, h$n_probes, h$concordance,
        h$p_unit_adj %||% NA_real_,
        if (h$exact_confirmed %in% TRUE)
          " - an exact, balanced family-wise claim on top of the FDR naming" else "",
        h$p_omnibus, h$p_omnibus_adj %||% NA_real_,
        if (h$omnibus_confirmed %in% TRUE)
          " - an exact family-wise claim from the three lenses COMBINED (ACAT rewards agreement between lenses but does not require it; one dominant lens can carry it)"
        else "",
        min(h$p_coherence_adj, h$p_composite_adj, h$p_diffuse_adj,
            na.rm = TRUE), frame$correction,
        if (h$prereg_any_lens %in% TRUE) "also named" else "not named",
        .rp_dir_sentence(frame, h$outcome, h$direction, h$unit)), "")
      out <- c(out, sprintf("  (%s carries it: %s.)", h$best_lens,
                            .rp_lens_gloss(h$best_lens)), "")
    }
  }
  if (nrow(modhits)) {
    if (!nrow(hits)) out <- c(out, "### Results", "")
    for (i in seq_len(nrow(modhits))) {
      h <- modhits[i, ]
      fam <- sum(RES$level == "module" & RES$system == h$system &
                   RES$outcome == h$outcome & RES$n_probes > 0, na.rm = TRUE)
      out <- c(out, sprintf(
        paste0("> At the module level, **%s** (%s system) was named at FDR ",
               "%.2g within its %d-module level-local family ",
               "(dependence-calibrated q = %.3f, unit p = %.4f; B = %s, %s ",
               "engine), carried by the **%s** lens. Its %d direction-called ",
               "CpGs showed sign concordance %.2f; the ACAT omnibus was ",
               "p = %.4f raw, %.4f family-corrected, and the joint union ",
               "test (family-corrected) = %.4f%s. %s A module pools the ",
               "probes of several genes, so this is not a claim about any ",
               "one of them."),
        h$unit, h$system, frame$fdr %||% .10, fam, h$q_unit_dep, h$p_unit,
        format(frame$B, big.mark = ","), frame$weighting %||% "combined",
        h$best_lens, h$n_probes, h$concordance,
        h$p_omnibus, h$p_omnibus_adj %||% NA_real_,
        h$p_unit_adj %||% NA_real_,
        if (h$exact_confirmed %in% TRUE)
          ", so the finding also stands at exact family-wise alpha" else "",
        .rp_dir_sentence(frame, h$outcome, h$direction, "this module")), "")
    }
  }
  out <- c(out, "### Methods", "",
    sprintf(paste0("> Set-level association was tested with Directional ",
      "Methylation Set Analysis (DMSA; dmsa R package v%s). Each probe effect ",
      "was aligned to its predicted consequence for gene expression through ",
      "m_j = d_j x w_g x (2p+ - 1) before aggregation, where d_j is the ",
      "probe's methylation-to-expression direction, w_g the gene's polarity ",
      "with respect to its system's activation tone, and p+ the confidence in ",
      "the direction call. Each unit was tested through three pre-specified ",
      "lenses - coherence, composite and diffuse - on one shared permutation ",
      "stream (B = %s), combined by an ACAT omnibus. Each unit's evidence was ",
      "its permutation p for the best of the three lenses, calibrated ",
      "against the permutation null of that maximum so the cross-lens ",
      "dependence is priced exactly. A unit was named when its ",
      "Benjamini-Hochberg q-value (Benjamini & Hochberg, 1995) within its ",
      "level-local family was at or below %.2g, with the BH level calibrated ",
      "on the same permutation stream so that a family without signal is ",
      "named in exactly that fraction of draws whatever the dependence ",
      "between units (resampling-calibrated FDR; Yekutieli & Benjamini, ",
      "1999). Named units that also passed the family-corrected joint ",
      "union test at alpha are marked exact-confirmed. The pre-registered ",
      "any-lens rule (any single lens's family-adjusted p below alpha, ",
      "Westfall-Young step-down %s within each lens) is reported beside the ",
      "naming rule, with its realized family-wise error measured from the ",
      "same stream. All corrections operate **within level-local ",
      "families only**: the named systems form one family, the genes of a ",
      "system another, and so on down the %s hierarchy, so a gene is never ",
      "corrected against genes in systems that were not selected. Analyses ",
      "used n = %d, the %s weighting engine, and the %s direction map. %s%s"),
      utils::packageVersion("dmsa"), format(frame$B, big.mark = ","),
      frame$fdr %||% .10, frame$correction,
      paste(names(frame$levels)[frame$levels], collapse = " > "),
      n_used, frame$weighting %||% "combined", frame$cpg_map,
      ## The moderated Methods block has always named the covariates and the
      ## block structure; the linear one never did. A reader of the main-effect
      ## report could not tell what was adjusted for - and in a run where a
      ## mistyped covariate was silently dropped, or where covariates were
      ## passed explicitly and all of them were unknown, the Methods text gave
      ## no hint that the model was unadjusted. State it, always.
      .rp_adjust_words(frame),
      if (!is.null(SHAPE) && nrow(SHAPE))
        " These unit-level tests are linear in the aligned score; departure from linearity was assessed separately and is reported below."
      else ""), "",
    paste("Direction is reported as expression tone, not methylation level: a",
          "negative direction means the observed methylation pattern is",
          "consistent with lower expression of the unit, whichever way the",
          "individual probes moved."), "")
  c(out, .rp_shape_report(frame, SHAPE, n_used), .rp_column_glossary(frame))
}

## spec 26/28: what every statistic column in tables/units.csv IS, which
## multiplicity it is corrected for, and - the point - whether it NAMES a
## finding. Re-ruled 2026-10-15: only the dependence-calibrated q names; the
## pre-registered any-lens rule is reported beside it; the union and omnibus
## tiers confirm; everything else is synthesis a reader may quote but must
## not mistake for a verdict.
.rp_column_glossary <- function(frame) {
  a <- frame$alpha %||% .05
  f <- frame$fdr %||% .10
  c("", "### The statistic columns, and which of them name a finding", "",
    paste0("`tables/units.csv` carries one row per unit and level. Read the ",
           "p-value columns with these roles in mind (FDR level = ",
           format(f), ", alpha = ", format(a), "):"),
    "",
    "| column | what it tests | corrected for | names a finding? |",
    "|---|---|---|---|",
    "| `p_coherence`, `p_composite`, `p_diffuse` | one lens each, raw permutation p | nothing (per-unit) | no - raw |",
    "| `p_unit` | the unit's best lens, calibrated against the permutation null of the best of three | the three lenses, exactly (measured dependence) | no - the evidence the naming rule ranks |",
    "| `q_unit` | Benjamini-Hochberg q-value of `p_unit` | the unit's level-local family (BH) | no - shown for reference |",
    "| `q_unit_dep` | `q_unit` with its level calibrated on the family's own permutation stream | lenses AND family, at the measured dependence | **YES - the naming statistic.** A unit is named when it is at or below the FDR level |",
    "| `selected` | `q_unit_dep` at or below the FDR level | - | `selected` IS the verdict |",
    "| `null_alarm_bh` | how often ordinary (uncalibrated) BH at the FDR level would name something in this family with no signal, measured from this run's permutations | - | no - shows what the calibration corrected |",
    "| `p_coherence_adj`, `p_composite_adj`, `p_diffuse_adj` | the same lens, family-adjusted (maxT/minP within that lens) | the unit's level-local family, within one lens | no - the pre-registered any-lens rule; `best_lens` says which lens carries the unit |",
    "| `selected_coherence`, `selected_composite`, `selected_diffuse`, `n_lenses_hit`, `any_lens_hit`, `prereg_any_lens` | the pre-registered any-lens rule spelled out per lens: named when ANY of the three is below alpha | - | no - reported beside the naming rule |",
    "| `fwer_realized` | the realized family-wise error of the pre-registered any-lens rule in this family, measured from this run's permutations | - | no - the disclosure that accompanies the pre-registered rule |",
    "| `p_unit_adj` | the joint max-over-lenses statistic, family-adjusted (the exact-alpha union test) | lenses AND family, exactly | no - CONFIRMS a named unit (dagger, `exact_confirmed`) |",
    "| `p_union_exact` | second-level Westfall-Young calibration of the any-lens rule's best lens | lenses AND family, exactly | no - optional synthesis; it does not gate naming |",
    "| `p_omnibus` | ACAT combination of the three raw lens p's | lenses only (one valid combined test); NOT the family | no - optional synthesis; never a verdict on its own |",
    "| `p_omnibus_adj` | permutation minP of the omnibus across the family | lenses AND family, exactly | no - CONFIRMS a named unit (double dagger, `omnibus_confirmed`) |",
    "",
    paste("So: name by `selected` (calibrated q at or below the FDR level",
          "within the unit's own family), report the carrying lens, the q and",
          "the unit p, add the dagger tiers where they hold, and report the",
          "pre-registered any-lens result beside it with `fwer_realized`. The",
          "omnibus and the other synthesis columns - quote them, do not",
          "name by them."), "")
}

## What the model was actually adjusted for, in one sentence. Used by the
## linear Methods block; the moderated block builds its own equivalent.
.rp_adjust_words <- function(frame) {
  cv <- setdiff(as.character(frame$covariates %||% character(0)), "")
  ## chip enters as a constructed factor, not as a named covariate
  if (.nz(frame$chip))
    cv <- unique(c(cv, "chip_f"))
  a <- if (!length(cv))
    paste0("**No covariates were entered, so the reported association is ",
           "unadjusted.**")
  else paste0("Covariates: ", paste(cv, collapse = ", "), ".")
  b <- if (length(frame$block_cols))
    paste0(" Cluster structure (", paste(frame$block_cols, collapse = ", "),
           ") was handled by restricting the permutation to exchangeable ",
           "blocks rather than by a random intercept - no mixed model was ",
           "fitted, so there is no variance component to report.")
  else " Observations were treated as independent; no block structure was used."
  d <- if (!is.null(frame$corrections) &&
           any(grepl("absent from data", frame$corrections$issue %||% "")))
    paste0(" NOTE: one or more requested covariates were not found in the data ",
           "and were dropped - see Design notes; the model is NOT adjusted for ",
           "them.")
  else ""
  paste0(a, b, d)
}

## ---- module evidence block for summary.md ---------------------------------
## Printed for every run that has module annotation, and expanded for the
## modules a run actually reported, so a reader sees the caveat attached to the
## finding rather than to the appendix.
.report_evidence_md <- function(frame, RES = NULL) {
  ev <- frame$module_evidence
  if (is.null(ev) || !nrow(ev) || all(is.na(ev$evidence_strength))) return(NULL)
  n_hi <- sum(ev$evidence_strength %in% "High", na.rm = TRUE)
  n_mo <- sum(ev$evidence_strength %in% "Moderate", na.rm = TRUE)
  flag <- grepl("heterogeneous|measurement_defined", ev$audit_status)
  out <- c("## Module evidence", "",
    sprintf("Module definitions come from %s: %d of %d module(s) in the selected system(s) rest on High evidence, %d on Moderate%s.",
            frame$selection$name, n_hi, nrow(ev), n_mo,
            if (any(flag, na.rm = TRUE))
              sprintf(", and %d carry a membership flag (heterogeneous or measurement-defined)",
                      sum(flag, na.rm = TRUE)) else ""),
    "")
  ## which modules did this run actually report on?
  rep_mods <- character(0)
  if (!is.null(RES) && "level" %in% names(RES))
    rep_mods <- unique(RES$unit[RES$level == "module"])
  show <- ev[ev$module_id %in% rep_mods | ev$module %in% rep_mods, , drop = FALSE]
  if (!nrow(show))
    show <- ev[ev$evidence_strength %in% "Moderate" | flag, , drop = FALSE]
  if (nrow(show)) {
    out <- c(out, "| module | evidence | status | references |", "|---|---|---|---|")
    for (i in seq_len(nrow(show)))
      out <- c(out, sprintf("| %s %s | %s | %s | %s |",
        show$module_id[i], show$module[i],
        ifelse(is.na(show$evidence_strength[i]), "-", show$evidence_strength[i]),
        ifelse(is.na(show$audit_status[i]), "-",
               gsub("_", " ", show$audit_status[i])),
        ifelse(is.na(show$citation_keys[i]), "-", show$citation_keys[i])))
    out <- c(out, "")
    note <- unique(stats::na.omit(show$evidence_note))
    if (length(note))
      out <- c(out, paste0("- ", utils::head(note, 4)), "")
  }
  u <- unique(stats::na.omit(ev$deep_search_url))
  if (length(u)) out <- c(out, paste("Evidence search:", u[1]), "")
  c(out, "The full table is `dmsa_evidence(frame)`.", "")
}
