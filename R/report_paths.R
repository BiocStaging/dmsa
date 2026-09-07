# ============================================================================
# dmsa_report() FOR THE FOUR STRUCTURAL METHODS
#
# The frame path (dmsa_frame -> dmsa_report) delivers figures/, tables/ and
# summary.md. The structural methods - dmsa_mediate(), dmsa_modmed(),
# dmsa_transmit(), dmsa_concordance() - printed to the console only, so a
# result could be read but not filed. dmsa_report() is now generic: handed
# any of those objects it writes the SAME three things, from the numbers the
# object already carries. No new statistics live in this file: every value
# written here was computed by the method that produced the object.
# ============================================================================

## a minimal stand-in for the frame's output options, so the writers of the
## frame report (.rp_write_table, .rp_dev, .rp_pal) can be reused verbatim
.mr_stub <- function(plot_type = "png", table_type = "docx",
                     palette = "viridis")
  list(plot_type = plot_type, table_type = table_type, palette = palette)

.mr_p <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "< .001",
                            formatC(p, format = "f", digits = 3)))
.mr_f <- function(x, d = 3) formatC(x, format = "f", digits = d)

## ---- pieces of a mediate / modmed object as data frames --------------------
## The engine estimates on the raw metric of the observed variables (a latent
## role already has variance 1). The registered runs and the methods paper
## print the standardized path, so both are given: b_std = b x SD(term) /
## SD(outcome), with the corrected SDs the moments carry.
.mr_paths_df <- function(x) {
  if (!length(x$paths)) return(NULL)
  d <- data.frame(path = names(x$paths),
                  b = vapply(x$paths, function(p) p$b %||% NA_real_, numeric(1)),
                  se = vapply(x$paths, function(p) p$se %||% NA_real_, numeric(1)),
                  p = vapply(x$paths, function(p) p$p %||% NA_real_, numeric(1)),
                  row.names = NULL, stringsAsFactors = FALSE)
  s <- x$scales
  if (!is.null(s) && length(s)) {
    inner <- sub("^[abc]'?\\[(.*)\\]$", "\\1", d$path)
    ## with several mediators an a-path is written a[X->M]; otherwise the
    ## a-path predicts the (single) mediator and every other path the outcome
    arrow <- grepl("->", inner, fixed = TRUE)
    term <- ifelse(arrow, sub("->.*$", "", inner), inner)
    out <- ifelse(arrow, sub("^.*->", "", inner),
                  ifelse(grepl("^a\\[", d$path), x$m[1], x$y %||% x$m[1]))
    num <- unname(s[term]); den <- unname(s[out])
    d$b_std <- ifelse(is.na(num) | is.na(den) | den == 0, NA_real_,
                      d$b * num / den)
  }
  d
}

.mr_indices_df <- function(x) {
  if (!length(x$indices)) return(NULL)
  g <- function(f) vapply(x$indices, function(i) i[[f]] %||% NA_real_,
                          numeric(1))
  lab <- vapply(x$indices, function(i)
                  as.character(i$kind %||% i$label %||% ""), character(1))
  data.frame(index = names(x$indices), estimate = g("est"), p_js = g("p_js"),
             p_holm = g("p_holm"), what = lab, row.names = NULL,
             stringsAsFactors = FALSE)
}

.mr_meas_df <- function(x) {
  if (!length(x$measurement)) return(NULL)
  do.call(rbind, lapply(names(x$measurement), function(s) {
    m <- x$measurement[[s]]
    data.frame(system = s, modules = m$p %||% NA_integer_,
               rung = m$model %||% NA_character_, omega = m$omega %||% NA_real_,
               error_variance = m$V %||% NA_real_,
               determinacy = m$determinacy %||% NA_real_,
               stepped = length(m$notes) > 0L, row.names = NULL,
               stringsAsFactors = FALSE)
  }))
}

.mr_cond_df <- function(x) {
  cd <- x$conditional
  if (is.null(cd)) return(NULL)
  if (is.data.frame(cd)) {
    names(cd)[names(cd) %in% c("w_sd", "w", "at")] <- "W_SD"
    return(cd)
  }
  if (is.list(cd) && !is.null(cd$indirect))
    return(data.frame(W_SD = cd$at %||% seq_along(cd$indirect),
                      indirect = as.numeric(cd$indirect), row.names = NULL))
  NULL
}

.mr_ci_df <- function(x) {
  ci <- x$ci
  if (is.null(ci) || is.null(ci$quantiles)) return(NULL)
  q <- as.data.frame(ci$quantiles)
  q$quantity <- rownames(ci$quantiles)
  q[, c(ncol(q), seq_len(ncol(q) - 1L))]
}

## ---- figures ---------------------------------------------------------------
## one forest: estimate with +/- 1.96 SE where an SE exists, filled when
## p < .05. The reading is the package's: existence from the permutation p,
## magnitude from the estimate.
.mr_fig_forest <- function(df, est, se, lab, p, file, stub, title, xlab) {
  if (is.null(df) || !nrow(df)) return(NULL)
  d <- data.frame(lab = rev(as.character(df[[lab]])), est = rev(df[[est]]),
                  se = if (is.null(se)) NA_real_ else rev(df[[se]]),
                  p = rev(df[[p]]))
  pal <- .rp_pal(stub, 4)
  f <- .rp_dev(stub, file, width = 7.5,
               height = max(2.6, 0.42 * nrow(d) + 1.5))
  on.exit(grDevices::dev.off(), add = TRUE)
  lo <- d$est - 1.96 * d$se; hi <- d$est + 1.96 * d$se
  lo[is.na(lo)] <- d$est[is.na(lo)]; hi[is.na(hi)] <- d$est[is.na(hi)]
  xr <- range(c(lo, hi, 0), na.rm = TRUE); pad <- diff(xr) * .08
  op <- graphics::par(mar = c(4.2, max(6, min(14, max(nchar(d$lab)) * .40)),
                              2.6, 1))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(NA, xlim = c(xr[1] - pad, xr[2] + pad),
                 ylim = c(.5, nrow(d) + .5), yaxt = "n", xlab = xlab,
                 ylab = "", main = title, bty = "n")
  graphics::abline(v = 0, col = "grey60")
  graphics::axis(2, at = seq_len(nrow(d)), labels = d$lab, las = 1,
                 cex.axis = .8, tick = FALSE)
  hit <- !is.na(d$p) & d$p < .05
  graphics::segments(lo, seq_len(nrow(d)), hi, seq_len(nrow(d)),
                     col = pal[2], lwd = 2)
  graphics::points(d$est, seq_len(nrow(d)), pch = ifelse(hit, 19, 21),
                   col = pal[2], bg = "white", cex = 1.3)
  graphics::legend("topright", c("p < .05", "p >= .05"), pch = c(19, 21),
                   col = pal[2], pt.bg = "white", bty = "n", cex = .8)
  f
}

.mr_fig_conditional <- function(cd, file, stub, wname) {
  if (is.null(cd) || !nrow(cd)) return(NULL)
  pal <- .rp_pal(stub, 4)
  f <- .rp_dev(stub, file, width = 6.5, height = 4.6)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- graphics::par(mar = c(4.4, 4.6, 2.6, 1)); on.exit(graphics::par(op),
                                                          add = TRUE)
  graphics::plot(cd$W_SD, cd$indirect, type = "b", pch = 19, lwd = 2,
                 col = pal[2], xlab = paste0(wname, " (SD from the mean)"),
                 ylab = "conditional indirect effect",
                 main = "The indirect effect at each level of the moderator",
                 bty = "n")
  graphics::abline(h = 0, col = "grey60")
  graphics::text(cd$W_SD, cd$indirect, labels = .mr_f(cd$indirect, 3),
                 pos = 3, cex = .75, col = pal[3])
  f
}

## transmission: mother and father slope side by side, per unit
.mr_fig_transmit <- function(x, file, stub, top = 20L) {
  d <- x$named
  if (is.null(d) || !nrow(d)) {
    d <- x$step0
    if (is.null(d)) return(NULL)
  }
  d <- utils::head(d[order(pmin(d$p_mother_maxT, d$p_father_maxT)), ], top)
  if (!nrow(d)) return(NULL)
  pal <- .rp_pal(stub, 5)
  f <- .rp_dev(stub, file, width = 7.5, height = max(3, .45 * nrow(d) + 1.6))
  on.exit(grDevices::dev.off(), add = TRUE)
  lab <- rev(paste0(d$unit, " (", d$level, ")"))
  bm <- rev(d$b_mother); bf <- rev(d$b_father)
  pm <- rev(d$p_mother_maxT); pf <- rev(d$p_father_maxT)
  xr <- range(c(bm, bf, 0), na.rm = TRUE); pad <- diff(xr) * .1
  op <- graphics::par(mar = c(4.2, max(6, min(16, max(nchar(lab)) * .40)),
                              2.6, 1))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(NA, xlim = c(xr[1] - pad, xr[2] + pad),
                 ylim = c(.5, length(lab) + .5), yaxt = "n",
                 xlab = "slope per latent SD of the child's unit", ylab = "",
                 main = "Transmission: mother and father, same unit", bty = "n")
  graphics::abline(v = 0, col = "grey60")
  graphics::axis(2, at = seq_along(lab), labels = lab, las = 1, cex.axis = .75,
                 tick = FALSE)
  graphics::points(bm, seq_along(lab) + .12, pch = ifelse(pm < .05, 19, 21),
                   col = pal[2], bg = "white", cex = 1.3)
  graphics::points(bf, seq_along(lab) - .12, pch = ifelse(pf < .05, 19, 21),
                   col = pal[4], bg = "white", cex = 1.3)
  graphics::legend("topright", c("mother", "father"), pch = 19,
                   col = c(pal[2], pal[4]), bty = "n", cex = .8)
  f
}

## concordance: the largest within-couple r at each level against its null
.mr_fig_concordance <- function(x, file, stub) {
  s <- x$summary
  if (is.null(s) || !nrow(s)) return(NULL)
  pal <- .rp_pal(stub, 5)
  f <- .rp_dev(stub, file, width = 7, height = 4.8)
  on.exit(grDevices::dev.off(), add = TRUE)
  lab <- paste(s$system, s$level, sep = " / ")
  op <- graphics::par(mar = c(7.5, 4.6, 2.6, 1)); on.exit(graphics::par(op),
                                                          add = TRUE)
  yr <- range(c(s$max_r, s$null_at_max, 0), na.rm = TRUE)
  graphics::plot(NA, xlim = c(.5, nrow(s) + .5),
                 ylim = c(min(yr) - .05, max(yr) + .08), xaxt = "n",
                 xlab = "", ylab = "within-couple r",
                 main = "Concordance against the batch-matched null",
                 bty = "n")
  graphics::abline(h = 0, col = "grey60")
  graphics::axis(1, at = seq_len(nrow(s)), labels = lab, las = 2, cex.axis = .75)
  graphics::segments(seq_len(nrow(s)), s$null_at_max, seq_len(nrow(s)),
                     s$max_r, col = "grey70", lwd = 2)
  graphics::points(seq_len(nrow(s)), s$null_at_max, pch = 4, col = pal[4],
                   cex = 1.1)
  graphics::points(seq_len(nrow(s)), s$max_r,
                   pch = ifelse(s$min_p_maxT < .05, 19, 21), col = pal[2],
                   bg = "white", cex = 1.4)
  graphics::legend("topleft", c("largest r (filled = named)", "null at that unit"),
                   pch = c(19, 4), col = c(pal[2], pal[4]), bty = "n", cex = .8)
  f
}

## ---- the shared writer -----------------------------------------------------
.mr_md_table <- function(df, digits = 3) {
  if (is.null(df) || !nrow(df)) return("(none)")
  d <- as.data.frame(df)
  for (j in seq_along(d)) {
    if (is.numeric(d[[j]]))
      d[[j]] <- if (grepl("^p($|_)", names(d)[j])) .mr_p(d[[j]])
                else formatC(d[[j]], format = "g", digits = digits)
    d[[j]] <- as.character(d[[j]])
  }
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

## spec 54, as on the frame path: never silently replace an earlier analysis.
## Checked BEFORE anything is drawn, so a refusal leaves the directory as it
## was rather than half-written.
.mr_check_outdir <- function(outdir, overwrite) {
  if (!isTRUE(overwrite)) {
    prev <- c(list.files(file.path(outdir, "figures"),
                         pattern = "\\.(png|pdf)$"),
              list.files(file.path(outdir, "tables"),
                         pattern = "\\.(csv|html|docx|rtf)$"),
              if (file.exists(file.path(outdir, "summary.md"))) "summary.md")
    if (length(prev))
      stop("`outdir` already contains a DMSA report: ",
           normalizePath(outdir, winslash = "/", mustWork = FALSE), "\n",
           length(prev), " existing output file(s). Refusing to replace an ",
           "earlier analysis. Point `outdir` at a new directory, or call ",
           "dmsa_report(x, outdir, overwrite = TRUE).\nNo report was written.",
           call. = FALSE)
  }
  dir.create(file.path(outdir, "figures"), recursive = TRUE,
             showWarnings = FALSE)
  dir.create(file.path(outdir, "tables"), recursive = TRUE,
             showWarnings = FALSE)
  invisible(outdir)
}

.mr_deliver <- function(x, outdir, tables, figures, head_lines, notes,
                        stub, overwrite) {
  written_t <- character()
  for (nm in names(tables)) {
    df <- tables[[nm]]
    if (is.null(df) || !nrow(as.data.frame(df))) next
    .rp_write_table(as.data.frame(df), file.path(outdir, "tables", nm), stub,
                    title = nm)
    written_t <- c(written_t, nm)
  }
  written_f <- unlist(Filter(Negate(is.null), figures), use.names = FALSE)
  md <- c(head_lines, "")
  for (nm in written_t) {
    md <- c(md, paste0("### ", nm), "",
            .mr_md_table(as.data.frame(tables[[nm]])), "")
  }
  if (length(written_f))
    md <- c(md, "### figures", "",
            paste0("- figures/", basename(written_f)), "")
  if (length(notes)) md <- c(md, "### how to read this", "", notes, "")
  md <- c(md, "", paste0("_dmsa ", utils::packageVersion("dmsa"), ", ",
                         format(Sys.time(), "%Y-%m-%d %H:%M"), "_"))
  writeLines(md, file.path(outdir, "summary.md"))
  out <- list(outdir = normalizePath(outdir, winslash = "/", mustWork = FALSE),
              tables = written_t, figures = basename(written_f),
              summary = "summary.md")
  message("dmsa report written to ", out$outdir, ": ",
          length(written_t), " table(s), ", length(out$figures),
          " figure(s), summary.md")
  invisible(structure(out, class = "dmsa_path_report"))
}

#' @export
print.dmsa_path_report <- function(x, ...) {
  cat("DMSA report -", x$outdir, "\n")
  cat("  tables : ", paste(x$tables, collapse = ", "), "\n")
  cat("  figures: ", paste(x$figures, collapse = ", "), "\n")
  cat("  summary: summary.md\n")
  invisible(x)
}

.mr_call_lines <- function(x, title) {
  cl <- if (is.null(x$call)) "(call not recorded)"
        else paste(deparse(x$call), collapse = "\n")
  c(paste0("# ", title), "", "### the call", "", "```r", cl, "```", "")
}

## ---- the generic -----------------------------------------------------------

#' Deliver a DMSA analysis: figures, tables and summary.md
#'
#' The second call of the two-call interface. Handed a \code{dmsa_frame} it
#' runs the frame and writes the report. Handed the result of one of the
#' structural methods - \code{\link{dmsa_mediate}}, \code{\link{dmsa_modmed}},
#' \code{\link{dmsa_transmit}}, \code{\link{dmsa_concordance}} - it writes the
#' same three things for that result: \code{tables/} (csv, html and, with
#' \pkg{gt}, docx), \code{figures/} and \code{summary.md}, which opens with
#' the call that produced them. Nothing is recomputed: every number written
#' is one the method already estimated.
#'
#' @param frame A \code{dmsa_frame}, or the result of a structural method.
#' @param outdir Where to write (methods only; the frame carries its own
#'   \code{outdir}). Refuses to overwrite an existing report unless
#'   \code{overwrite = TRUE}.
#' @param overwrite Replace an existing report in \code{outdir}.
#' @param plot_type,table_type,palette Output options for the method reports
#'   (the frame supplies its own).
#' @param ... Passed to the method.
#' @return For a frame, the report object it has always returned; for a
#'   structural result, an invisible list of what was written.
#' @export
dmsa_report <- function(frame, ...) UseMethod("dmsa_report")

#' @export
dmsa_report.default <- function(frame, ...)
  stop("dmsa_report() takes a dmsa_frame or the result of dmsa_mediate(), ",
       "dmsa_modmed(), dmsa_transmit() or dmsa_concordance(); got ",
       paste(class(frame), collapse = "/"), call. = FALSE)

#' @rdname dmsa_report
#' @export
dmsa_report.dmsa_modmed <- function(frame, outdir = ".", overwrite = FALSE,
                                    plot_type = c("png", "pdf"),
                                    table_type = c("docx", "html", "rtf"),
                                    palette = "viridis", ...) {
  x <- frame
  stub <- .mr_stub(match.arg(plot_type), match.arg(table_type), palette)
  .mr_check_outdir(outdir, overwrite)
  paths <- .mr_paths_df(x); idx <- .mr_indices_df(x); cond <- .mr_cond_df(x)
  tabs <- list(paths = paths, indices = idx, conditional_indirect = cond,
               measurement = .mr_meas_df(x), bootstrap = .mr_ci_df(x),
               coverage = x$coverage)
  dir.create(file.path(outdir, "figures"), recursive = TRUE,
             showWarnings = FALSE)
  figs <- list(
    .mr_fig_forest(paths, "b", "se", "path", "p",
                   file.path(outdir, "figures", "paths"), stub,
                   "Structural paths (Freedman-Lane p)",
                   "b (per latent SD where the role is latent)"),
    if (!is.null(idx))
      .mr_fig_forest(idx, "estimate", NULL, "index", "p_js",
                     file.path(outdir, "figures", "indices"), stub,
                     "Indices (joint significance)", "estimate"),
    if (!is.null(cond) && !is.null(x$w))
      .mr_fig_conditional(cond, file.path(outdir, "figures", "conditional"),
                          stub, x$w))
  ttl <- if (is.null(x$w)) "Mediation through a DMSA tone"
         else "Moderated mediation through a DMSA tone"
  head_lines <- c(.mr_call_lines(x, ttl), "### design", "",
    sprintf("- n = %d, B = %d permutations, blocks: %s, pivot: %s", x$n, x$B,
            x$perm_blocks, x$pivot %||% "homoscedastic"),
    sprintf("- roles: %s", paste(sprintf("%s [%s]", names(x$roles),
                                         unlist(x$roles)), collapse = ", ")),
    sprintf("- map: %s", x$map_source %||% "not recorded"),
    sprintf("- measurement rule: %s", x$rule %||% "not recorded"),
    if (!is.null(x$moments) && isTRUE(x$moments$lambda_triggered))
      sprintf("- NOTE: offsets shrunk to %.2f of the full correction to keep the moment matrix positive definite",
              x$moments$shrink))
  notes <- c(
    "- Existence is read from the permutation p, magnitude from the corrected estimate with its interval.",
    "- An index is tested by joint significance: every factor must reject.",
    "- Latent roles carry the measurement offsets; a path between two tones is on the true-score scale.",
    if (!is.null(x$w))
      "- An interaction claim should also survive the quadratic safeguard (`quadratic = TRUE`).")
  .mr_deliver(x, outdir, tabs, figs, head_lines, notes, stub, overwrite)
}

#' @rdname dmsa_report
#' @export
dmsa_report.dmsa_mediate <- dmsa_report.dmsa_modmed

#' @rdname dmsa_report
#' @export
dmsa_report.dmsa_transmit <- function(frame, outdir = ".", overwrite = FALSE,
                                      plot_type = c("png", "pdf"),
                                      table_type = c("docx", "html", "rtf"),
                                      palette = "viridis", ...) {
  x <- frame; d <- x$design
  stub <- .mr_stub(match.arg(plot_type), match.arg(table_type), palette)
  .mr_check_outdir(outdir, overwrite)
  tabs <- list(per_level = x$summary, named_units = x$named,
               measurement = x$measurement, all_units = x$step0)
  dir.create(file.path(outdir, "figures"), recursive = TRUE,
             showWarnings = FALSE)
  figs <- list(.mr_fig_transmit(x, file.path(outdir, "figures", "transmission"),
                                stub))
  head_lines <- c(.mr_call_lines(x, "Parent-to-child transmission at every level"),
    "### design", "",
    sprintf("- %d triads (%d families), B = %d, blocks: family-exchange, pivot: hc",
            d$n, d$families, d$B),
    sprintf("- maxT within the level-local family; p_tree across the levels of the system"),
    sprintf("- map: %s", d$map_source),
    sprintf("- rule (child): %s", d$rule_child),
    sprintf("- rule (parents): %s", d$rule_parents))
  notes <- c(
    "- A unit is named when its maxT p within its level-local family is below .05; p_tree corrects across every level of the system.",
    "- A father slope of the mother's size or larger argues against an intrauterine route: genotype, household, or the paternal germline.",
    "- The probe level is observed; gene, module and system carry the measurement offsets.")
  .mr_deliver(x, outdir, tabs, figs, head_lines, notes, stub, overwrite)
}

#' @rdname dmsa_report
#' @export
dmsa_report.dmsa_concordance <- function(frame, outdir = ".",
                                         overwrite = FALSE,
                                         plot_type = c("png", "pdf"),
                                         table_type = c("docx", "html", "rtf"),
                                         palette = "viridis", ...) {
  x <- frame; d <- x$design
  stub <- .mr_stub(match.arg(plot_type), match.arg(table_type), palette)
  .mr_check_outdir(outdir, overwrite)
  tabs <- list(per_level = x$summary, named_units = x$named,
               measurement = x$measurement, step0 = x$step0,
               surface = x$surface, exposure = x$exposure)
  dir.create(file.path(outdir, "figures"), recursive = TRUE,
             showWarnings = FALSE)
  figs <- list(.mr_fig_concordance(x, file.path(outdir, "figures",
                                                "concordance"), stub))
  head_lines <- c(.mr_call_lines(x, "Dyadic concordance at every level"),
    "### design", "",
    sprintf("- %d couples, NP = %d re-pairings within %s (%d batches)",
            d$couples, d$NP, d$batch, d$batches),
    sprintf("- maxT within the level-local family; p_tree across the levels of the system"),
    sprintf("- map: %s", d$map_source), sprintf("- rule: %s", d$rule))
  notes <- c(
    "- The null re-pairs partners WITHIN the processing batch: partners who share a batch make a batch effect look like concordance.",
    "- The re-pairing is the permutation stream, so the level-local maxT and p_tree come from the same draws.",
    if (!is.null(x$surface))
      "- Surface: existence from the uncorrected a4, magnitude from the corrected one with its interval (a second-order product of two error-prone scores is hard to identify at these n).",
    if (!is.null(x$exposure))
      "- Exposure: gamma3 is the change in the within-couple slope with the exposure; the exposure is never permuted.")
  .mr_deliver(x, outdir, tabs, figs, head_lines, notes, stub, overwrite)
}
