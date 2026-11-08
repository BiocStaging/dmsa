# dmsa gene-to-system polarity: w_g
# --------------------------------------------------------------------------
# w_g answers one question: does higher expression of this gene RAISE or LOWER
# its system's activation tone? It enters the alignment chain as a multiplier
#
#     m_j = d_j * w_g * (2 p_plus - 1)
#
# so a wrong w_g does not weaken a finding, it INVERTS it. Everything in this
# file exists to make that risk visible rather than convenient: every sign
# carries its source, its confidence and a review flag, and the print methods
# report the evidence grade of the units actually being tested instead of
# leaving it in a data file.
#
# The rule that generates the signs, stated once because it is the one thing
# most easily got wrong: a gene is signed against the SYSTEM'S DECLARED TONE, not
# against its immediate molecular partner. FKBP5 restrains GR; GR mediates the
# HPA negative feedback; so under axis-drive tone FKBP5 is +1, a brake on a
# brake. Reading FKBP5 as "-1 because it inhibits GR" scores it against the wrong
# referent and inverts the module. Hence the role vocabulary carries
# brake-of-brake and feedback-enabler as first-class labels.
#
# Two resources were tested against this table and rejected, which is worth
# recording so nobody spends the effort again. The Pathway Commons GMT
# (direction words in pathway names) agreed with the existing signs on only 65%
# of signed rows - barely above chance for a two-class call - so it is not used
# as a source. And chaining a signed transcription-factor edge to a signed target
# inside the same system inverts on negative-feedback loops: TRRUST makes ARNTL,
# CLOCK and AR come out negative, and all three are ANCHORS of their own systems.
# A chain that contradicts an anchor has run through the loop the system is built
# on; those verdicts are discarded rather than flagged, because the failure is in
# the method and not in the sign.
#
# That rule earned its keep. All 38 database-vs-derived contradictions were then
# resolved against the primary literature, and the databases lost 34 of them.
# Four were traced to defects in the database record itself: TRRUST codes a
# dominant-negative construct's blocking as IKBKB repressing VCAM1 (sign inverted
# at source), carries no RORA-ARNTL edge at all, sources its MBD1 rows from two
# retracted papers, and attributes the PAX8-PPARG fusion oncoprotein's effect to
# wild-type PAX8. Had any been applied automatically it would have inverted a
# module.
#
# Design rule, stated once: a database may raise the EVIDENCE GRADE of a sign
# that already agrees with it, but it may not flip a sign on its own. Signed
# interaction resources default to "activating" when they do not know (CollecTRI
# and DoRothEA both document this), and pairwise gene-gene signs are not
# gene-to-system signs. Contradictions are therefore surfaced for adjudication,
# never applied silently. See dmsa_polarity_review().

.POL_REQ  <- c("system_id", "gene", "w_g")
.POL_OPT  <- c("system_short", "system", "module_id", "module", "role",
               "confidence", "w_g_source", "anchor", "evidence", "citation",
               "review_flag")
## sources ordered by how much weight a reader should give them
.POL_GRADE <- c(curated = "curated", gtopdb = "database", signor = "database",
                pathwaycommons = "heuristic",
                omnipath = "database", trrust = "database",
                uniprot = "database", go = "database",
                literature = "literature", rule = "heuristic",
                unresolved = "none")

## ---- role vocabulary and the direction each role commits to ---------------
## The point of this table is the last group. An accessory subunit, a downstream
## transducer, a capacity marker and a readout all take their sign from whatever
## they serve or report, so their role word commits to no direction of its own -
## GNAI2 is an "accessory" at -1 because the D2 receptor it serves is -1, and
## TERT is "capacity" at -1 because more telomerase means LESS ageing tone.
## Only the first two groups let a role contradict a sign, and when one does it
## means the gene was scored against its own sub-process instead of the system.
.POL_ROLE_POS <- c("driver", "synthesis", "brake-of-brake", "driver-adjacent",
                   "cotransmitter", "specificity-control", "feedback-enabler-of-brake")
.POL_ROLE_NEG <- c("brake", "clearance", "catabolism", "antagonist", "autoreceptor",
                   "feedback-enabler", "buffer", "brake-of-driver")
.POL_ROLE_ZERO <- c("off-axis", "unresolved")
## sign inherited from what the gene serves, reports, or transduces
.POL_ROLE_INHERIT <- c("accessory", "transducer", "capacity", "readout",
                       "feedback", "ambiguous", "axis-dependent",
                       "opposing-transducer", "transducer-with-caveat")

## A gene must be signed toward the SYSTEM'S tone, never toward its module's own
## sub-process. Two independent signatures betray the latter, and both are checked.
.pol_toward_system <- function(p) {
  role <- tolower(trimws(ifelse(is.na(p$role), "", p$role)))
  ## (1) a direction-committing role that disagrees with its own sign
  bad <- (role %in% .POL_ROLE_POS  & p$w_g < 0) |
         (role %in% .POL_ROLE_NEG  & p$w_g > 0) |
         (role %in% .POL_ROLE_ZERO & p$w_g != 0)
  unknown <- setdiff(unique(role[nzchar(role)]),
                     c(.POL_ROLE_POS, .POL_ROLE_NEG, .POL_ROLE_ZERO,
                       .POL_ROLE_INHERIT))
  ## (2) a module whose label names a braking process but whose every signed gene
  ## is +1, or the reverse. Word boundaries matter: "sex deTERMINATION" must not
  ## match "termination".
  flip <- character(0)
  if (all(c("module_id", "module") %in% names(p))) {
    br <- "\\b(clearance|degradation|catabolism|inactivation|reuptake|feedback|negative|antagonism|buffering|disposal|breakdown|efflux)\\b"
    dr <- "\\b(synthesis|biosynthesis|drive|release|activation|production|secretion|reception)\\b"
    for (m in unique(p$module_id)) {
      d <- p[p$module_id == m, , drop = FALSE]
      lbl <- d$module[1]; sg <- d$w_g[d$w_g != 0]
      if (length(sg) < 2 || is.na(lbl)) next
      uni <- all(sg > 0) || all(sg < 0)
      if (!uni) next
      isbr <- grepl(br, lbl, ignore.case = TRUE)
      isdr <- grepl(dr, lbl, ignore.case = TRUE)
      if (isbr && !isdr && sg[1] > 0) flip <- c(flip, m)
      if (isdr && !isbr && sg[1] < 0) flip <- c(flip, m)
    }
  }
  list(role_conflict = which(bad), unknown_roles = unknown,
       module_flip = unique(flip))
}

.pol_builtin <- function() {
  f <- system.file("extdata", "alpha_polarity_2026c.csv", package = "dmsa")
  if (!nzchar(f) || !file.exists(f)) NULL else f
}

#' Load a gene-to-system polarity table
#'
#' @param x \code{"alpha"} for the bundled Project Alpha 2026c polarity table, a
#'   path to a CSV in the same shape, a \code{data.frame}, or a
#'   \code{dmsa_sets} / \code{dmsa_selection} / \code{dmsa_frame} object whose
#'   attached polarity should be used - so a narrowed selection reviews only its
#'   own genes. \code{NULL}
#'   returns \code{NULL}, which is how a caller opts out of polarity entirely
#'   (system-level scores then weight every gene +1 and say so).
#' @param sets Optional \code{dmsa_sets} or \code{dmsa_selection}. When given,
#'   coverage is reported against it: which of its genes have a sign, which do
#'   not, and at what evidence grade.
#' @return An object of class \code{dmsa_polarity}, or \code{NULL}.
#' @seealso \code{dmsa_polarity_review()} for the rows needing adjudication,
#'   \code{dmsa_polarity_fetch()} to draft signs for a panel of your own,
#'   \code{dmsa_polarity_check()} to validate a candidate table.
#' @details Signs are read against the system's declared activation tone, not
#'   against a gene's immediate partner. A gene that inhibits a brake therefore
#'   carries \code{w_g = +1}: FKBP5 restrains GR, GR mediates HPA negative
#'   feedback, so more FKBP5 means more axis drive. The \code{role} vocabulary
#'   names this explicitly (\code{brake-of-brake},
#'   \code{feedback-enabler}) because it is the error that inverts a module.
#' @examples
#' pol <- dmsa_polarity()
#' pol
#' dmsa_polarity_review(pol)
#' @export
dmsa_polarity <- function(x = "alpha", sets = NULL) {
  if (is.null(x)) return(NULL)
  ## A user holding a sets or selection object should be able to hand it straight
  ## over - the print banner tells them to call dmsa_polarity_review(), and the
  ## obvious argument is the object in their hand.
  if (inherits(x, "dmsa_frame") && !is.null(x$polarity_table))
    return(x$polarity_table)
  if (inherits(x, c("dmsa_sets", "dmsa_selection"))) {
    if (is.null(x$polarity))
      stop("no polarity table is attached to this ",
           class(x)[1], "; pass one via dmsa_sets(polarity = ...)", call. = FALSE)
    return(x$polarity)
  }
  if (inherits(x, "dmsa_polarity")) return(x)
  src <- NULL; name <- NULL
  if (is.character(x) && length(x) == 1L && identical(x, "alpha")) {
    f <- .pol_builtin()
    if (is.null(f)) {
      warning("the bundled Alpha polarity table is not installed; ",
              "system-level scores would weight every gene +1", call. = FALSE)
      return(NULL)
    }
    p <- utils::read.csv(f, stringsAsFactors = FALSE,
                         colClasses = c(system_id = "character",
                                        module_id = "character"))
    src <- f; name <- "Project Alpha 2026c (audited)"
  } else if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("polarity file not found: ", x, call. = FALSE)
    p <- utils::read.csv(x, stringsAsFactors = FALSE,
                         colClasses = c(system_id = "character",
                                        module_id = "character"))
    src <- x; name <- basename(x)
  } else {
    p <- as.data.frame(x, stringsAsFactors = FALSE); name <- "user data.frame"
  }

  miss <- setdiff(.POL_REQ, names(p))
  if (length(miss))
    stop("polarity table is missing required column(s): ",
         paste(miss, collapse = ", "),
         "\nrequired: system_id, gene, w_g", call. = FALSE)
  p$system_id <- as.character(p$system_id)
  p$gene <- as.character(p$gene)
  p$w_g <- suppressWarnings(as.numeric(p$w_g))
  if (any(!is.finite(p$w_g)))
    stop(sum(!is.finite(p$w_g)), " row(s) have a non-numeric w_g", call. = FALSE)
  if (any(abs(p$w_g) > 1))
    stop("w_g must lie in [-1, 1]; ", sum(abs(p$w_g) > 1), " row(s) do not",
         call. = FALSE)
  for (v in .POL_OPT) if (!v %in% names(p)) p[[v]] <- NA
  p$anchor <- .pol_logical(p$anchor)
  p$w_g_source[is.na(p$w_g_source)] <- "unstated"
  p$grade <- unname(.POL_GRADE[p$w_g_source])
  p$grade[is.na(p$grade)] <- "unstated"

  dup <- duplicated(p[c("system_id", "gene")])
  if (any(dup))
    stop(sum(dup), " duplicate system_id + gene row(s); a gene may appear in ",
         "several systems but only once within one", call. = FALSE)

  cov <- NULL
  if (!is.null(sets)) cov <- .pol_coverage(p, sets)

  structure(list(polarity = p, coverage = cov, source = src, name = name),
            class = "dmsa_polarity")
}

.pol_logical <- function(x) {
  if (is.logical(x)) return(x)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
}

## Coverage against a cascade, per system. Genes with no sign are the honest
## denominator: they are the ones a system-level score has to weight +1 by
## default, which is an assumption and is reported as one.
.pol_coverage <- function(p, sets) {
  cas <- if (is.list(sets) && !is.null(sets$cascade)) sets$cascade
         else as.data.frame(sets, stringsAsFactors = FALSE)
  keep <- unique(cas[, intersect(c("system_id", "system_short", "gene"),
                                 names(cas))])
  keep$system_id <- as.character(keep$system_id)
  m <- match(paste(keep$system_id, keep$gene),
             paste(p$system_id, p$gene))
  keep$w_g <- p$w_g[m]
  keep$grade <- p$grade[m]
  keep$grade[is.na(keep$grade)] <- "missing"
  agg <- lapply(split(keep, keep$system_id), function(d) data.frame(
    system_id = d$system_id[1],
    system_short = if ("system_short" %in% names(d)) d$system_short[1] else NA,
    n_genes = nrow(d),
    n_signed = sum(!is.na(d$w_g) & d$w_g != 0),
    n_zero = sum(!is.na(d$w_g) & d$w_g == 0),
    n_missing = sum(is.na(d$w_g)),
    n_curated = sum(d$grade == "curated"),
    n_database = sum(d$grade == "database"),
    n_literature = sum(d$grade == "literature"),
    n_heuristic = sum(d$grade == "heuristic"),
    stringsAsFactors = FALSE))
  out <- do.call(rbind, agg)
  out <- out[order(suppressWarnings(as.numeric(out$system_id))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @export
print.dmsa_polarity <- function(x, ...) {
  p <- x$polarity
  cat("dmsa gene-to-system polarity:", x$name, "\n")
  cat(sprintf("  %d gene-system rows | %d systems | +1: %d, -1: %d, 0: %d | %d anchor(s)\n",
              nrow(p), length(unique(p$system_id)), sum(p$w_g > 0),
              sum(p$w_g < 0), sum(p$w_g == 0), sum(p$anchor, na.rm = TRUE)))
  g <- table(factor(p$grade, levels = c("curated", "database", "literature",
                                        "heuristic", "none", "unstated")))
  g <- g[g > 0]
  cat("  evidence grade:",
      paste(sprintf("%s %d (%.0f%%)", names(g), g, 100 * g / nrow(p)),
            collapse = " | "), "\n")
  nf <- sum(nzchar(stats::na.omit(p$review_flag)))
  if (nf)
    cat(sprintf("  %d row(s) flagged for review - dmsa_polarity_review()\n", nf))
  if (!is.null(x$coverage)) {
    cv <- x$coverage
    cat(sprintf("  coverage against the selected cascade: %d of %d gene-system pairs signed, %d unsigned\n",
                sum(cv$n_signed), sum(cv$n_genes), sum(cv$n_missing)))
    if (sum(cv$n_missing))
      cat("   unsigned genes are weighted +1 in a system score - that is an assumption, not a measurement\n")
  }
  h <- sum(p$grade == "heuristic")
  if (h / nrow(p) > .25)
    cat(sprintf("  NOTE %.0f%% of signs rest on a functional-class heuristic rather than a database or paper.\n",
                100 * h / nrow(p)),
        "       Treat system-level direction claims in those systems as provisional.\n",
        sep = "")
  invisible(x)
}

#' Rows of a polarity table that need a human decision
#'
#' Three kinds of row should never be used without someone looking at them: a
#' sign the curator and a database disagree about, a sign no evidence could
#' settle, and a non-zero sign held at low confidence. This returns them.
#'
#' @param x A \code{dmsa_polarity}, or anything \code{dmsa_polarity()} accepts -
#'   including a \code{dmsa_selection}, which scopes the review to that
#'   selection's genes.
#' @param which \code{"all"}, \code{"disagreement"} (a database contradiction or
#'   a framing choice only the analyst can settle), \code{"unresolved"} or
#'   \code{"low_confidence"}.
#' @return data.frame of class \code{dmsa_polarity_review}.
#' @examples
#' ## The three kinds of row nobody should use unseen: a database contradiction,
#' ## a sign no evidence settles, and a non-zero sign held at low confidence.
#' flagged <- dmsa_polarity_review(dmsa_select(systems = "hpa"))
#' table(flagged$review_flag)
#' ## a disagreement is a framing call the analyst has to make, not a bug
#' dmsa_polarity_review(which = "disagreement")$gene
#' @export
dmsa_polarity_review <- function(x = "alpha",
                                 which = c("all", "disagreement",
                                           "unresolved", "low_confidence")) {
  which <- match.arg(which)
  pol <- dmsa_polarity(x)
  if (is.null(pol)) return(invisible(NULL))
  p <- pol$polarity
  fl <- ifelse(is.na(p$review_flag), "", p$review_flag)
  sel <- switch(which,
    all = nzchar(fl),
    disagreement = grepl("adjudicate|framing", fl),
    unresolved = grepl("unresolved", fl),
    low_confidence = grepl("low_confidence", fl))
  out <- p[sel, c("system_id", "system_short", "module_id", "gene", "w_g",
                  "role", "confidence", "w_g_source", "review_flag",
                  "evidence", "citation")]
  out <- out[order(match(sub("_.*", "", out$review_flag), c("pi", "unresolved", "low")),
                   suppressWarnings(as.numeric(out$system_id)), out$gene), ,
             drop = FALSE]
  rownames(out) <- NULL
  structure(out, class = c("dmsa_polarity_review", "data.frame"))
}

#' @export
print.dmsa_polarity_review <- function(x, n = 20, ...) {
  if (!nrow(x)) { cat("no polarity rows flagged for review\n"); return(invisible(x)) }
  tb <- table(x$review_flag)
  cat(sprintf("polarity rows needing a decision: %d\n", nrow(x)))
  for (i in seq_along(tb))
    cat(sprintf("  %-38s %d\n", names(tb)[i], tb[i]))
  adj <- x[grepl("adjudicate|framing", x$review_flag), , drop = FALSE]
  if (nrow(adj)) {
    cat("\nthe ones that change a sign if you rule the other way:\n")
    for (i in seq_len(nrow(adj)))
      cat(sprintf(" %-9s sys %-3s w_g=%+d  %s\n     %s\n", adj$gene[i],
                  adj$system_id[i], as.integer(adj$w_g[i]),
                  adj$w_g_source[i], substr(adj$evidence[i], 1, 150)))
  }
  rest <- x[!grepl("adjudicate|framing", x$review_flag), , drop = FALSE]
  if (nrow(rest)) {
    k <- min(nrow(rest), n)
    cat(sprintf("\nfirst %d of %d lower-priority rows:\n", k, nrow(rest)))
    for (i in seq_len(k))
      cat(sprintf(" %-9s sys %-3s w_g=%+d %-8s %-10s %s\n", rest$gene[i],
                  rest$system_id[i], as.integer(rest$w_g[i]),
                  rest$confidence[i], rest$w_g_source[i],
                  substr(rest$evidence[i], 1, 60)))
  }
  invisible(x)
}

## ---- validation ----------------------------------------------------------

#' Validate a candidate polarity table
#'
#' @param x Path, \code{data.frame}, or \code{dmsa_polarity}.
#' @param sets Optional cascade to check coverage against.
#' @param verbose Print the report. Default \code{TRUE}.
#' @return Invisibly, a list with \code{ok} and one entry per check.
#' @examples
#' ## An anchor gene is what defines "more activation" for its system, so an
#' ## anchor that brakes its own system is a contradiction, not a small error.
#' cand <- data.frame(system_id = "HPA", gene = c("CRH", "POMC", "NR3C1"),
#'                    w_g = c(-1, 1, -1), anchor = c(TRUE, FALSE, FALSE),
#'                    w_g_source = "curated")
#' chk <- dmsa_polarity_check(cand)
#' chk$ok
#' dmsa_polarity_check("alpha", verbose = FALSE)$ok
#' @export
dmsa_polarity_check <- function(x = "alpha", sets = NULL, verbose = TRUE) {
  pol <- try(dmsa_polarity(x), silent = TRUE)
  if (inherits(pol, "try-error")) {
    if (verbose) cat("  [FAIL] table could not be loaded:",
                     conditionMessage(attr(pol, "condition")), "\n")
    return(invisible(list(ok = FALSE)))
  }
  p <- pol$polarity
  res <- list(); fail <- character()
  add <- function(nm, ok, detail = "") {
    res[[nm]] <<- list(ok = ok, detail = detail)
    if (!ok) fail <<- c(fail, nm)
  }
  add("w_g in [-1, 1]", all(abs(p$w_g) <= 1))
  add("no duplicate gene within a system",
      !anyDuplicated(p[c("system_id", "gene")]))
  ## An anchor defines what "more activation" means. An anchor that lowers its
  ## own system's tone is not a subtle error, it is a contradiction in terms.
  bad <- p$anchor & p$w_g <= 0
  add("every anchor has w_g > 0", !any(bad, na.rm = TRUE),
      if (any(bad, na.rm = TRUE))
        paste("offenders:", paste(p$gene[which(bad)], collapse = ", ")) else "")
  ## A system with no anchor has no defined tone, so its system-level sign is
  ## uninterpretable even when every gene is signed.
  noanch <- setdiff(unique(p$system_id), unique(p$system_id[p$anchor]))
  add("every system has at least one anchor", !length(noanch),
      if (length(noanch)) paste("systems without an anchor:",
                                paste(noanch, collapse = ", ")) else "")
  ## A system whose genes all point the same way cannot produce a two-sided
  ## aggregate; that is worth knowing before it is interpreted as a finding.
  onesided <- names(which(vapply(split(p$w_g, p$system_id),
    function(v) all(v >= 0) || all(v <= 0), TRUE)))
  add("no system is entirely one-sided", !length(onesided),
      if (length(onesided)) paste("one-sided systems:",
                                  paste(onesided, collapse = ", ")) else "")
  ## Is every sign pointing at the system rather than at its module's own
  ## sub-process? This is the failure that inverts a module, so it is checked
  ## two independent ways rather than assumed.
  ts <- .pol_toward_system(p)
  add("no role contradicts its own sign", !length(ts$role_conflict),
      if (length(ts$role_conflict))
        paste0(length(ts$role_conflict), " row(s): ",
               paste(utils::head(paste0(p$gene[ts$role_conflict], " (",
                                        p$role[ts$role_conflict], ", w_g ",
                                        sprintf("%+d", as.integer(p$w_g[ts$role_conflict])),
                                        ")"), 4), collapse = "; ")) else "")
  add("role vocabulary is closed", !length(ts$unknown_roles),
      if (length(ts$unknown_roles))
        paste("roles outside the vocabulary:",
              paste(ts$unknown_roles, collapse = ", ")) else "")
  add("no module signed toward its own sub-process", !length(ts$module_flip),
      if (length(ts$module_flip))
        paste("module(s) whose label points against every sign inside it:",
              paste(ts$module_flip, collapse = ", ")) else "")
  add("source stated for every row", !any(p$w_g_source == "unstated"),
      paste(sum(p$w_g_source == "unstated"), "row(s) without w_g_source"))
  hp <- mean(p$grade == "heuristic")
  add("most signs rest on evidence, not heuristics", hp <= 0.5,
      sprintf("%.0f%% heuristic", 100 * hp))

  cov <- if (!is.null(sets)) .pol_coverage(p, sets) else NULL
  if (!is.null(cov))
    add("every cascade gene has a sign", sum(cov$n_missing) == 0,
        paste(sum(cov$n_missing), "gene-system pair(s) unsigned"))

  if (verbose) {
    cat(sprintf("polarity: %d rows | %d systems | %.0f%% heuristic | %d flagged\n",
                nrow(p), length(unique(p$system_id)), 100 * hp,
                sum(nzchar(stats::na.omit(p$review_flag)))))
    for (nm in names(res))
      cat(sprintf("  [%s] %-44s %s\n", if (res[[nm]]$ok) "ok" else "WARN", nm,
                  res[[nm]]$detail))
    cat(if (length(fail))
      paste0("  -> ", length(fail), " check(s) raised a warning; none is fatal, ",
             "but each changes how a system-level sign should be read\n")
      else "  -> clean\n")
  }
  invisible(list(ok = !length(fail), checks = res, coverage = cov))
}

## ---- drafting polarity for a user's own panel -----------------------------

.POL_SOURCES <- list(
  go = list(
    name = "Gene Ontology",
    url = function(sym) paste0(
      "https://api.geneontology.org/api/bioentity/gene/HGNC%3A", sym,
      "/function"),
    note = "positive/negative regulation of <process> terms give a gene-to-process sign directly"),
  omnipath = list(
    name = "OmniPath",
    url = function(sym) paste0(
      "https://omnipathdb.org/interactions?genesymbols=yes&partners=", sym,
      "&organisms=9606&fields=sources,references,curation_effort&directed=1"),
    note = "signed pairwise edges; consensus_stimulation at curation_effort 1 is often a default, not a finding"),
  trrust = list(
    name = "TRRUST v2",
    url = function(sym) "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
    note = "explicit Activation / Repression / Unknown per TF-target pair; the Unknown is the reason to prefer it"),
  signor = list(
    name = "SIGNOR 3.0",
    url = function(sym) "https://signor.uniroma2.it/API/getHumanData.php",
    note = "rows with TYPEB = phenotype and EFFECT = up-regulates are gene-to-phenotype signs, no path chaining needed")
)

#' Draft polarity signs for a panel from public databases
#'
#' Queries public signed-relationship resources for the genes of a cascade and
#' returns a DRAFT polarity table for a curator to check. It does not decide
#' anything on its own, and it is deliberately unwilling to guess.
#'
#' What each resource can and cannot do is worth knowing before you rely on it.
#' \code{go} is the only one that gives a gene-to-process sign directly, through
#' \code{positive/negative regulation of X} terms, so it needs no path chaining -
#' but its coverage is uneven gene to gene. \code{signor} carries explicit
#' gene-to-phenotype rows, with a cell-biological vocabulary that fits metabolic
#' and immune systems better than behavioural ones. \code{omnipath} has the
#' broadest coverage and the weakest signs: CollecTRI and DoRothEA, two of its
#' inputs, both document assigning an activating mode BY DEFAULT when they do not
#' know, so a \code{consensus_stimulation} at \code{curation_effort = 1} should be
#' read as sign-unknown. \code{trrust} is the most honest transcription-factor
#' resource because it emits \code{Unknown} rather than guessing.
#'
#' Nothing returned here is a polarity. It is a proposal, and the
#' \code{needs_review} column is \code{TRUE} for every row.
#'
#' @param sets A \code{dmsa_sets} or \code{dmsa_selection} whose genes to draft
#'   signs for, or a character vector of gene symbols.
#' @param anchors Named list mapping \code{system_id} to the anchor genes that
#'   define that system's activation tone. Required for anything except
#'   \code{sources = "go"}: without an anchor there is no tone for a sign to
#'   point at.
#' @param sources Which resources to query, in order. See details.
#' @param max_genes Stop after this many genes. Drafting is one HTTP call per
#'   gene for \code{go} and \code{omnipath}, so raise it deliberately.
#' @param quiet Suppress progress messages.
#' @return A data.frame of drafted rows with \code{w_g}, \code{w_g_source},
#'   \code{evidence}, \code{citation} and \code{needs_review}, plus an
#'   \code{attr(, "unreachable")} listing resources that did not respond.
#' @examples
#' \donttest{
#' ## Needs network access: each gene is one HTTP call to the chosen resources.
#' sel <- dmsa_select(systems = "hpa")
#' draft <- dmsa_polarity_fetch(sel, anchors = list("2" = c("CRH", "POMC")))
#' utils::write.csv(draft, file.path(tempdir(), "polarity_draft.csv"),
#'                  row.names = FALSE)
#' }
#' @export
dmsa_polarity_fetch <- function(sets, anchors = NULL,
                                sources = c("go", "omnipath"),
                                max_genes = 200L, quiet = FALSE) {
  sources <- intersect(sources, names(.POL_SOURCES))
  if (!length(sources))
    stop("`sources` must name at least one of: ",
         paste(names(.POL_SOURCES), collapse = ", "), call. = FALSE)
  if (inherits(sets, c("dmsa_sets", "dmsa_selection"))) {
    gs <- unique(sets$cascade[, c("system_id", "gene")])
  } else {
    gs <- data.frame(system_id = NA_character_, gene = unique(as.character(sets)),
                     stringsAsFactors = FALSE)
  }
  if (nrow(gs) > max_genes) {
    if (!quiet)
      message("drafting the first ", max_genes, " of ", nrow(gs),
              " gene-system pairs; raise `max_genes` to go further")
    gs <- gs[seq_len(max_genes), , drop = FALSE]
  }
  if (is.null(anchors) && !identical(sources, "go"))
    warning("no `anchors` given: without the genes that define each system's ",
            "activation tone, a pairwise sign cannot be turned into a ",
            "gene-to-system sign. Only `go` results will be usable.",
            call. = FALSE)

  unreachable <- character(0)
  out <- vector("list", nrow(gs))
  for (i in seq_len(nrow(gs))) {
    sym <- gs$gene[i]; sid <- gs$system_id[i]
    row <- data.frame(system_id = sid, gene = sym, w_g = 0,
                      role = "unresolved", confidence = "low",
                      w_g_source = "unresolved", anchor = FALSE,
                      evidence = "no resource returned a usable sign",
                      citation = "", needs_review = TRUE,
                      stringsAsFactors = FALSE)
    for (s in setdiff(sources, unreachable)) {
      got <- .pol_query(s, sym, anchors[[as.character(sid)]])
      if (is.null(got)) { unreachable <- unique(c(unreachable, s)); next }
      if (is.na(got$w_g)) next
      row$w_g <- got$w_g; row$w_g_source <- s; row$role <- got$role
      row$confidence <- got$confidence; row$evidence <- got$evidence
      row$citation <- got$citation
      break
    }
    out[[i]] <- row
    if (!quiet && i %% 25 == 0) message("  drafted ", i, "/", nrow(gs))
  }
  res <- do.call(rbind, out)
  attr(res, "unreachable") <- unreachable
  if (length(unreachable) && !quiet)
    message("no response from: ", paste(unreachable, collapse = ", "),
            " - those rows are unresolved rather than guessed")
  if (!quiet)
    message(sum(res$w_g_source != "unresolved"), " of ", nrow(res),
            " drafted with a sign; every row needs review before use")
  res
}

## One resource, one gene. Returns NULL when the resource itself is unreachable
## (so the caller can stop asking it) and w_g = NA when it simply has nothing.
.pol_query <- function(source, sym, anchor) {
  spec <- .POL_SOURCES[[source]]
  txt <- try(suppressWarnings(
    readLines(url(spec$url(sym), open = "rb"), warn = FALSE)), silent = TRUE)
  if (inherits(txt, "try-error") || !length(txt)) return(NULL)
  txt <- paste(txt, collapse = "\n")
  none <- list(w_g = NA_real_, role = NA, confidence = NA, evidence = "",
               citation = "")
  if (source == "go") {
    pos <- grepl("positive regulation of", txt, fixed = TRUE)
    neg <- grepl("negative regulation of", txt, fixed = TRUE)
    if (pos && !neg) return(list(w_g = 1, role = "driver", confidence = "medium",
      evidence = "GO: positive regulation term(s) only", citation = "GO API"))
    if (neg && !pos) return(list(w_g = -1, role = "brake", confidence = "medium",
      evidence = "GO: negative regulation term(s) only", citation = "GO API"))
    if (pos && neg) return(list(w_g = 0, role = "ambiguous", confidence = "low",
      evidence = "GO carries both positive and negative regulation terms",
      citation = "GO API"))
    return(none)
  }
  if (source == "omnipath") {
    if (is.null(anchor) || !length(anchor)) return(none)
    ln <- strsplit(txt, "\n")[[1]]
    hit <- ln[vapply(ln, function(z) any(vapply(anchor, function(a)
      grepl(paste0("\t", a, "\t"), z, fixed = TRUE), TRUE)), TRUE)]
    if (!length(hit)) return(none)
    stim <- sum(grepl("\tTrue\tFalse\t", hit, fixed = TRUE))
    inh  <- sum(grepl("\tFalse\tTrue\t", hit, fixed = TRUE))
    if (stim == inh) return(list(w_g = 0, role = "ambiguous", confidence = "low",
      evidence = "OmniPath: stimulation and inhibition edges to the anchor tie",
      citation = "omnipathdb.org"))
    return(list(w_g = if (stim > inh) 1 else -1,
      role = if (stim > inh) "driver" else "brake", confidence = "low",
      evidence = paste0("OmniPath: ", stim, " stimulation vs ", inh,
                        " inhibition edge(s) to anchor - VERIFY, a default-to-activating",
                        " input can produce this"),
      citation = "omnipathdb.org"))
  }
  none
}

#' Which public resources DMSA can draft polarity from, and what each is good for
#'
#' @return data.frame of resource name, endpoint template and caveat.
#' @examples
#' src <- dmsa_polarity_sources()
#' src[c("source", "resource")]
#' ## the trap worth knowing about: two of OmniPath's inputs record an
#' ## activating sign by default when the direction is actually unknown
#' src$caveat[src$source == "omnipath"]
#' @export
dmsa_polarity_sources <- function() {
  d <- data.frame(
    source = names(.POL_SOURCES),
    resource = vapply(.POL_SOURCES, `[[`, "", "name"),
    caveat = vapply(.POL_SOURCES, `[[`, "", "note"),
    stringsAsFactors = FALSE)
  rownames(d) <- NULL
  d
}
