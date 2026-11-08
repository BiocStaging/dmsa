# dmsa selection cascade: system > module > gene > probe
# --------------------------------------------------------------------------
# The cascade is the DECLARED STRUCTURE of the analysis: which probes belong to
# which gene, which genes to which module, which modules to which system. It is
# deliberately separate from the direction map (which supplies d and p_plus per
# probe) and from the polarity table (which supplies w_g per gene), because the
# three answer different questions and are curated by different evidence.
#
# Why this matters for inference, not just convenience: because families are
# declared here and never harvested from whatever the analyst happened to load,
# a unit's family is a property of the map rather than of the session. That is
# the mechanism behind DMSA's invariance to set selection - see dmsa_pdance().
#
# The bundled cascade is Project Alpha 2026c, literature-audited at module level
# (30 systems, 188 modules, 1,234 genes, 16,823 CpGs). Any panel in the same
# shape can replace it; dmsa_sets_template() writes the schema and
# dmsa_sets_check() validates a candidate file before it is used.

## The BIOLOGICAL cascade is system > module > gene. A CpG/probe column is
## OPTIONAL provenance, not part of the biological definition: which CpGs
## exist is a property of the user's data, not of the biology. Cascades that
## still carry `cpg` keep working and gain the CpG counts below.
.CAS_REQ <- c("system_id", "system", "module_id", "module", "gene")
.CAS_OPT_CPG <- c("cpg", "probe_id")
.CAS_ID  <- c("system_id", "module_id")          # never parse these as numeric
.CAS_COLPREFIX <- "col_"

## ---- short names ---------------------------------------------------------
## Auto-derived when a cascade carries no `system_short`: lower case, first
## token before a comma/slash/ampersand, non-alphanumerics to underscore.
## Deterministic, and disambiguated by suffix if two systems collide.
.cas_slug <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s <- sub("[,/(].*$", "", s)
  s <- sub("\\s*&.*$", "", s)
  s <- gsub("[^a-z0-9]+", "_", s)
  s <- gsub("^_+|_+$", "", s)
  s[!nzchar(s)] <- "set"
  d <- duplicated(s) | duplicated(s, fromLast = TRUE)
  if (any(d)) s[d] <- paste0(s[d], "_", seq_len(sum(d)))
  s
}

.cas_read <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else path
  hdr <- utils::read.csv(con, nrows = 1L, stringsAsFactors = FALSE,
                         check.names = FALSE)
  con2 <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else path
  want <- c("system_id", "module_id", "cpg", "gene", "probe_id",
            "system_short", "system", "module")
  cc <- stats::setNames(rep("character", length(intersect(want, names(hdr)))),
                        intersect(want, names(hdr)))
  utils::read.csv(con2, stringsAsFactors = FALSE, colClasses = cc,
                  check.names = FALSE)
}

.cas_builtin <- function(what = c("cascade", "audit")) {
  what <- match.arg(what)
  f <- if (what == "cascade")
    system.file("extdata", "alpha_cascade_2026c.csv.gz", package = "dmsa")
  else
    system.file("extdata", "alpha_module_audit_2026c.csv", package = "dmsa")
  if (!nzchar(f) || !file.exists(f)) NULL else f
}

#' Load a selection cascade (system > module > gene > probe)
#'
#' The cascade declares the analysis structure DMSA corrects within. Loading it
#' is separate from loading a direction map: the cascade says which probes form
#' which unit, the direction map says which way each probe points.
#'
#' @param x Either \code{"alpha"} for the bundled, literature-audited Project
#'   Alpha 2026c cascade, a path to a CSV (optionally \code{.gz}) in the same
#'   shape, or a \code{data.frame} already in that shape. See
#'   \code{dmsa_sets_template()} for the schema and
#'   \code{dmsa_sets_check()} to validate a file before use.
#' @param audit Optional module-level evidence table: a path or
#'   \code{data.frame} with one row per \code{module_id} and any of
#'   \code{evidence_strength}, \code{audit_status}, \code{citation_keys},
#'   \code{evidence_note}, \code{deep_search_url}. Defaults to the bundled audit
#'   when \code{x = "alpha"}, and to \code{NULL} otherwise. When a user cascade
#'   already carries these columns they are used and no separate table is
#'   needed.
#' @param polarity Optional gene-to-system polarity table (path,
#'   \code{data.frame} or \code{dmsa_polarity}) supplying \code{w_g}. Defaults
#'   to the bundled Alpha polarity when \code{x = "alpha"}, to the cascade's own
#'   \code{w_g} column when it has one, and otherwise to none - in which case a
#'   system-level score weights every gene +1 and the print method says so. Use
#'   \code{dmsa_polarity_fetch()} to draft one for your own panel from public
#'   databases.
#' @param name Optional label, printed and carried into reports.
#' @return An object of class \code{dmsa_sets}.
#' @seealso \code{dmsa_systems()} for the short names,
#'   \code{dmsa_select()} to choose systems, \code{dmsa_evidence()} for the
#'   module evidence table.
#' @examples
#' cas <- dmsa_sets()          # bundled Alpha 2026c
#' cas
#' dmsa_systems()                 # the short names to use in `systems =`
#' @export
dmsa_sets <- function(x = "alpha", audit = NULL, polarity = NULL,
                      name = NULL) {
  src <- NULL
  is_alpha <- is.character(x) && length(x) == 1L && identical(x, "alpha")
  if (is_alpha) {
    f <- .cas_builtin("cascade")
    if (is.null(f))
      stop("the bundled Alpha cascade is not installed with this copy of dmsa; ",
           "pass a cascade CSV instead (see dmsa_sets_template())",
           call. = FALSE)
    cas <- .cas_read(f); src <- f
    if (is.null(audit)) audit <- .cas_builtin("audit")
    if (is.null(name)) name <- "Project Alpha 2026c (module-audited)"
  } else if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("cascade file not found: ", x, call. = FALSE)
    cas <- .cas_read(x); src <- x
    if (is.null(name)) name <- basename(x)
  } else {
    cas <- as.data.frame(x, stringsAsFactors = FALSE)
    if (is.null(name)) name <- "user data.frame"
  }

  miss <- setdiff(.CAS_REQ, names(cas))
  if (length(miss))
    stop("cascade is missing required column(s): ", paste(miss, collapse = ", "),
         "\nrun dmsa_sets_template() to see the schema", call. = FALSE)
  for (v in intersect(c(.CAS_ID, "gene", "cpg", "probe_id"), names(cas)))
    cas[[v]] <- as.character(cas[[v]])
  if (!"probe_id" %in% names(cas) && "cpg" %in% names(cas))
    cas$probe_id <- cas$cpg
  ## and the mirror: a cascade declaring only probe_id still HAS a probe layer.
  ## Every downstream "does this cascade carry CpGs" test keys on `cpg`, so
  ## derive it (canonical cg id = the probe id up to the first underscore)
  ## rather than teaching each of those sites about both spellings.
  if (!"cpg" %in% names(cas) && "probe_id" %in% names(cas))
    cas$cpg <- sub("_.*$", "", tolower(as.character(cas$probe_id)))
  if (!"system_short" %in% names(cas)) {
    tab <- unique(cas[, c("system_id", "system")])
    tab <- tab[order(.cas_num(tab$system_id), tab$system), , drop = FALSE]
    tab$system_short <- .cas_slug(tab$system)
    cas$system_short <- tab$system_short[match(cas$system_id, tab$system_id)]
  }
  cas$system_short <- tolower(as.character(cas$system_short))

  ## ---- module-level evidence ---------------------------------------------
  ev_cols <- c("evidence_strength", "audit_status", "citation_keys",
               "evidence_note", "deep_search_url")
  ad <- NULL
  if (!is.null(audit)) {
    ad <- if (is.character(audit) && length(audit) == 1L) {
      if (!file.exists(audit)) stop("audit file not found: ", audit, call. = FALSE)
      .cas_read(audit)
    } else as.data.frame(audit, stringsAsFactors = FALSE)
    if (!"module_id" %in% names(ad))
      stop("audit table needs a `module_id` column", call. = FALSE)
    ad$module_id <- as.character(ad$module_id)
  }
  ## a cascade that already carries the evidence columns wins over a table
  inline <- intersect(ev_cols, names(cas))
  mods <- unique(cas[, c("system_id", "system_short", "system",
                         "module_id", "module")])
  for (v in ev_cols) {
    mods[[v]] <- if (v %in% inline)
      cas[[v]][match(mods$module_id, cas$module_id)]
    else if (!is.null(ad) && v %in% names(ad))
      ad[[v]][match(mods$module_id, ad$module_id)]
    else NA_character_
  }
  mods$n_genes <- as.integer(tapply(cas$gene, cas$module_id,
                                    function(g) length(unique(g)))[mods$module_id])
  mods$n_cpgs  <- if ("cpg" %in% names(cas))
    as.integer(tapply(cas$cpg, cas$module_id,
                      function(g) length(unique(g)))[mods$module_id])
  else NA_integer_
  mods <- mods[order(.cas_num(mods$system_id), .cas_modnum(mods$module_id)), ,
               drop = FALSE]
  rownames(mods) <- NULL

  sysd <- unique(cas[, c("system_id", "system_short", "system")])
  sysd$n_modules <- as.integer(tapply(cas$module_id, cas$system_id,
                                      function(m) length(unique(m)))[sysd$system_id])
  sysd$n_genes <- as.integer(tapply(cas$gene, cas$system_id,
                                    function(g) length(unique(g)))[sysd$system_id])
  sysd$n_cpgs <- if ("cpg" %in% names(cas))
    as.integer(tapply(cas$cpg, cas$system_id,
                      function(g) length(unique(g)))[sysd$system_id])
  else NA_integer_
  sysd <- sysd[order(.cas_num(sysd$system_id)), , drop = FALSE]
  rownames(sysd) <- NULL

  ## Polarity travels with the cascade: the bundled Alpha table for the bundled
  ## cascade, and for a user cascade whatever polarity columns it carries. A
  ## cascade with no polarity is usable - system scores then weight every gene
  ## +1 - but the print method says so rather than letting it pass unnoticed.
  ## spec 41: POLARITY THE USER SUPPLIED IS NEVER SILENTLY DISCARDED.
  ## try()/inherits("try-error") -> NULL turned a malformed polarity table into
  ## no polarity at all, and a cascade with no polarity weights every gene +1 -
  ## so a typo in a w_g column became an unannounced all-activating system
  ## score. Anything the USER provided now errors and names the problem; only
  ## the bundled table may degrade, and then it warns.
  pol <- NULL
  if (!is.null(polarity)) {
    pol <- tryCatch(dmsa_polarity(polarity), error = function(e)
      stop("the `polarity` you supplied could not be read: ",
           conditionMessage(e),
           "\nIt is not being ignored: a cascade with no polarity weights ",
           "every gene +1, which would silently make every system look purely ",
           "activating.\nFix the table, or pass polarity = NULL to declare ",
           "that you have none.\nNo cascade was built.", call. = FALSE))
    ## E13 (2026-09-02): on the BUNDLED Alpha cascade a user table MERGES
    ## with the bundled curation - user rows override gene for gene, the
    ## curation fills every gene the user did not name - exactly as
    ## dmsa_align(polarity =) has done since E4. Replacing the whole table
    ## silently left every unnamed gene without polarity, and the report
    ## then weighted those genes 0 at the system level with nothing said.
    ## A user cascade keeps the old meaning: the table IS its polarity.
    if (isTRUE(is_alpha) && !inherits(polarity, "dmsa_polarity")) {
      bnd <- tryCatch(dmsa_polarity("alpha"), error = function(e) NULL)
      if (!is.null(bnd)) {
        usr <- pol$polarity; bp <- bnd$polarity
        usr$w_g_source[is.na(usr$w_g_source) | usr$w_g_source == "unstated"] <- "user"
        usr$grade <- "user"
        key_u <- paste(usr$system_id, usr$gene)
        fill <- bp[!paste(bp$system_id, bp$gene) %in% key_u, , drop = FALSE]
        cols <- union(names(bp), names(usr))
        for (cn in setdiff(cols, names(usr))) usr[[cn]] <- NA
        for (cn in setdiff(cols, names(fill))) fill[[cn]] <- NA
        merged <- rbind(usr[cols], fill[cols])
        merged <- merged[order(merged$system_id, merged$gene), , drop = FALSE]
        rownames(merged) <- NULL
        pol$polarity <- merged
        pol$name <- sprintf("%s + user overrides (%d gene(s))", bnd$name, nrow(usr))
        message(sprintf(paste0("dmsa: polarity merged - %d gene(s) from your ",
                               "table override the bundled Alpha curation, ",
                               "%d filled from it."), nrow(usr), nrow(fill)))
      }
    }
  } else if (isTRUE(is_alpha)) {
    ## keyed on WHERE the cascade came from, not on the display `name` - a
    ## user renaming the bundled cascade must not silently lose its polarity
    pol <- tryCatch(dmsa_polarity("alpha"), error = function(e) {
      warning("the bundled Alpha polarity table could not be read (",
              conditionMessage(e), "); this cascade carries no polarity, so ",
              "signed system-level analysis is unavailable")
      NULL })
  } else if ("w_g" %in% names(cas)) {
    keep <- intersect(c("system_id", "system_short", "system", "module_id",
                        "gene", "w_g", "role", "confidence", "w_g_source",
                        "anchor", "evidence", "citation"), names(cas))
    pol <- tryCatch(dmsa_polarity(unique(cas[, keep])), error = function(e)
      stop("this cascade carries a `w_g` column, but it could not be read as ",
           "polarity: ", conditionMessage(e),
           "\nIgnoring it would weight every gene +1 and report a purely ",
           "activating system score without saying so.\nFix the w_g column, ",
           "or remove it to declare that the cascade has no polarity.\n",
           "No cascade was built.", call. = FALSE))
  }

  structure(list(cascade = cas, systems = sysd, modules = mods,
                 polarity = pol, source = src, name = name,
                 columns = grep(paste0("^", .CAS_COLPREFIX), names(cas),
                                value = TRUE)),
            class = "dmsa_sets")
}

.cas_num <- function(x) suppressWarnings(as.numeric(x))
## sort 24.10 after 24.9, not after 24.1
.cas_modnum <- function(x) {
  p <- strsplit(as.character(x), ".", fixed = TRUE)
  a <- suppressWarnings(as.numeric(vapply(p, function(z) z[1], "")))
  b <- suppressWarnings(as.numeric(vapply(p, function(z)
    if (length(z) > 1) z[2] else "0", "")))
  a * 1000 + ifelse(is.na(b), 0, b)
}

#' @export
print.dmsa_sets <- function(x, ...) {
  cat("dmsa selection cascade:", x$name, "\n")
  .ncpg <- if ("cpg" %in% names(x$cascade))
    length(unique(x$cascade$cpg)) else NA_integer_
  if (is.na(.ncpg))
    cat(sprintf("  %d systems | %d modules | %d genes | %d rows (biological; no CpG column)\n",
                nrow(x$systems), nrow(x$modules),
                length(unique(x$cascade$gene)), nrow(x$cascade)))
  else
    cat(sprintf("  %d systems | %d modules | %d genes | %d CpGs | %d rows\n",
                nrow(x$systems), nrow(x$modules),
                length(unique(x$cascade$gene)), .ncpg, nrow(x$cascade)))
  if (length(x$columns))
    cat("  data-column keys available:", paste(x$columns, collapse = ", "), "\n")
  cat("\n  use these short names in systems = c(...):\n")
  s <- x$systems
  w <- max(nchar(s$system_short))
  for (i in seq_len(nrow(s)))
    if (is.na(s$n_cpgs[i]))
      cat(sprintf("   %3s  %-*s  %-44s  %2d mod  %4d gene\n",
                  s$system_id[i], w, s$system_short[i],
                  substr(s$system[i], 1, 44), s$n_modules[i], s$n_genes[i]))
    else
      cat(sprintf("   %3s  %-*s  %-44s  %2d mod  %4d gene  %5d cpg\n",
                  s$system_id[i], w, s$system_short[i],
                  substr(s$system[i], 1, 44), s$n_modules[i], s$n_genes[i],
                  s$n_cpgs[i]))
  .cas_evidence_banner(x$modules, indent = "  ")
  .cas_polarity_banner(x, indent = "  ")
  cat("\n  dmsa_select(systems = c(\"", s$system_short[1], "\")) ",
      "-> modules, genes and probes default to \"full\"\n", sep = "")
  invisible(x)
}

## The evidence banner is printed wherever a cascade or selection is shown.
## Module labels in the bundled cascade were literature-audited; a user who
## reports a module-level finding should be able to see, without asking, how
## strong the evidence for that module's DEFINITION is and what it rests on.
.cas_evidence_banner <- function(mods, indent = "") {
  es <- mods$evidence_strength
  if (all(is.na(es))) {
    cat(indent, "module evidence: not annotated in this cascade\n", sep = "")
    return(invisible(NULL))
  }
  n_hi <- sum(es %in% "High", na.rm = TRUE)
  n_mo <- sum(es %in% "Moderate", na.rm = TRUE)
  cat(sprintf("%smodule evidence: %d High, %d Moderate", indent, n_hi, n_mo))
  st <- mods$audit_status
  flag <- !is.na(st) & (grepl("heterogeneous|measurement_defined", st))
  if (any(flag)) cat(sprintf(", %d flagged (heterogeneous or measurement-defined)",
                             sum(flag)))
  cat("\n")
  if (n_mo || any(flag))
    cat(indent, "  dmsa_evidence() lists them with their citation keys\n", sep = "")
  invisible(NULL)
}

#' The systems available in a cascade, with the short names to select them by
#'
#' @param x A \code{dmsa_sets} or \code{dmsa_selection}, or anything
#'   \code{dmsa_sets()} accepts. A bare string that is neither \code{"alpha"}
#'   nor an existing file path is treated as \code{pattern}, so
#'   \code{dmsa_systems("hpa")} searches the bundled cascade.
#' @param pattern Optional case-insensitive regular expression; matched against
#'   the short name, the full system name and the system id.
#' @return A data.frame with \code{system_id}, \code{system_short},
#'   \code{system} and the module/gene/CpG counts. Printed as a table.
#' @examples
#' dmsa_systems()             # all 30
#' dmsa_systems("stress|hpa") # search
#' @export
dmsa_systems <- function(x = "alpha", pattern = NULL) {
  ## dmsa_systems("hpa") should search, not hunt for a file called "hpa". A bare
  ## string that is neither "alpha" nor an existing path is a pattern - which is
  ## what any reader of the examples would assume, and what they now get.
  if (is.null(pattern) && is.character(x) && length(x) == 1L &&
      !identical(x, "alpha") && !file.exists(x)) {
    pattern <- x; x <- "alpha"
  }
  cas <- if (inherits(x, c("dmsa_sets", "dmsa_selection"))) x else dmsa_sets(x)
  s <- cas$systems
  if (!is.null(pattern)) {
    hit <- grepl(pattern, s$system_short, ignore.case = TRUE) |
           grepl(pattern, s$system, ignore.case = TRUE) |
           grepl(pattern, s$system_id, ignore.case = TRUE)
    if (!any(hit)) {
      message("no system matches '", pattern, "'; all short names: ",
              paste(s$system_short, collapse = ", "))
      return(invisible(s[0, , drop = FALSE]))
    }
    s <- s[hit, , drop = FALSE]
  }
  s
}

## ---- resolver ------------------------------------------------------------
## Accepts, in this order: exact short name, exact system id, exact full name,
## unique short-name prefix, unique substring of the full name. Case-insensitive
## throughout. Ambiguity and misses both error with the candidates named, because
## silently selecting the wrong system is the one failure mode that would not
## announce itself downstream.
.cas_resolve_systems <- function(cas, systems) {
  s <- cas$systems
  if (is.null(systems)) return(s)
  q <- as.character(systems)
  out <- integer(0); bad <- character(0)
  for (k in q) {
    kk <- tolower(trimws(k))
    i <- which(tolower(s$system_short) == kk)
    if (!length(i)) i <- which(tolower(s$system_id) == kk)
    if (!length(i)) i <- which(tolower(s$system) == kk)
    if (!length(i)) i <- which(startsWith(tolower(s$system_short), kk))
    if (!length(i)) i <- grep(kk, tolower(s$system), fixed = TRUE)
    if (!length(i)) { bad <- c(bad, k); next }
    if (length(i) > 1L)
      stop("`systems = \"", k, "\"` is ambiguous - it matches ",
           paste(sprintf("'%s'", s$system_short[i]), collapse = ", "),
           "; use one of those short names", call. = FALSE)
    out <- c(out, i)
  }
  if (length(bad))
    stop("`systems` did not match: ", paste(sprintf("'%s'", bad), collapse = ", "),
         "\navailable short names: ", paste(s$system_short, collapse = ", "),
         "\n(dmsa_systems() prints them with the full names)", call. = FALSE)
  s[unique(out), , drop = FALSE]
}

.cas_resolve_below <- function(vals, ask, what, extra = NULL) {
  if (is.null(ask) || (length(ask) == 1L && identical(tolower(as.character(ask)),
                                                      "full")))
    return(list(all = TRUE, keep = vals))
  ask <- as.character(ask)
  pool <- unique(c(vals, extra))
  kk <- tolower(ask); vv <- tolower(pool)
  i <- match(kk, vv)
  bad <- ask[is.na(i)]
  if (length(bad))
    stop("`", what, "` did not match within the selected system(s): ",
         paste(sprintf("'%s'", bad), collapse = ", "),
         "\nfirst few available: ", paste(utils::head(sort(vals), 8),
                                          collapse = ", "),
         if (length(vals) > 8) " ..." else "", call. = FALSE)
  list(all = FALSE, keep = pool[i])
}

#' Select a set of units from a cascade
#'
#' Choose one or more systems by short name. Everything below the system -
#' modules, genes, probes - defaults to \code{"full"}, meaning every unit the
#' cascade assigns to the chosen system(s). Narrow a level only when the
#' question is genuinely narrower, and note that narrowing changes the declared
#' family and therefore the multiplicity toll.
#'
#' @param x A \code{dmsa_sets}, or anything \code{dmsa_sets()} accepts.
#' @param systems Short names (\code{c("hpa", "oxytocin")}), system ids, or full
#'   system names. \code{NULL} selects every system in the cascade. Matching is
#'   case-insensitive and accepts a unique prefix of a short name.
#' @param modules,genes,probes \code{"full"} (default) or explicit ids/symbols
#'   to restrict to. \code{modules} takes module ids such as \code{"2.6"},
#'   \code{genes} takes symbols, \code{probes} takes CpG ids or probe ids.
#' @param columns Which cascade column holds the data-matrix column names for
#'   this analysis - e.g. \code{"col_parent_T1"}. \code{"auto"} picks the
#'   \code{col_*} column with the most overlap against \code{data}, when
#'   \code{data} is supplied. \code{NULL} leaves it unresolved.
#' @param data Optional data.frame or matrix used by \code{columns = "auto"}.
#' @return An object of class \code{dmsa_selection}.
#' @examples
#' sel <- dmsa_select(systems = c("hpa", "oxytocin"))
#' sel
#' dmsa_select(systems = "hpa", modules = c("2.6", "2.8"))
#' @export
dmsa_select <- function(x = "alpha", systems = NULL, modules = "full",
                        genes = "full", probes = "full", columns = NULL,
                        data = NULL) {
  cas <- if (inherits(x, "dmsa_sets")) x else dmsa_sets(x)
  sysd <- .cas_resolve_systems(cas, systems)
  d <- cas$cascade[cas$cascade$system_id %in% sysd$system_id, , drop = FALSE]

  m <- .cas_resolve_below(unique(d$module_id), modules, "modules")
  d <- d[d$module_id %in% m$keep, , drop = FALSE]
  g <- .cas_resolve_below(unique(d$gene), genes, "genes")
  d <- d[d$gene %in% g$keep, , drop = FALSE]
  ## A BIOLOGICAL cascade (system > module > gene) carries no CpG column, so
  ## there is nothing to restrict at probe level; asking for one is a user
  ## error worth naming rather than silently ignoring.
  p <- list(keep = character(), all = TRUE)   # no probe layer by default
  if ("cpg" %in% names(d)) {
    p <- .cas_resolve_below(unique(d$cpg), probes, "probes",
                            extra = unique(d$probe_id))
    d <- d[d$cpg %in% p$keep | d$probe_id %in% p$keep, , drop = FALSE]
  } else if (!is.null(probes) && !identical(probes, "full")) {
    stop("`probes` was supplied but this cascade is biological ",
         "(system > module > gene) and carries no CpG column. Which CpGs ",
         "exist is a property of your methylation data, not of the ",
         "reference; restrict at gene level instead.", call. = FALSE)
  }
  if (!nrow(d)) stop("selection is empty after applying the restrictions",
                     call. = FALSE)

  colkey <- NULL
  if (!is.null(columns)) {
    if (identical(columns, "auto")) {
      if (is.null(data))
        stop("columns = 'auto' needs `data` to match against", call. = FALSE)
      nm <- if (is.matrix(data)) colnames(data) else names(data)
      cand <- cas$columns
      if (!length(cand))
        stop("this cascade carries no col_* column to match", call. = FALSE)
      hits <- vapply(cand, function(cc) sum(d[[cc]] %in% nm), 0L)
      if (max(hits) == 0L)
        stop("no col_* column in the cascade matches any column of `data`",
             call. = FALSE)
      colkey <- names(which.max(hits))
    } else {
      if (!columns %in% names(d))
        stop("`columns = \"", columns, "\"` is not a column of the cascade; ",
             "available: ", paste(cas$columns, collapse = ", "), call. = FALSE)
      colkey <- columns
    }
  }

  mods <- cas$modules[cas$modules$module_id %in% unique(d$module_id), ,
                      drop = FALSE]
  ## recount within the selection, so a narrowed gene set is reported honestly
  mods$n_genes_selected <- as.integer(tapply(d$gene, d$module_id,
      function(z) length(unique(z)))[mods$module_id])
  mods$n_cpgs_selected <- if ("cpg" %in% names(d))
    as.integer(tapply(d$cpg, d$module_id,
        function(z) length(unique(z)))[mods$module_id])
  else NA_integer_

  pol <- cas$polarity
  if (!is.null(pol)) {
    pp <- pol$polarity
    pol$polarity <- pp[pp$system_id %in% sysd$system_id &
                         pp$gene %in% unique(d$gene), , drop = FALSE]
    pol$coverage <- .pol_coverage(pol$polarity, list(cascade = d))
  }

  structure(list(
    cascade = d, systems = sysd, modules = mods, polarity = pol,
    column_key = colkey,
    spec = list(systems = systems, modules = modules, genes = genes,
                probes = probes),
    full = c(modules = m$all, genes = g$all, probes = p$all),
    source = cas$source, name = cas$name
  ), class = "dmsa_selection")
}

#' @export
print.dmsa_selection <- function(x, ...) {
  d <- x$cascade
  cat("dmsa selection from:", x$name, "\n")
  lv <- c(modules = if (x$full["modules"]) "full" else "restricted",
          genes   = if (x$full["genes"])   "full" else "restricted",
          probes  = if (x$full["probes"])  "full" else "restricted")
  cat(sprintf("  %d system(s) selected | modules %s | genes %s | probes %s\n",
              nrow(x$systems), lv[1], lv[2], lv[3]))
  for (i in seq_len(nrow(x$systems))) {
    sid <- x$systems$system_id[i]
    dd <- d[d$system_id == sid, , drop = FALSE]
    if ("cpg" %in% names(dd))
      cat(sprintf("   %-16s %-40s %2d mod  %4d gene  %5d cpg\n",
                  x$systems$system_short[i], substr(x$systems$system[i], 1, 40),
                  length(unique(dd$module_id)), length(unique(dd$gene)),
                  length(unique(dd$cpg))))
    else
      cat(sprintf("   %-16s %-40s %2d mod  %4d gene\n",
                  x$systems$system_short[i], substr(x$systems$system[i], 1, 40),
                  length(unique(dd$module_id)), length(unique(dd$gene))))
  }
  if ("cpg" %in% names(d))
    cat(sprintf("  total: %d modules, %d genes, %d CpGs\n",
                length(unique(d$module_id)), length(unique(d$gene)),
                length(unique(d$cpg))))
  else
    cat(sprintf("  total: %d modules, %d genes\n",
                length(unique(d$module_id)), length(unique(d$gene))))
  if (!is.null(x$column_key))
    cat("  data columns keyed by:", x$column_key, "\n")
  .cas_evidence_banner(x$modules, indent = "  ")
  .cas_polarity_banner(x, indent = "  ")
  mo <- x$modules[x$modules$evidence_strength %in% "Moderate" |
                    grepl("heterogeneous|measurement_defined",
                          x$modules$audit_status), , drop = FALSE]
  if (nrow(mo)) {
    cat("  modules to report with their evidence caveat:\n")
    for (i in seq_len(min(nrow(mo), 8L)))
      cat(sprintf("   %-7s %-44s %-8s %s\n", mo$module_id[i],
                  substr(mo$module[i], 1, 44),
                  ifelse(is.na(mo$evidence_strength[i]), "-",
                         mo$evidence_strength[i]),
                  ifelse(is.na(mo$citation_keys[i]), "", mo$citation_keys[i])))
    if (nrow(mo) > 8L) cat(sprintf("   ... and %d more (dmsa_evidence())\n",
                                   nrow(mo) - 8L))
  }
  invisible(x)
}

#' Module-level evidence and citations for a cascade, selection or frame
#'
#' Every module label in the bundled cascade was checked against the literature.
#' This returns what that check concluded, so a module-level finding can be
#' reported with the strength of its own definition attached: whether the label
#' was retained or revised, whether the locked membership is homogeneous, the
#' evidence tier, and the citation keys behind it.
#'
#' @param x A \code{dmsa_sets}, \code{dmsa_selection}, or \code{dmsa_frame}.
#' @param which \code{"all"}, \code{"moderate"} (Moderate evidence only), or
#'   \code{"flagged"} (Moderate, heterogeneous, or measurement-defined).
#' @return data.frame of class \code{dmsa_evidence}, one row per module.
#' @examples
#' dmsa_evidence(dmsa_select(systems = "immune"), which = "flagged")
#' @export
dmsa_evidence <- function(x, which = c("all", "moderate", "flagged")) {
  which <- match.arg(which)
  mods <- if (inherits(x, c("dmsa_sets", "dmsa_selection"))) x$modules
          else if (!is.null(x$module_evidence)) x$module_evidence
          else stop("no module evidence attached to this object", call. = FALSE)
  if (which == "moderate")
    mods <- mods[mods$evidence_strength %in% "Moderate", , drop = FALSE]
  if (which == "flagged")
    mods <- mods[mods$evidence_strength %in% "Moderate" |
                   grepl("heterogeneous|measurement_defined",
                         mods$audit_status), , drop = FALSE]
  structure(mods, class = c("dmsa_evidence", "data.frame"))
}

#' @export
print.dmsa_evidence <- function(x, n = 25, ...) {
  if (!nrow(x)) { cat("no modules with an evidence caveat\n"); return(invisible(x)) }
  cat(sprintf("module evidence (%d module%s)\n", nrow(x),
              if (nrow(x) == 1L) "" else "s"))
  k <- min(nrow(x), n)
  for (i in seq_len(k)) {
    cat(sprintf(" %-7s %-8s %-42s %s\n", x$module_id[i],
                ifelse(is.na(x$evidence_strength[i]), "-", x$evidence_strength[i]),
                substr(x$module[i], 1, 42),
                ifelse(is.na(x$citation_keys[i]), "", x$citation_keys[i])))
    st <- x$audit_status[i]
    if (!is.na(st) && !identical(st, "supported_as_named"))
      cat("          status:", st, "\n")
  }
  if (nrow(x) > k) cat(sprintf(" ... %d more rows\n", nrow(x) - k))
  u <- unique(stats::na.omit(x$deep_search_url))
  if (length(u)) cat(" evidence searches:\n", paste("  ", u, collapse = "\n"), "\n")
  invisible(x)
}

## ---- user cascades -------------------------------------------------------

#' Write a cascade template and the instructions for building one
#'
#' @param path Where to write the template CSV. The instructions are written
#'   alongside it as a \code{.md} file with the same stem.
#' @param example Include two filled example rows. Default \code{TRUE}.
#' @return The paths written, invisibly.
#' @examples
#' p <- tempfile(fileext = ".csv")
#' dmsa_sets_template(p)
#'
#' ## the .md alongside it states the rules the validator enforces
#' basename(c(p, sub("\\.csv$", ".md", p)))
#' utils::read.csv(p)[, c("system_short", "module_id", "gene", "cpg")]
#' @export
dmsa_sets_template <- function(path = "dmsa_sets_template.csv",
                                  example = TRUE) {
  cols <- c("system_id", "system_short", "system", "module_id", "module",
            "gene", "cpg", "probe_id", "evidence_strength", "audit_status",
            "citation_keys", "evidence_note")
  tem <- as.data.frame(matrix(character(), 0, length(cols),
                              dimnames = list(NULL, cols)),
                       stringsAsFactors = FALSE)
  if (example)
    tem <- data.frame(
      system_id = c("1", "1"), system_short = c("hpa", "hpa"),
      system = c("HPA axis", "HPA axis"),
      module_id = c("1.1", "1.1"),
      module = c("Corticosteroid receptors", "Corticosteroid receptors"),
      gene = c("NR3C1", "FKBP5"),
      cpg = c("cg01234567", "cg07654321"),
      probe_id = c("cg01234567", "cg07654321"),
      evidence_strength = c("High", "High"),
      audit_status = c("supported_as_named", "supported_as_named"),
      citation_keys = c("[Fri17]; [Fad23]", "[Fri17]; [Fad23]"),
      evidence_note = rep("receptor and chaperone complex; curated", 2),
      stringsAsFactors = FALSE)
  utils::write.csv(tem, path, row.names = FALSE, na = "")
  md <- sub("\\.csv$", ".md", path)
  if (identical(md, path)) md <- paste0(path, ".md")
  writeLines(.cas_instructions(), md)
  message("wrote ", path, " and ", md)
  invisible(c(csv = path, instructions = md))
}

.cas_instructions <- function() c(
"# Supplying your own selection cascade to DMSA",
"",
"DMSA corrects multiplicity inside *declared* families: the genes of a named",
"system, the probes of a named gene. The cascade is where you declare them. It",
"is one long CSV, one row per CpG, and it answers four questions at once:",
"which system, which module inside it, which gene, which CpG.",
"",
"    system  >  module  >  gene  >  probe",
"",
"This file is not the direction map. The cascade says what belongs together;",
"the direction map (`map =`) says which way each probe points, and the polarity",
"table says how each gene relates to its system's activation tone. Keep them",
"separate: they are curated from different evidence and revised on different",
"schedules.",
"",
"## Required columns",
"",
"| column | what it holds |",
"|---|---|",
"| `system_id` | stable id, e.g. `2`. Treated as text - `24.10` must not become `24.1` |",
"| `system` | the printed system name |",
"| `module_id` | stable id nested under the system, e.g. `2.6` |",
"| `module` | the printed module name |",
"| `gene` | gene symbol |",
"| `cpg` | CpG id, e.g. `cg04712664` |",
"",
"## Optional but recommended",
"",
"| column | what it buys you |",
"|---|---|",
"| `system_short` | the handle you type in `systems = c(\"hpa\")`. Derived automatically if absent, but deriving it from a long name gives you something clumsy - set it yourself |",
"| `probe_id` | when one CpG is measured by more than one probe. Defaults to `cpg` |",
"| `col_*` | one column per analysis holding the *data-matrix* column name for that probe, e.g. `col_parent_T1`. Lets `dmsa_select(columns = \"auto\")` key the cascade to your matrix instead of you renaming columns |",
"| `evidence_strength` | `High` or `Moderate` - how well supported the module's *definition* is |",
"| `audit_status` | e.g. `supported_as_named`, `renamed_for_mechanistic_precision`, `heterogeneous_locked_membership`, `measurement_defined_module` |",
"| `citation_keys` | the references behind the module label, e.g. `[Fri17]; [Fad23]` |",
"| `evidence_note` | one sentence a reader can check |",
"| `deep_search_url` | where the literature check lives |",
"",
"The evidence columns may instead live in a separate one-row-per-module table",
"passed as `dmsa_sets(x, audit = \"my_audit.csv\")`. Either way DMSA surfaces",
"them whenever a module-level result is printed, because a module-level finding",
"is only as good as the module's definition.",
"",
"## Rules the validator enforces",
"",
"1. Every `module_id` belongs to exactly one `system_id`.",
"2. Every `gene` belongs to exactly one `module_id`. A gene in two modules makes",
"   its family ambiguous and there is no principled way to charge the toll.",
"3. `system_short` is unique across systems.",
"4. Module metadata (`module`, `evidence_strength`, ...) is constant within a",
"   `module_id`.",
"5. No duplicate `module_id` + `gene` + `cpg` rows.",
"",
"A CpG *may* appear under two genes - overlapping annotation is real, and the",
"bundled Alpha cascade has 228 such CpGs. The validator reports them rather than",
"rejecting them, so you can decide.",
"",
"## Using it",
"",
"```r",
"dmsa_sets_check(\"my_sets.csv\")        # validate before you trust it",
"cas <- dmsa_sets(\"my_sets.csv\")       # load; prints the short names",
"dmsa_systems(cas)                          # the catalogue",
"",
"frame <- dmsa_frame(data = dat, cascade = cas,",
"                    systems = c(\"hpa\", \"oxytocin\"),   # modules/genes/probes = full",
"                    outcome = \"symptoms\")",
"```",
"",
"Narrowing a level below the system is allowed but it changes the declared",
"family, and therefore the multiplicity toll. Declare the narrower question in",
"advance if you intend to narrow it.")

#' Validate a candidate cascade before using it
#'
#' @param x Path, \code{data.frame}, or \code{dmsa_sets}.
#' @param verbose Print the report. Default \code{TRUE}.
#' @return Invisibly, a list with \code{ok} and one entry per check.
#' @examples
#' sets <- data.frame(
#'   system_id = "1", system_short = "hpa", system = "HPA axis",
#'   module_id = "1.1", module = "Corticosteroid receptors",
#'   gene = c("NR3C1", "NR3C1", "FKBP5"),
#'   cpg = c("cg01234567", "cg07654321", "cg11223344"))
#'
#' ## a gene in two modules, or a duplicated row, would make the family
#' ## ambiguous; missing module evidence is only a caveat
#' chk <- dmsa_sets_check(sets)
#' chk$ok
#' @export
dmsa_sets_check <- function(x, verbose = TRUE) {
  ## Evidence may live in a companion audit table rather than in the cascade
  ## itself, which is how the bundled pair is shipped. Validating the raw CSV
  ## would then report the bundled set as un-annotated, so "alpha" is validated
  ## as it is actually loaded.
  mods <- NULL
  if (inherits(x, "dmsa_sets")) {
    cas <- x$cascade; mods <- x$modules
  } else if (is.character(x) && length(x) == 1L && identical(x, "alpha")) {
    s <- dmsa_sets("alpha"); cas <- s$cascade; mods <- s$modules
  } else if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("file not found: ", x, call. = FALSE)
    cas <- .cas_read(x)
  } else cas <- as.data.frame(x, stringsAsFactors = FALSE)

  res <- list(); fail <- character(); warn <- character()
  add <- function(name, ok, detail = "", level = c("fail", "warn")) {
    level <- match.arg(level)
    res[[name]] <<- list(ok = ok, detail = detail, level = level)
    if (!ok) { if (level == "warn") warn <<- c(warn, name)
               else fail <<- c(fail, name) }
  }
  miss <- setdiff(.CAS_REQ, names(cas))
  add("required columns", !length(miss),
      if (length(miss)) paste("missing:", paste(miss, collapse = ", ")) else "")
  if (length(miss)) {
    if (verbose) .cas_report(res, fail)
    return(invisible(list(ok = FALSE, checks = res)))
  }
  for (v in intersect(c(.CAS_ID, "gene", "cpg"), names(cas)))
    cas[[v]] <- as.character(cas[[v]])

  t1 <- tapply(cas$system_id, cas$module_id, function(z) length(unique(z)))
  add("each module in one system", all(t1 == 1),
      paste("offenders:", paste(utils::head(names(t1)[t1 > 1], 5), collapse = ", ")))
  t2 <- tapply(cas$module_id, cas$gene, function(z) length(unique(z)))
  add("each gene in one module", all(t2 == 1),
      paste("offenders:", paste(utils::head(names(t2)[t2 > 1], 5), collapse = ", ")))
  if ("system_short" %in% names(cas)) {
    u <- unique(cas[, c("system_id", "system_short")])
    add("system_short unique", !anyDuplicated(u$system_short),
        paste("duplicated:", paste(unique(u$system_short[duplicated(u$system_short)]),
                                   collapse = ", ")))
  } else add("system_short present", TRUE, "absent - will be derived")
  meta <- intersect(c("module", "evidence_strength", "audit_status",
                      "citation_keys"), names(cas))
  const <- vapply(meta, function(v)
    all(tapply(cas[[v]], cas$module_id, function(z) length(unique(z))) == 1), TRUE)
  add("module metadata constant within module", all(const),
      paste("varying:", paste(names(const)[!const], collapse = ", ")))
  ## The unit of measurement is the probe, not the CpG site. EPIC v2 ships
  ## replicate designs for the same site (cg09834343_BC21 and _BC22 are two
  ## columns measuring one CpG), so a repeated module+gene+cpg triple is the
  ## manifest working as intended. Keying uniqueness on cpg called 238 such rows
  ## duplicates and declared the bundled set unusable.
  key <- if ("probe_id" %in% names(cas)) "probe_id"
         else if ("cpg" %in% names(cas)) "cpg" else NULL
  if (is.null(key)) {
    dup <- anyDuplicated(cas[, c("module_id", "gene")])
    add("no duplicate module+gene rows", dup == 0,
        if (dup) paste("first at row", dup) else "")
  } else {
    dup <- anyDuplicated(cas[, c("module_id", "gene", key)])
    add(paste0("no duplicate module+gene+", key, " rows"), dup == 0,
        if (dup) paste("first at row", dup) else "")
  }
  ev <- if (!is.null(mods) && "evidence_strength" %in% names(mods))
    sum(!is.na(mods$evidence_strength))
  else if ("evidence_strength" %in% names(cas))
    sum(!is.na(unique(cas[, c("module_id", "evidence_strength")])$evidence_strength))
  else 0L
  ## A cascade with no module evidence is usable; it just prints without the
  ## evidence banner. That is a caveat on how results read, not a defect.
  add("module evidence annotated", ev > 0,
      if (ev > 0) paste(ev, "module(s) annotated")
      else "no evidence_strength - results will print without an evidence banner",
      level = "warn")

  has_cpg <- "cpg" %in% names(cas)
  shared <- if (has_cpg)
    sum(tapply(cas$gene, cas$cpg, function(z) length(unique(z))) > 1) else 0L
  rep_probes <- if (has_cpg && "probe_id" %in% names(cas))
    sum(tapply(cas$probe_id, cas$cpg, function(z) length(unique(z))) > 1)
  else 0L

  if (verbose) {
    if (has_cpg)
      cat(sprintf("cascade: %d rows | %d systems | %d modules | %d genes | %d CpGs\n",
                  nrow(cas), length(unique(cas$system_id)),
                  length(unique(cas$module_id)), length(unique(cas$gene)),
                  length(unique(cas$cpg))))
    else
      cat(sprintf("cascade: %d rows | %d systems | %d modules | %d genes (biological)\n",
                  nrow(cas), length(unique(cas$system_id)),
                  length(unique(cas$module_id)), length(unique(cas$gene))))
    .cas_report(res, fail, warn)
    if (shared)
      cat(sprintf("note: %d CpG(s) are annotated to more than one gene - allowed, "
                  , shared),
          "but their probe-level families overlap\n", sep = "")
    if (rep_probes)
      cat(sprintf(paste0("note: %d CpG(s) carry more than one probe design ",
                         "(EPIC v2 replicates) - allowed, but the replicates ",
                         "are correlated measurements of one site\n"), rep_probes))
  }
  invisible(list(ok = !length(fail), checks = res, shared_cpgs = shared,
                 replicate_cpgs = rep_probes, warnings = warn))
}

.cas_report <- function(res, fail, warn = character()) {
  for (nm in names(res)) {
    lev <- if (isTRUE(res[[nm]]$ok)) "ok"
           else if (identical(res[[nm]]$level, "warn")) "WARN" else "FAIL"
    cat(sprintf("  [%s] %-40s %s\n", lev, nm,
                if (nzchar(res[[nm]]$detail)) res[[nm]]$detail else ""))
  }
  cat(if (length(fail))
        paste0("  -> not usable: ", length(fail), " check(s) failed\n")
      else if (length(warn))
        paste0("  -> usable; ", length(warn),
               " caveat(s) to carry into how results read\n")
      else "  -> usable\n")
}

## Polarity banner. w_g is the multiplier that decides which way a system score
## points, so the evidence behind it belongs next to the selection, not in a
## separate file the reader has to think to open.
.cas_polarity_banner <- function(x, indent = "") {
  pol <- x$polarity
  if (is.null(pol)) {
    cat(indent, "polarity: none attached - a system score would weight every gene +1\n",
        sep = "")
    return(invisible(NULL))
  }
  p <- pol$polarity
  if (!nrow(p)) return(invisible(NULL))
  g <- table(factor(p$grade, levels = c("user", "curated", "database",
                                        "literature", "heuristic", "none",
                                        "unstated")))
  g <- g[g > 0]
  cat(sprintf("%spolarity: %d signed, %d off-axis (%s)\n", indent,
              sum(p$w_g != 0), sum(p$w_g == 0),
              paste(sprintf("%s %d", names(g), g), collapse = ", ")))
  nf <- sum(nzchar(stats::na.omit(p$review_flag)))
  if (nf) cat(indent, sprintf("  %d row(s) need a decision - dmsa_polarity_review()\n", nf),
              sep = "")
  invisible(NULL)
}
