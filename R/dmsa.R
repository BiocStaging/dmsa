# dmsa: Directional Methylation Set Analysis
# Two-layer sign chain: s_j = d_j (CpG -> expression, from cpgdirection)
#                             x w_g (gene -> system activation, curated polarity)
# Gene-level analysis uses d_j alone; system-level uses d_j * w_g.

# Session-scoped memo so the all-one-direction note in dmsa_align() is said
# once per distinct probe set, not once per lens battery (three lenses would
# otherwise repeat it verbatim three times per report).
.dmsa_align_once <- new.env(parent = emptyenv())

#' Alpha panel gene -> system map (bundled)
#'
#' Read the Project Alpha gene codebook shipped with the package:
#' 549 genes assigned to 30 biological systems.
#' @return data.frame with system_id, system, gene, n_kept, usability
#' @examples
#' gs <- alpha_gene_systems()
#' dim(gs)
#' length(unique(gs$system_id))
#' head(gs[, c("system_id", "system", "gene", "n_kept", "usability")], 4)
#' @export
alpha_gene_systems <- function() {
  p <- system.file("extdata", "alpha_gene_systems.csv", package = "dmsa")
  utils::read.csv(p, stringsAsFactors = FALSE)
}

#' Alpha panel gene -> system-activation polarity (bundled)
#'
#' Curated w_g per gene: +1 gene product increases the system's activation
#' tone, -1 decreases it (e.g. NR3C1 in the HPA axis: more GR = stronger
#' negative feedback = LOWER axis tone), 0 = not on the activation axis
#' (readouts, specificity controls, ambiguous).
#'
#' Returns the audited 2026c table: 1,234 genes across all 30 systems. A
#' legacy 115-gene draft covering six systems also ships, and is used only as
#' a fallback if the 2026c file is missing from the installation.
#'
#' The 2026c table is why a system-level alignment works at all: an
#' \code{alpha_polarity()} that covered only six systems made
#' \code{dmsa_align(level = "system")} fail on most of the panel, since a
#' gene with no w_g entry stops the alignment by design.
#' @return data.frame with system_id, system, gene, w_g, role, confidence and,
#'   for the 2026c table, module, w_g_source, anchor, evidence and citation.
#'   The \code{"status"} attribute is \code{"audited"} or \code{"draft"}.
#' @examples
#' pol <- alpha_polarity()
#' nrow(pol)
#' head(pol[, c("system", "gene", "w_g", "role")], 4)
#' # CRH drives the HPA axis (w_g = +1); NR3C1 is its brake, so more GR means
#' # stronger negative feedback and LOWER axis tone (w_g = -1)
#' pol[pol$gene %in% c("CRH", "NR3C1") & pol$system_id == 2,
#'     c("gene", "w_g", "role")]
#' @export
alpha_polarity <- function() {
  p <- system.file("extdata", "alpha_polarity_2026c.csv", package = "dmsa")
  audited <- nzchar(p)
  if (!audited)
    p <- system.file("extdata", "alpha_polarity.csv", package = "dmsa")
  if (!nzchar(p))
    stop("no bundled Alpha polarity table found; reinstall dmsa", call. = FALSE)
  out <- utils::read.csv(p, stringsAsFactors = FALSE)
  names(out)[names(out) == "w_g_draft"] <- "w_g"
  out$system_id <- as.character(out$system_id)
  attr(out, "status") <- if (audited) "audited" else "draft"
  out
}

#' Build the aligned sign s_j for a set of probes
#'
#' @param direction Either a plain character vector of CpG identifiers, in
#'   which case the direction calls and the gene mapping are taken from the
#'   map bundled with dmsa (see \code{\link{dmsa_directions}}), or a table of
#'   calls you supply yourself. A supplied table may be the output of
#'   \code{cpgdirection::cpg_expression_direction()} or
#'   \code{cpgdirection::cpg_gene_pairs()} (needs \code{cpg_id} or
#'   \code{input}, \code{best_direction}, and if available
#'   \code{best_accuracy}), or any data.frame with \code{cpg},
#'   \code{d} (+1/-1) and optionally \code{p_plus} = P(d = +1). When the
#'   table carries \code{best_accuracy} (cpgdirection's November 2026 builds: the
#'   validated probability that \code{best_direction} is right, in the
#'   call's tissue), P(d = +1) is \code{best_accuracy} for a +1 call and
#'   \code{1 - best_accuracy} for a -1 call, ahead of every other column; a
#'   direction with no \code{best_accuracy} value is then unusable. Without
#'   that column, P(d = +1) comes from \code{p_plus}, then
#'   \code{probability_plus1}, then \code{best_confidence} /
#'   \code{confidence}, then the SMR tier.
#' @param genes The set-membership gene of each probe - which gene this probe
#'   is being read against. Required, and the same length as the probes, when
#'   \code{direction} is a table. Optional when \code{direction} is a
#'   character vector of probe IDs: leave it \code{NULL} to take every gene
#'   the bundled map maps each probe to, or supply one gene per probe to pin
#'   specific pairs.
#' @param level \code{"gene"} (alignment = d_j only) or \code{"system"}
#'   (alignment = d_j * w_g).
#' @param polarity For \code{level = "system"}: where w_g comes from.
#'   \code{"alpha"} (default) uses the bundled Alpha polarity table.
#'   Or supply your own data.frame with columns \code{gene} and \code{w_g}
#'   (+1 / -1 / 0). User-supplied values OVERRIDE the bundled table where
#'   both cover a gene.
#' @param system_id Optional integer: restrict the bundled polarity lookup
#'   to one Alpha system (a gene can carry different roles in different
#'   systems).
#' @param tissue Which bundled direction layer to read when \code{direction}
#'   is a character vector of probe IDs: \code{"blood"} (default) or
#'   \code{"epithelium"}. Choose by what the sample is made of - adult saliva
#'   is about 82% immune and takes blood calls; neonatal and infant buccal
#'   swabs are .89 to 1.00 epithelial and take epithelial calls. Getting this
#'   wrong fails quietly: coverage looks full while every direction comes from
#'   the wrong tissue. Ignored when \code{direction} is a table you supply.
#'   Further tissues live in \code{cpgdirection}.
#' @param missing_polarity What to do with a called probe whose gene has no
#'   polarity entry: \code{"error"} (default - stop and list the genes so
#'   the user can specify each, per protocol), \code{"zero"} (treat as
#'   off-axis, probe drops from the system pool but is reported), or
#'   \code{"drop"} (silently exclude).
#' @return data.frame: probe, gene, d, p_plus, w_g, s (= d * w_g at system
#'   level, d at gene level), p_s_plus = P(s = +1) (chained), usable
#'   (logical), and reason for every non-usable probe.
#' @examples
#' # The short way: probe IDs alone. Directions and gene mapping come from
#' # the bundled map, so nothing has to be obtained first.
#' al <- dmsa_align(c("cg00052046", "cg00176879", "cg00308631"))
#' al[, c("probe", "gene", "d", "s", "usable")]
#'
#' # Which map produced this alignment travels with it.
#' attr(al, "direction_map")
#'
#' # Pin each probe to one gene when the panel says which gene it belongs to.
#' dmsa_align(c("cg00052046", "cg00176879"), genes = c("AVP", "AVP"))$s
#'
#' # Direction is tissue-specific, so which layer to read is part of the
#' # analysis, not a default to leave alone. An adult saliva sample is about
#' # 82% immune and takes blood calls; a neonatal buccal swab is near-pure
#' # epithelium and takes epithelial ones.
#' bu <- dmsa_align("cg00176879", tissue = "epithelium")
#' attr(bu, "direction_map")
#'
#' # d is the CpG -> expression sign, w_g the gene -> system-tone sign.
#' # A promoter CpG silences NR3C1 (d = -1) and GR is the HPA brake (w_g = -1),
#' # so the probe votes for HIGHER axis tone: s = d * w_g = +1.
#' g <- c("CRH", "NR3C1")
#' dcall <- data.frame(cpg = c("cg_crh", "cg_nr3c1"), d = c(1, -1),
#'                     p_plus = c(0.9, 0.1))
#' pol <- data.frame(gene = g, w_g = c(1, -1))
#' dmsa_align(dcall, genes = g, level = "gene")
#' dmsa_align(dcall, genes = g, level = "system", polarity = pol)
#' @export
dmsa_align <- function(direction, genes = NULL,
                       level = c("gene", "system"),
                       polarity = "alpha",
                       system_id = NULL,
                       missing_polarity = c("error", "zero", "drop"),
                       tissue = "blood") {
  level <- match.arg(level)
  missing_polarity <- match.arg(missing_polarity)

  ## -- probe IDs alone: take the calls from the bundled map ---------------
  ## The commonest thing a user has after normalising is a set of probe
  ## names and nothing else. Requiring them to go and find a direction
  ## resource first is the barrier this branch removes.
  map_meta <- NULL
  if (is.character(direction) && is.null(dim(direction))) {
    asked <- unique(as.character(direction))
    got <- dmsa_directions(probes = as.character(direction),
                           genes = genes, tissue = tissue)
    map_meta <- attr(got, "dmsa_direction_map")
    covered <- length(intersect(asked, got$probe))
    if (!nrow(got))
      stop("none of the ", length(asked), " probe(s) carry a direction call ",
           "in the bundled ", tissue, " map. Supply your own direction table, ",
           "or see cpgdirection for other tissues and the full resource.",
           call. = FALSE)
    message(sprintf(
      "dmsa: %d of %d probes carry a call in the bundled %s map (%.0f%%), %d probe-gene pairs [map %s]",
      covered, length(asked), tissue, 100 * covered / length(asked),
      nrow(got), map_meta$version))
    ## E9: rows come back in the order of the probes the CALLER passed, not
    ## the bundled map's storage order - every consumer joins by position
    got <- got[order(match(got$probe, asked)), , drop = FALSE]
    genes <- got$gene
    direction <- got[, c("probe", "d", "p_plus")]
    names(direction)[1] <- "cpg_id"
  }

  ## -- normalise the direction input ------------------------------------
  dd <- as.data.frame(direction)
  cpg_col <- intersect(c("cpg_id", "cpg", "input"), names(dd))[1]
  if (is.na(cpg_col)) stop("direction needs a 'cpg_id' / 'cpg' / 'input' column")
  d_col <- intersect(c("best_direction", "d", "direction"), names(dd))[1]
  if (is.na(d_col)) stop("direction needs 'best_direction' or 'd'")
  out <- data.frame(probe = as.character(dd[[cpg_col]]),
                    d     = suppressWarnings(as.numeric(dd[[d_col]])),
                    stringsAsFactors = FALSE)
  ## ---- one-way instruments -----------------------------------------------
  ## cpgdirection's distance curves require unanimity across three tissues,
  ## and the blood curve never exceeds 0.449, so `distance_only` cannot
  ## return +1 for any CpG at any distance - every direction it contributes
  ## is -1 by construction (see cpg_expression_direction's documentation).
  ## Aligning on a block of probes that all point the same way by
  ## construction manufactures coherence that is not in the data.
  if ("best_evidence" %in% names(dd)) {
    oneway <- c("distance_only", "targeted_last_resort")
    hit <- as.character(dd$best_evidence) %in% oneway
    if (any(hit))
      warning(sum(hit), " probe(s) carry a ONE-WAY direction call (",
              paste(intersect(oneway, unique(as.character(dd$best_evidence))),
                    collapse = ", "),
              "). Those layers cannot return +1 by construction, so every ",
              "call they contribute is -1. Filter them out before aligning ",
              "unless you have a specific reason not to.", call. = FALSE)
  }
  ## A large block of probes sharing one direction is either real biology or
  ## an instrument artefact - but WHICH depends on where the calls came from
  ## (PI, 2026-08-29: "if cpgdirection works, it works; if the user CAN do
  ## something, tell them what and where"). Calls resolved by cpgdirection
  ## carry evidence tiers, and an all-one-way set there is usually genuine
  ## curation: say so once, as a message, with the exact places to verify.
  ## Calls from a user map or a bare `d` vector carry no evidence, and an
  ## all-one-way block is the signature of a mis-built map: keep the
  ## warning, and make it name the check. Either way it is said ONCE per
  ## set per session, not once per lens battery.
  if (sum(!is.na(out$d)) >= 50) {
    fp <- mean(out$d[!is.na(out$d)] > 0)
    if (fp == 0 || fp == 1) {
      .dirtxt <- if (fp == 1) "+1" else "-1"
      .has_ev <- any(c("best_evidence", "direction_tier") %in% names(dd))
      .sig <- paste(sum(!is.na(out$d)), .dirtxt, .has_ev,
                    if (!is.null(out$probe)) out$probe[1] else "",
                    if (!is.null(out$probe)) out$probe[nrow(out)] else "")
      if (!exists(.sig, envir = .dmsa_align_once, inherits = FALSE)) {
        assign(.sig, TRUE, envir = .dmsa_align_once)
        if (.has_ev)
          message("dmsa_align(): all ", sum(!is.na(out$d)),
                  " direction calls in this set are ", .dirtxt,
                  ". With cpgdirection-resolved calls this is usually ",
                  "genuine curation (the whole set's methylation ",
                  if (fp == 1) "raises" else "lowers",
                  " expression). To verify, read the per-probe calls with ",
                  "their evidence: tables/analysis_set.csv (columns ",
                  "direction, direction_tier, evidence), or frame$map. ",
                  "One-way evidence layers are flagged separately. ",
                  "(Shown once per set.)")
        else
          warning("every one of ", sum(!is.na(out$d)),
                  " direction calls is ", .dirtxt,
                  ", and these calls carry NO evidence tiers (a user map ",
                  "or bare `d` vector). A biological set usually mixes ",
                  "+1 and -1, so check the map you supplied - ",
                  "table(map$best_direction), or frame$map after the ",
                  "build - before trusting an aligned test.",
                  call. = FALSE)
      }
    }
  }

  ## P(d = +1). A probe whose confidence cannot be recovered here does not
  ## fail loudly - it acquires an NA weight and evaporates several steps
  ## later inside dmsa_triangulate(), which is the worst possible place for
  ## it to happen. So take every source cpgdirection offers, in order of how
  ## direct each one is, and say plainly what is left over.
  .conf_to_pplus <- function(conf, d) {
    conf <- pmin(pmax(as.numeric(conf), 0), 1)
    ifelse(is.na(d), NA_real_, ifelse(d > 0, 0.5 + conf / 2, 0.5 - conf / 2))
  }
  ## E17 (PI ruling, 2026-11-07). cpgdirection (November 2026) gives every direction
  ## its validated accuracy in the call's tissue - best_accuracy, the
  ## probability that best_direction is right (the lower 95% bound by
  ## default; the point estimate under options(cpgdirection.accuracy =
  ## "median")). The confidence columns below are not accuracies: measured
  ## calls read 1 whatever their validated accuracy, SMR rows carry none,
  ## tumour-consensus rows carry none, distance reads a curve distance from
  ## 0.5. So when the table carries best_accuracy it decides P(d = +1) for
  ## every row, ahead of every other column, and the ladder below serves only
  ## tables without it (older cpgdirection, user maps).
  use_acc <- "best_accuracy" %in% names(dd)
  if (use_acc) {
    a <- pmin(pmax(suppressWarnings(as.numeric(dd$best_accuracy)), 0), 1)
    out$p_plus <- ifelse(is.na(out$d), NA_real_, ifelse(out$d > 0, a, 1 - a))
  } else if ("p_plus" %in% names(dd)) {
    out$p_plus <- suppressWarnings(as.numeric(dd$p_plus))
  } else {
    out$p_plus <- NA_real_
  }
  ## the single-tissue schema
  if (!use_acc && anyNA(out$p_plus) && "probability_plus1" %in% names(dd)) {
    k <- is.na(out$p_plus)
    out$p_plus[k] <- suppressWarnings(as.numeric(dd$probability_plus1))[k]
  }
  ## confidence in the CALLED direction (single tissue, or the winning layer
  ## of a tissue = "all" query)
  for (cn in c("best_confidence", "confidence")) {
    if (!use_acc && anyNA(out$p_plus) && cn %in% names(dd)) {
      k <- is.na(out$p_plus)
      out$p_plus[k] <- .conf_to_pplus(dd[[cn]], out$d)[k]
    }
  }
  ## cpgdirection's tissue = "all" schema fills best_confidence only for its
  ## catalogue layers, so SMR rows arrive with none. Their confidence is
  ## published per tier, with validation counts, and that is what to use -
  ## NOT the universal/distance probability, which would substitute a
  ## 0.60-0.65 prior for a 0.95-0.97 validated causal call.
  if (!use_acc && anyNA(out$p_plus) && "smr_tier" %in% names(dd)) {
    acc <- c(S1 = 0.96, S2 = 0.85, S3 = 0.705)   # midpoints of the published ranges
    a <- acc[as.character(dd$smr_tier)]
    k <- is.na(out$p_plus) & !is.na(a)
    if (any(k)) out$p_plus[k] <- .conf_to_pplus(2 * a - 1, out$d)[k]
  }
  ## Absence of a confidence COLUMN and absence of a confidence VALUE are
  ## different things. A plain d-only table has no column, and treating those
  ## calls as certain is the documented behaviour. A table that carries
  ## confidence columns but leaves this row empty is a gap, and filling it
  ## with certainty would hand maximum weight to the least-supported rows -
  ## which is exactly backwards. Those are marked unusable instead.
  conf_cols <- c("best_accuracy", "p_plus", "probability_plus1", "best_confidence",
                 "confidence", "smr_tier")
  has_conf_col <- any(conf_cols %in% names(dd))
  if (anyNA(out$p_plus)) {
    k <- is.na(out$p_plus) & !is.na(out$d)
    if (any(k) && !has_conf_col) {
      out$p_plus[k] <- ifelse(out$d[k] > 0, 1, 0)
    } else if (any(k)) {
      warning(sum(k), " probe(s) carry a direction but no confidence value, ",
              "in a table that does have confidence columns. They are marked ",
              "unusable rather than assumed certain - assuming certainty ",
              "would give the least-supported rows the most weight.",
              call. = FALSE)
    }
  }
  if (is.null(genes))
    stop("`genes` is required when `direction` is a table. Pass a character ",
         "vector of probe IDs instead and dmsa will take both the genes and ",
         "the direction calls from the bundled map.", call. = FALSE)
  if (length(genes) != nrow(out))
    stop("genes must have one entry per probe row (got ", length(genes),
         " for ", nrow(out), " probes)")
  out$gene <- as.character(genes)

  ## -- polarity layer -----------------------------------------------------
  if (level == "system") {
    pol <- NULL
    if (is.character(polarity) && identical(polarity, "alpha")) {
      pol <- alpha_polarity()
      pol_status <- attr(pol, "status")     # subsetting drops attributes
      if (!is.null(system_id))
        pol <- pol[pol$system_id == as.character(system_id), , drop = FALSE]
      in_sys <- if (is.null(system_id)) "" else
        paste0(" in system ", system_id)
      message(sprintf(
        "dmsa: using bundled Alpha polarity table (%s, %d genes%s).",
        pol_status, nrow(pol), in_sys))
      if (!nrow(pol))
        stop("no polarity rows for system_id = ", system_id,
             "; check the id against dmsa_systems()", call. = FALSE)
    } else {
      ## E4 fix (2026-08-29, PI-approved): a user table MERGES with the
      ## bundled curation, exactly as the documentation has always promised -
      ## user rows override gene-for-gene, the bundled table fills every gene
      ## the user did not name. The old code silently REPLACED the bundled
      ## table, so polarity = data.frame(gene = "NR3C1", w_g = +1) built a
      ## system in which every other gene had no polarity at all.
      usr <- as.data.frame(polarity)
      if (!all(c("gene", "w_g") %in% names(usr)))
        stop("user polarity needs columns 'gene' and 'w_g' (+1/-1/0)")
      usr$gene <- as.character(usr$gene)
      if (anyDuplicated(usr$gene[if (!is.null(usr$system_id) &&
                                     !is.null(system_id))
        as.character(usr$system_id) == as.character(system_id) else
          rep(TRUE, nrow(usr))]))
        stop("user polarity lists the same gene more than once",
             call. = FALSE)
      bnd <- tryCatch(alpha_polarity(), error = function(e) NULL)
      if (!is.null(bnd) && !is.null(system_id))
        bnd <- bnd[bnd$system_id == as.character(system_id), , drop = FALSE]
      if (!is.null(usr$system_id) && !is.null(system_id))
        usr <- usr[as.character(usr$system_id) == as.character(system_id), ,
                   drop = FALSE]
      if (!is.null(bnd) && nrow(bnd)) {
        fill <- bnd[!bnd$gene %in% usr$gene, , drop = FALSE]
        pol <- rbind(usr[intersect(c("gene", "w_g"), names(usr))],
                     fill[c("gene", "w_g")])
        n_here <- length(intersect(out$gene, usr$gene))
        message(sprintf(paste0("dmsa: polarity merged - %d gene(s) from ",
                               "your table (%d used here), %d filled from ",
                               "the bundled Alpha curation."),
                        nrow(usr), n_here, nrow(fill)))
      } else pol <- usr
    }
    idx <- match(out$gene, pol$gene)
    out$w_g <- pol$w_g[idx]

    unknown <- sort(unique(out$gene[is.na(out$w_g) & !is.na(out$d)]))
    if (length(unknown) > 0) {
      if (missing_polarity == "error") {
        stop("No polarity (w_g) for gene(s): ", paste(unknown, collapse = ", "),
             "\nSpecify each via the 'polarity' argument, e.g.\n",
             "  polarity = data.frame(gene = c(",
             paste0('"', unknown[seq_len(min(2, length(unknown)))], '"',
                    collapse = ", "),
             ", ...), w_g = c(...))\n",
             "(+1 activates the system, -1 opposes it - e.g. NR3C1 in HPA, ",
             "0 = off-axis)", call. = FALSE)
      } else if (missing_polarity == "zero") {
        warning("w_g set to 0 (off-axis) for: ", paste(unknown, collapse = ", "))
        out$w_g[is.na(out$w_g)] <- 0
      } else {
        out <- out[!is.na(out$w_g) | is.na(out$d), , drop = FALSE]
      }
    }
    out$s <- out$d * out$w_g
    ## chained probability: P(s=+1) = p_d*p_w + (1-p_d)(1-p_w); w_g known -> p_w in {0,1}
    p_w <- ifelse(out$w_g > 0, 1, ifelse(out$w_g < 0, 0, NA_real_))
    out$p_s_plus <- out$p_plus * p_w + (1 - out$p_plus) * (1 - p_w)
  } else {
    out$w_g <- NA_real_
    out$s <- out$d
    out$p_s_plus <- out$p_plus
  }

  out$p_s_plus[is.na(out$d)] <- NA_real_   # an abstained probe carries no sign certainty
  ## A probe with no confidence is not usable: its weight would be NA, and an
  ## NA weight does not fail here - it disappears silently inside
  ## dmsa_triangulate(). Refuse it now, with a reason, so the loss is visible.
  out$usable <- !is.na(out$s) & out$s != 0 & !is.na(out$p_plus)
  out$reason <- ifelse(is.na(out$d), "no_direction_call",
                ifelse(is.na(out$p_plus), "no_confidence",
                ifelse(!is.na(out$w_g) & out$w_g == 0, "off_axis_gene",
                ifelse(is.na(out$s), "no_polarity", ""))))
  ## record which direction map produced this alignment. A finding that moved
  ## because the annotation moved is the failure this package exists to name,
  ## so the provenance travels with the result rather than being remembered.
  attr(out, "direction_map") <- if (!is.null(map_meta))
    sprintf("dmsa bundled %s map v%s (from %s, %s)", map_meta$tissue,
            map_meta$version, map_meta$source, map_meta$source_doi)
  else "user-supplied direction table"
  class(out) <- c("dmsa_alignment", class(out))
  out
}

#' Fixed-sign or expected-sign DMSA pooled test
#'
#' Pools probe-level coefficients after alignment. Two estimators:
#' \code{"fixed"} multiplies each b_j by the hard sign s_j; expected-sign
#' (\code{"expected"}) multiplies by E[s_j] = 2 P(s_j = +1) - 1, the
#' first-order latent-sign approximation - low-tier calls contribute in
#' proportion to their calibrated certainty, which is why no tier filtering
#' is needed upstream.
#'
#' The z test here is descriptive; the reportable error control is the
#' family-wise permutation (\code{dmsa_perm_pvalue}) on the same statistic.
#'
#' @param b,se Probe-level coefficient and SE from YOUR model (fitted on
#'   M-values standardised within probe; same estimator observed and null).
#' @param alignment A \code{dmsa_align()} result covering the same probes,
#'   or vectors via \code{s} / \code{p_s_plus}.
#' @param method \code{"fixed"} or \code{"expected"}.
#' @param w Optional per-probe reliability weight (see
#'   \code{\link{dmsa_relweights}}); scales each aligned multiplier. \code{NULL}
#'   (default) or all ones is the flat engine.
#' @return list: estimate (pooled aligned slope), se, z, p_normal,
#'   n_used, n_excluded, table (per-probe contributions).
#' @examples
#' ## Two promoter CpGs in NR3C1: methylation silences the gene (d = -1), and
#' ## GR is the HPA brake (w_g = -1), so both probes vote for HIGHER axis tone.
#' calls <- data.frame(cpg = c("cg01", "cg02"), d = c(-1, -1), p_plus = c(.1, .2))
#' al <- dmsa_align(calls, genes = c("NR3C1", "NR3C1"), level = "system",
#'                  polarity = data.frame(gene = "NR3C1", w_g = -1))
#'
#' tt <- dmsa_test(b = c(0.28, 0.19), se = c(0.09, 0.08), al, method = "expected")
#' tt[c("estimate", "se", "z", "p_normal", "n_used")]
#' tt$table
#' @export
dmsa_test <- function(b, se, alignment,
                      method = c("expected", "fixed"), w = NULL) {
  method <- match.arg(method)
  al <- as.data.frame(alignment)
  if (nrow(al) != length(b) || length(b) != length(se))
    stop("b, se and alignment must cover the same probes, in order")
  mult <- if (method == "fixed") al$s else (2 * al$p_s_plus - 1)
  if (!is.null(w)) mult <- mult * w
  use <- !is.na(mult) & mult != 0 & !is.na(b) & !is.na(se) & se > 0
  ## GLS estimator of the common aligned effect mu under b_j ~ N(mu * m_j, se_j^2):
  ## mu_hat = sum(m_j b_j / se_j^2) / sum(m_j^2 / se_j^2)
  ## (reduces to standard inverse-variance pooling when |m_j| = 1, i.e. fixed-sign;
  ##  under expected-sign, low-certainty probes contribute ~nothing rather than dominating)
  num <- sum(mult[use] * b[use] / se[use]^2)
  den <- sum(mult[use]^2 / se[use]^2)
  est <- num / den
  pse <- sqrt(1 / den)
  z   <- est / pse
  structure(list(
    estimate = est, se = pse, z = z,
    p_normal = 2 * stats::pnorm(-abs(z)),
    method = method, n_used = sum(use), n_excluded = sum(!use),
    table = data.frame(probe = al$probe, gene = al$gene, b = b, se = se,
                       s = al$s, p_s_plus = al$p_s_plus,
                       multiplier = mult, b_aligned = b * mult,
                       used = use)
  ), class = "dmsa_test")
}

#' @export
print.dmsa_test <- function(x, ...) {
  cat("DMSA pooled aligned test (", x$method, "-sign)\n", sep = "")
  cat(sprintf("  estimate %.4f  se %.4f  z %.2f  p(normal) %.4g\n",
              x$estimate, x$se, x$z, x$p_normal))
  cat(sprintf("  probes used %d, excluded %d\n", x$n_used, x$n_excluded))
  cat("  NOTE: report the permutation p (dmsa_perm_pvalue), not p(normal).\n")
  invisible(x)
}

#' Permutation p-value for a DMSA statistic
#'
#' @param observed The observed statistic (e.g. z from \code{dmsa_test}).
#' @param null_stats Vector of the same statistic computed on B
#'   family-wise-permuted datasets (shuffle the exposure BY FAMILY, refit
#'   the probe models with the SAME estimator, re-run \code{dmsa_test}).
#' @return two-sided permutation p with the +1 correction (resolution
#'   floor 1/(B+1)).
#' @examples
#' set.seed(1)
#' null_z <- rnorm(999)          # z from B family-wise permuted datasets
#' dmsa_perm_pvalue(3.1, null_z)
#' ## the +1 correction puts a floor of 1/(B+1) on any permutation p
#' dmsa_perm_pvalue(50, null_z)
#' @export
dmsa_perm_pvalue <- function(observed, null_stats) {
  B <- length(null_stats)
  (1 + sum(abs(null_stats) >= abs(observed))) / (B + 1)
}

#' Per-system coverage and polarity balance
#'
#' The honest-reporting companion: how many probes each system can actually
#' align, on which side of the activation axis. A system whose brake genes
#' have no usable probes silently degrades to an activation-gene test -
#' this table is how that is caught and reported.
#' @param alignment result of dmsa_align(level = "system")
#' @return one-row data.frame per call
#' @examples
#' dcall <- data.frame(cpg = paste0("p", 1:4), d = c(1, -1, 1, NA),
#'                     p_plus = c(0.9, 0.1, 0.8, NA))
#' genes <- c("CRH", "POMC", "AVPR2", "NR3C1")
#' pol <- data.frame(gene = genes, w_g = c(1, 1, 0, -1))
#' al <- dmsa_align(dcall, genes, level = "system", polarity = pol)
#' # the only brake gene abstained, so this "system" test is really an
#' # activation-gene test - which is what the balance table makes visible
#' dmsa_balance(al)
#' @export
dmsa_balance <- function(alignment) {
  al <- as.data.frame(alignment)
  data.frame(
    probes            = nrow(al),
    called            = sum(!is.na(al$d)),
    usable            = sum(al$usable),
    activation_probes = sum(al$usable & al$w_g > 0, na.rm = TRUE),
    brake_probes      = sum(al$usable & al$w_g < 0, na.rm = TRUE),
    off_axis          = sum(al$reason == "off_axis_gene", na.rm = TRUE),
    no_direction      = sum(al$reason == "no_direction_call", na.rm = TRUE),
    activation_genes  = length(unique(al$gene[al$usable & al$w_g > 0])),
    brake_genes       = length(unique(al$gene[al$usable & al$w_g < 0]))
  )
}
