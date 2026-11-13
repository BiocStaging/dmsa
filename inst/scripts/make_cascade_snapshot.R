## Builds inst/extdata/coverage_v6_cpgdirection.csv.gz: cpgdirection's pair
## directions for the Alpha selection cascade in blood and nasal epithelium.
## It is the map of the path functions (dmsa_mediate(), dmsa_modmed(),
## dmsa_transmit(), dmsa_concordance()) when cpgdirection is not installed.
##
## Recipe (PI ruling 2026-11-12): the cascade's CpGs go through
## cpgdirection::cpg_gene_pairs() with the Alpha reference genes as the filter,
## exactly as dmsa_frame(direction_source = "cpgdirection") calls it, and
## through dmsa's own pair runtime (.frame_pairs_runtime), so every row is
## what the live path returns for that pair: best_direction, the source that
## decided it (best_evidence, best_rung), its validated accuracy in the
## tissue (best_accuracy) and p_plus as dmsa_align() resolves it (E17:
## best_accuracy for a +1 call, 1 - best_accuracy for a -1 call). Only the
## cascade's own CpG x gene assignments are kept, as in v5. Rows are ordered
## by tissue, system, module, gene, probe (radix order), which also keeps the
## compressed file small.
##
## v5 was built from cpgdirection 2.6.0. v6 needs cpgdirection 1.1.0 or
## later (the SMR accuracy bands and the distance bands) with
## cpgdirectionData 0.99.8 or later.
##
## Usage, from the dmsa package root:
##   Rscript inst/scripts/make_cascade_snapshot.R [--cpgdirection <source dir>]
##     [--data-dir <dir>] [--gene-tss <csv>] [--compare <old snapshot>]
##     [--out <file>] [--report <file.md>]
## --cpgdirection  a cpgdirection source checkout to load instead of the
##                 installed package;
## --data-dir      a local cpgdirection data directory
##                 (options(cpgdirection.data_dir = ));
## --gene-tss      a csv with gene, chr, tss for the distance layer, for a
##                 machine without the TxDb annotation packages (the released
##                 snapshot is built with them);
## --compare       an earlier snapshot to compare with (default: v5, if it
##                 is still in inst/extdata).

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
cpgd_src <- arg("--cpgdirection")
data_dir <- arg("--data-dir")
tss_file <- arg("--gene-tss")
out_file <- arg("--out", file.path("inst", "extdata", "coverage_v6_cpgdirection.csv.gz"))
cmp_file <- arg("--compare", file.path("inst", "extdata", "coverage_v5_cpgdirection.csv.gz"))
rep_file <- arg("--report")

if (!is.null(data_dir)) options(cpgdirection.data_dir = path.expand(data_dir))
if (!is.null(cpgd_src)) {
  suppressMessages(devtools::load_all(path.expand(cpgd_src), quiet = TRUE))
} else {
  suppressMessages(library(cpgdirection))
}
suppressMessages(devtools::load_all(".", quiet = TRUE))

## the build must rest on the 1.1.0 layers
smr <- cpgdirection::cpgd_smr_directions()
stopifnot("this cpgdirection has no per-pair SMR accuracy (needs 1.1.0 or later)" =
            "smr_accuracy" %in% names(smr))
cpgd_commit <- if (!is.null(cpgd_src)) {
  tryCatch(system2("git", c("-C", shQuote(path.expand(cpgd_src)), "rev-parse",
                            "--short=7", "HEAD"), stdout = TRUE, stderr = FALSE),
           error = function(e) NA_character_)
} else NA_character_

cas <- utils::read.csv(gzfile(system.file("extdata", "alpha_cascade_2026c.csv.gz",
                                          package = "dmsa")),
                       stringsAsFactors = FALSE)
ref <- alpha_reference()
cpgs <- unique(cas$cpg)
src <- NULL
if (!is.null(tss_file)) {
  gt <- utils::read.csv(tss_file, stringsAsFactors = FALSE)
  stopifnot(all(c("gene", "chr", "tss") %in% names(gt)))
  src <- list(gene_tss = data.table::as.data.table(gt[, c("gene", "chr", "tss")]))
}

keep <- c("probe", "column", "gene", "system_id", "system", "module_id", "module",
          "best_direction", "p_plus", "best_evidence", "direction_tier",
          "best_rung", "best_accuracy")
one <- function(tissue) {
  cl <- list(cpgs = cpgs, genes = unique(ref$systems$gene), gene_mode = "filter",
             annotation_mode = "augment", tissue = tissue, include_brain = FALSE,
             probe_qc = "exclude", verbose = FALSE)
  if (!is.null(src)) cl$sources <- src
  p <- suppressWarnings(suppressMessages(do.call(cpgdirection::cpg_gene_pairs, cl)))
  rt <- .frame_pairs_runtime(cpgs, ref, tissue = tissue, pairs = as.data.frame(p))
  f <- rt$full
  f <- f[paste(f$probe, f$gene) %in% paste(cas$cpg, cas$gene), keep, drop = FALSE]
  f$tissue <- tissue
  f[order(as.numeric(f$system_id), f$module_id, f$gene, f$probe, method = "radix"), ,
    drop = FALSE]
}
snap <- rbind(one("blood"), one("nasal_epithelium"))
rownames(snap) <- NULL
snap$best_direction <- as.integer(snap$best_direction)
snap$p_plus <- round(snap$p_plus, 4)
snap$best_accuracy <- round(snap$best_accuracy, 4)
stopifnot(all(snap$best_direction %in% c(-1L, 1L)), all(is.finite(snap$p_plus)),
          all(is.finite(snap$best_accuracy)),
          !anyDuplicated(snap[, c("tissue", "probe", "gene")]))

con <- gzfile(out_file, "wb", compression = 9)
utils::write.csv(snap, con, row.names = FALSE)
close(con)

## ---- report ---------------------------------------------------------------
lines <- character(0)
say <- function(...) lines <<- c(lines, paste0(...))
say("# Alpha cascade snapshot v6 (", format(Sys.Date()), ")")
say("")
say("- cpgdirection ", as.character(utils::packageVersion("cpgdirection")),
    if (!is.na(cpgd_commit)) paste0(" (source checkout at ", cpgd_commit, ")") else "",
    "; cpgdirectionData ",
    tryCatch(as.character(utils::packageVersion("cpgdirectionData")),
             error = function(e) "not installed"),
    "; SMR table ", format(nrow(smr), big.mark = ","), " pairs")
say("- distance layer: ", if (is.null(src)) "the package's own gene TSS table (TxDb)"
    else paste0("a supplied gene table (", basename(tss_file), "); NOT the release build"))
say("- cascade: ", format(length(cpgs), big.mark = ","), " CpGs, ",
    format(nrow(unique(cas[, c("cpg", "gene")])), big.mark = ","), " CpG x gene pairs")
say("- file: `", basename(out_file), "`, ", format(nrow(snap), big.mark = ","), " rows, ",
    format(file.size(out_file), big.mark = ","), " bytes, md5 ",
    unname(tools::md5sum(out_file)))
old <- if (file.exists(cmp_file)) utils::read.csv(gzfile(cmp_file), stringsAsFactors = FALSE) else NULL
for (tis in c("blood", "nasal_epithelium")) {
  n <- snap[snap$tissue == tis, ]
  say("")
  say("## ", tis)
  say("")
  say("- ", format(nrow(n), big.mark = ","), " usable cascade pairs, ",
      format(length(unique(n$probe)), big.mark = ","), " probes, ",
      length(unique(n$gene)), " genes; median weight |2p - 1| ",
      sprintf("%.2f", stats::median(abs(2 * n$p_plus - 1))))
  say("")
  say("| deciding source | pairs | median accuracy |")
  say("|---|---|---|")
  ev <- sort(table(n$best_evidence), decreasing = TRUE)
  for (e in names(ev))
    say("| ", e, " | ", format(ev[[e]], big.mark = ","), " | ",
        sprintf("%.3f", stats::median(n$best_accuracy[n$best_evidence == e])), " |")
  if (!is.null(old)) {
    o <- old[old$tissue == tis, ]
    ko <- paste(o$probe, o$gene); kn <- paste(n$probe, n$gene)
    both <- intersect(ko, kn)
    d_o <- o$best_direction[match(both, ko)]; d_n <- n$best_direction[match(both, kn)]
    say("")
    say("- against `", basename(cmp_file), "` (", format(nrow(o), big.mark = ","),
        " pairs): ", format(length(both), big.mark = ","), " in both, ",
        format(sum(d_o != d_n), big.mark = ","), " with a different sign (",
        sprintf("%.1f%%", 100 * mean(d_o != d_n)), "); ",
        format(sum(!ko %in% kn), big.mark = ","), " dropped; ",
        format(sum(!kn %in% ko), big.mark = ","), " added")
  }
}
cat(lines, sep = "\n")
if (!is.null(rep_file)) writeLines(lines, rep_file)
