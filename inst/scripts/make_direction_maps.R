## Builds the two direction maps bundled in inst/extdata
## (direction_blood_hg19.rds, direction_epithelium_hg19.rds) from the
## cpgdirection catalogue lookups served by cpgdirectionData.
##
## Recipe: every lookup row that carries a call (status PREDICTED or
## DIRECT_eQTM, i.e. a non-missing direction); ABSTAIN rows are dropped.
## Columns: probe = cpg_id, gene = target_gene, d = direction (-1/+1),
## p_plus = probability_plus1, tier = evidence_tier. Rows are ordered by
## probe, then gene (radix order); probe, gene and tier are factors. The
## provenance list "dmsa_direction_map" is attached with n_pairs and
## n_probes of the object itself.
##
## Run against the 0.99.4 lookups, this recipe regenerates the 2026.08 maps
## identically (identical() is TRUE for both). Map version 2026.11 is the
## same recipe on the cpgdirectionData 0.99.5 lookups (Zenodo record
## 24052448), which place every gene on its own chromosome at an hg19 TSS.
##
## Usage, from the dmsa package root:
##   Rscript inst/scripts/make_direction_maps.R [<dir with the lookup .rds>]
## Without a directory, the lookups are fetched with cpgdirectionData::cpgdData().

args <- commandArgs(trailingOnly = TRUE)
src_dir <- if (length(args) >= 1) path.expand(args[[1]]) else NULL

read_lookup <- function(name) {
  if (!is.null(src_dir)) return(readRDS(file.path(src_dir, paste0(name, ".rds"))))
  cpgdirectionData::cpgdData(name)
}

build_map <- function(lookup, meta) {
  L <- as.data.frame(lookup, stringsAsFactors = FALSE)
  r <- L[!is.na(L$direction), , drop = FALSE]
  r <- r[order(r$cpg_id, r$target_gene, method = "radix"), , drop = FALSE]
  out <- data.frame(probe = factor(r$cpg_id), gene = factor(r$target_gene),
                    d = as.integer(r$direction),
                    p_plus = as.numeric(r$probability_plus1),
                    tier = factor(r$evidence_tier))
  rownames(out) <- NULL
  meta$n_pairs <- nrow(out)
  meta$n_probes <- nlevels(out$probe)
  attr(out, "dmsa_direction_map") <- meta
  out
}

common <- list(build = "hg19",
               source = "cpgdirectionData 0.99.5 (catalogue lookups rebuilt 2026-11-08)",
               source_doi = "10.5281/zenodo.24052448",
               version = "2026.11")
layers <- list(
  blood = list(lookup = "lookup_blood_hg19",
               note = paste("Whole blood. The right proxy for an immune-dominated",
                            "sample - adult saliva by passive drool runs about",
                            "82% immune.")),
  epithelium = list(lookup = "lookup_nasal_epithelium_hg19",
                    note = paste("Nasal epithelium. The right proxy for an",
                                 "epithelium-dominated sample - neonatal buccal",
                                 "swabs run .89 to 1.00 epithelial, where blood",
                                 "calls would be the wrong tissue.")))

for (tissue in names(layers)) {
  meta <- c(list(tissue = tissue), common,
            list(n_pairs = NA_integer_, n_probes = NA_integer_,
                 note = layers[[tissue]]$note))
  m <- build_map(read_lookup(layers[[tissue]]$lookup), meta)
  f <- file.path("inst", "extdata", sprintf("direction_%s_hg19.rds", tissue))
  saveRDS(m, f, compress = "xz")
  cat(sprintf("%s: %d pairs, %d probes, %d genes -> %s\n", tissue,
              nrow(m), nlevels(m$probe), nlevels(m$gene), f))
}
