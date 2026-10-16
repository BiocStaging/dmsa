# `dmsa_change()` 0.99.18 statistical audit contract

This file records the inferential changes made after an adversarial audit of the
0.99.17 two-wave engine.

## Confirmatory guarantees intentionally narrowed

1. Before and after scores are built on one stacked-wave measurement scale.
2. Baseline adjustment is opt-in and is described as a different estimand.
3. Delta-cell/technical adjustment is optional and described as a
   composition/technical-adjusted estimand, not a cell-intrinsic effect.
4. Fixed two-person within-block roles are not sent through whole-block
   Freedman-Lane permutation. They require `paired = TRUE`.
5. Paired mode uses independent pair differences, HC3 studentisation, and a
   restricted-model Rademacher wild bootstrap. One multiplier stream is reused
   for every unit so maxT sees the observed cross-unit dependence.
6. Outcome-side change inference is flat-weight only until reliability weighting
   passes dedicated null-calibration simulations.
7. A CpG can map to multiple genes, but the same CpG x unit pair cannot be
   duplicated within a score.
8. Freedman-Lane is called permutation-calibrated, not unconditionally exact.

## Required validation before publication claims

Run `dmsa_change_null_calibration()` under at least:

- global null and partial null;
- n_pairs = 40, 60, 100;
- homoscedastic and role-specific heteroscedastic errors;
- pair ICC from near 0 to high;
- nuisance change correlated with role;
- baseline measurement error and baseline group differences;
- overlapping CpG-to-gene mappings;
- several family sizes and B values.

Report pointwise type-I error and maxT FWER with Monte Carlo confidence intervals.
Do not describe paired wild-bootstrap p-values as exact randomisation p-values;
mother/father roles are not randomized.
