# dmsa 0.99.17

## `dmsa_report()` now delivers the four structural methods

* **The same three things, for a result as for a frame.**
  `dmsa_mediate()`, `dmsa_modmed()`, `dmsa_transmit()` and
  `dmsa_concordance()` printed to the console and nothing else, so a
  result could be read but not filed. `dmsa_report()` is now an S3
  generic: handed any of those objects it writes `tables/` (csv, html
  and, with **gt**, docx), `figures/` and `summary.md` under `outdir`,
  exactly as the frame path has always done. The summary opens with the
  call that produced it, then the design, the tables inline, the figures,
  and how to read them. Nothing is recomputed: every number written was
  estimated by the method that produced the object.
* Tables: paths, indices, conditional indirect effects, measurement,
  bootstrap quantiles and module coverage (mediation / moderated
  mediation); per level, named units, measurement and every unit
  (transmission); per level, named units, step 0, surface and exposure
  (concordance). Figures: the path and index forests, the conditional
  indirect effect against the moderator, mother and father side by side
  per unit, and the largest within-couple r at each level against its
  batch-matched null. Viridis throughout.
* An existing report is never silently replaced (spec 54, as on the frame
  path), and the refusal now comes BEFORE anything is drawn, so a refused
  call leaves the directory as it was.
* `dmsa_report(frame)` is unchanged: the frame method is the function
  that was there, under a generic.

## Methylation columns are matched by canonical CpG id

* **A blocker, not a nicety.** The path functions matched a map's probe
  ids against `colnames()` verbatim. A frame whose columns carry study
  suffixes - `cg02187522_TC21_T1_AVP`, the shape of every Project Alpha
  build - matched NO pair against the live cpgdirection map (0 of
  12,498), and every one of the four methods stopped with "system 'hpa'
  has no usable module tone at 3+ probes". Only the bundled per-column
  snapshot happened to work, so the methods appeared to run while
  silently confined to it. `.paths_setup()` now reads every methylation
  column as its canonical CpG id and averages the replicate designs of
  one site (the frame's `collapse_site` policy) before any map is
  resolved; the map source line reports the collapse. On Alpha: 9,578
  cascade pairs, hpa 8 modules, oxytocin 741 probes over 9 modules.

## The registered recipe is now the shipped default (PI ruling, 2026-09-07)

* **`single_rel = 0.70`** in `dmsa_mediate()` and `dmsa_modmed()`, the
  registered choice (Savalei 2019), replacing `NULL` = estimate by
  Spearman-Brown. The estimate floors at .05 whenever a system's module
  tones do not cohere - which is the common case on real data - and a
  reliability of .05 gives an error variance of 19 and paths that cannot
  be read. `single_rel = NA` restores the estimated route.
* **`module_score = "majority"`** is now the default of `dmsa_modmed()`,
  `dmsa_transmit()` and `dmsa_concordance()`: the module's
  majority-direction probes, mean of the scaled M-values, signed by that
  direction - the v12 tone of the 2026i registered runs and of every
  number in the methods paper. `"aligned"`, the 0.99.12-16 default,
  remains available.
* Together these reproduce the registered runs from the package's own
  defaults. Registered run 3c re-estimated through `dmsa_modmed()`:
  measurement identical (hpa single, omega .70; noradrenaline
  congeneric, omega .57) and every standardized path within .0004 of the
  frozen result, every p within permutation noise.

## Two notes on the divergence this creates

* `dmsa_mediate()`'s all-systems branch keeps its own aligned tone builder
  (`dmsa_align` x `dmsa_scores`), so it and `dmsa_modmed()` no longer agree
  under their defaults; `dmsa_modmed(module_score = "aligned")` is its
  equivalent, and the reduction test now compares them that way. Unifying
  the two builders is a decision for the next version, not a silent change
  here: the engine (`.med_measure` / `.med_moments` / `.med_reg` /
  `.med_fl`) is untouched.
* `dmsa_report()` is a generic, so `progress` and `beep` are formals of
  `dmsa_report.dmsa_frame()` rather than of `dmsa_report()`. Every existing
  call is unaffected.

## Paths are reported on both scales

* The object now carries `scales`, the corrected SD of every modeled
  column, and the report's paths table adds `b_std` beside `b`: the
  engine estimates on the raw metric of the observed variables, while the
  registered runs and the paper print the standardized path. Both are
  now available without re-fitting.

# dmsa 0.99.16

## A CpG keeps every gene it pairs with (`dmsa_transmit()`, `dmsa_concordance()`, `dmsa_modmed()` / `dmsa_mediate()` maps; engine untouched)

* **Definition, not option.** DMSA takes a probe to every gene it pairs
  with in the direction map, because one methylation site can regulate
  several genes, in several systems (cpgdirection). The 0.99.14 path
  functions had kept one pair per probe within a system
  (`!duplicated()` in `.paths_tones()` and `.hier_scores()`), so a CpG
  paired with two cascade genes counted for the first only. Both now keep
  every pair at the gene level; the probe, module and system levels use
  each probe once (a probe is one unit, a module tone averages each of its
  probes once). `meta$pairs` (probe, gene, module) is returned; a probe's
  `probe_gene` label lists its genes as `A/B`. On the Alpha cascade this
  touches 2 noradrenergic and 30 immune pairs (0 in HPA, oxytocin,
  kynurenine, dopamine); the Alpha hierarchy results change at the immune
  gene level only (103 genes, was 101; nothing named before or after).
* Test: a synthetic map in which one CpG pairs with two genes scores the
  probe once and both genes (test-hierarchy.R).

# dmsa 0.99.15

## Planning ships with the engine: `dmsa_power_modmed()`, `dmsa_n_modmed()`, `dmsa_power_slope()` / `dmsa_n_slope()`, `dmsa_power_transmit()`, `dmsa_power_concordance()` / `dmsa_power_exposure()` / `dmsa_power_surface()` (POWER v2; engine untouched)

* **Through the estimator, not beside it.** `dmsa_power_modmed()` simulates
  the structural model (`M = a1 X + a2 W + a3 XW`, `Y = b1 M + c1 X + c2 W
  + b3 MW`, unit-variance latents, measurement error per role by its
  reliability, optional nuisance covariates and family blocks) and runs
  the package's own corrected-moments core (`.mm_core_v`) on every
  replicate. `mode = "t"` reads the engine's HC errors-in-variables t (the
  E15 pivot) against a normal reference - seconds per cell; `mode =
  "perm"` runs the shipped Freedman-Lane test with `B` permutations. Power
  per path, for the indirect effect (joint significance) and for the index
  of moderated mediation (first stage a3 b1, second stage a1 b3, dual stage
  the Holm union - the registered test); the mean corrected estimate
  against the truth and the conditioning-guard rate come with it.
  `quadratic = TRUE` prices the safeguard; `moderator = FALSE` is plain
  mediation. `dmsa_n_modmed()` bisects n for a target power of any
  estimand.
* **Closed forms** for a single (partial) slope between error-prone scores
  (Fisher z on the attenuated correlation, covariates cost df), and for the
  hierarchy: `dmsa_power_transmit()` (the paired child ~ mother + father
  slope at one level, with the power to NAME the unit within its
  level-local family - Bonferroni, the conservative side of maxT),
  `dmsa_power_concordance()` (within-couple r per unit, and named within
  the level), `dmsa_power_exposure()` (gamma3 as the Fisher-z difference of
  two within-couple correlations), `dmsa_power_surface()` (the corrected
  dyadic surface through the MATH 1.6 offsets: power of a4, the guard rate
  and the median inflation of the corrected curvature).
* **Reading.** Effect sizes are latent-scale; a .70 reliability
  attenuates the observed signal by sqrt(.70) per error-prone side, the
  estimator undoes it in the estimate and not in the information, so a
  latent r = .30 between two .70 tones is tested as an observed .21.
  POWER_v2.md carries the paper's tables at the 2026i sample sizes.
* **Tests.** `test-power.R`: attenuation / monotonicity / inversion of the
  closed forms, both simulation modes through the core, the level at the
  null, the hierarchy forms and the surface.

# dmsa 0.99.14

## The four path methods of paper 3, in the package: `dmsa_transmit()`, `dmsa_concordance()`, vector roles in `dmsa_modmed()` / `dmsa_mediate()`, measurement rule v1.2, the cpgdirection default (API proposal v1.1, PI rulings 0-5 of 2026-09-04; engine untouched)

* **`dmsa_transmit(child, parents, ...)` - system to system as the full
  DMSA hierarchy on both ends.** Not a system-score regression (PI rule
  2026-09-05): every unit of a system is scored at probe (the CpG M-value,
  observed), gene (the aligned `dmsa_scores()` expected score, composite
  at `composite_rel`), module (the module tone, composite) and system (the
  ladder) level on BOTH frames under the same measurement rule, and the
  estimand per unit is the paired slope `child_u ~ mother_u + father_u +
  covariates` (the child's contract plus the three arrays' global means).
  The test is the engine's HC errors-in-variables statistic (the E15
  pivot of `.med_fl()`) on ONE shared family-block permutation stream, so
  that maxT within the level-local family (all units of the system at
  that level, as `dmsa_levels()`) is exact - and, because the stream is
  shared, a tree-wide maxT over all units of all four levels, `p_tree`,
  gives one exact FWER per system for "anything named anywhere". Father
  beside mother throughout; the print says whether a named unit is a
  locus, a gene / module carried by one parent, or a tone, and reads a
  father effect of the mother's size or larger as "not intrauterine"
  (genotype, household or the paternal germline remain - the contrast
  cannot separate those in a pre-conception cohort, Kitaba 2023). Two
  frames enter as `dmsa_contract()`s (data, methylation, structural
  covariates, technical measure-covariates, sex, chip, tissue).
* **`dmsa_concordance(data, couple, partner_by, batch, ...)` - dyadic
  concordance at every level, and what it predicts.** Step 0: per unit at
  every level, the within-couple correlation against pseudo-couples
  re-paired WITHIN `batch`; the re-pairing is the permutation stream, so
  the level-local maxT and `p_tree` come from the same draws. `batch` has
  NO default (ruling 4): the partners of a couple share a processing
  plate, and a null that ignores it names the plate as concordance (the
  oxytocin "tone concordance" of the unstratified 2026i run was plate);
  `batch = FALSE` declares that none is shared. With `outcome`: the dyadic
  response surface per unit on stacked actor-partner rows and the
  incongruence-line contrast a4 = b3 - b4 + b5 by Freedman-Lane with the
  HC contrast pivot, observed at the probe level and corrected (MATH 1.6
  offsets, the conditioning guard, the uncorrected surface beside it) at
  gene / module / system. With `exposure` (a 0/1 column of the partner
  named by `exposed`): the MATH 7.3 slope change gamma3 per unit at every
  level, the exposure never permuted. Probe, module and system level
  reproduce the 2026i prototypes (`adv/hier/`) to the digit; at the gene
  level the package resolves the confidence of direction-only pairs
  through the pair ledger (as `dmsa_frame()`), where the prototype's
  interim map treated them as certain - the named genes are the same.
* **`dmsa_modmed()`: vector roles.** `x` and `m` accept vectors - several
  predictors and several PARALLEL mediators: one M-model per mediator
  (`M_j ~ all X + W + X:W ...`), one Y-model (`Y ~ all M + all X + W + M:W
  ...`), every `X_i -> M_j -> Y` product by joint significance with Holm
  over the products as the declared family (`p_holm_family`), the index
  of moderated mediation per pair, and - with two or more LATENT
  mediators - the connected measurement of the mediators (MATH 1.7: the
  correlated-error covariance theta_MM of the two Bartlett scores from a
  joint two-factor fit, lavaan, subtracted from the (M_i, M_j) cell when
  |theta| / sqrt(V_i V_j) > .10; `theta_mm = "none"` skips it). The vector
  core `.mm_core_v()` reproduces the scalar core to the digit when `x` and
  `m` are scalars (tested), so every 0.99.12-0.99.13 result is unchanged
  under the old rule. `dmsa_mediate()` hands several mediators over to
  `dmsa_modmed()` (several X were always native there).
* **Measurement rule v1.2 is the default (ruling 2).** In `dmsa_modmed()`,
  `dmsa_mediate()` and the two new functions, every score is chip
  quasi-demeaned and residualised on `measure_covariates` (the technical
  terms; `NULL` = `covariates`) PLUS the `sex` column PLUS the array's
  global mean M-value (`global_mean = TRUE`; over the `methylation` matrix,
  or over every cg-site column of `data` - the output says which), and sex
  and the global mean enter every equation as covariates. The saliva-wide
  common factor of the 2026i builds is sex + global mean (MAP_V5_NOTE 6;
  share .74 -> .31 after residualisation). `sex = NULL` /
  `global_mean = FALSE` reproduce the 0.99.13 rule; `$rule` prints what ran.
* **Direction source (ruling 3).** `map = "auto"` is the default of the path
  functions: cpgdirection pair directions when the package is installed
  (as `dmsa_frame()`), otherwise the bundled snapshot - now
  `inst/extdata/coverage_v5_cpgdirection.csv.gz`, regenerated from
  cpgdirection 2.6.0 for the whole Alpha cascade in blood and nasal
  epithelium (9,578 / 9,574 usable pairs of 16,823 cascade CpGs), so the
  fallback is no longer the 3.7 % interim snapshot. `$map_source` and
  `$map_note` record which. `pairs = "cascade"` (default) keeps only the
  selection cascade's CpG x gene assignments - a probe that pairs with a
  second reference gene is not counted twice in the paired tests;
  `pairs = "all"` keeps every pair as `dmsa_frame()` does.
  `dmsa_frame()`'s own `map = "alpha"` snapshot (`coverage_v4_full.csv`)
  is untouched in this release; swapping it is a separate item because a
  dozen test files build their synthetic frames from it.
* **The sensitivity battery in the call (ruling 5).** `dmsa_modmed(...,
  sensitivities = c("homoscedastic", "unrestricted", "composite_rel",
  "quadratic", "global_mean_off"), B_sens = 1999)` re-runs the same call
  with one change each and returns one table (`$sensitivities`: model,
  term, estimate, p) - the registered pattern of every 2026i PoC.
* **Module tone builder.** `module_score = "aligned"` (the package rule:
  the direction-weighted expected score, unchanged) or `"majority"` (the
  v12 tone of the 2026i proofs of concept: the module's majority-direction
  probes, mean of the scaled M-values, signed by the majority direction;
  `anchor` orders a reference module first). The registered 2026i runs
  used "majority"; the package default stays "aligned".
* **Engine.** `mediate.R`'s `.med_measure / .med_moments / .med_reg /
  .med_fl` are untouched. `.hier_fl()` re-computes `.med_fl()`'s HC
  statistic against an EXTERNAL permutation stream so that the units of a
  level share it; it reproduces `.med_fl()`'s b, t and p exactly under the
  same seed (tested).
* **Tests.** `test-hierarchy.R` (the shared-stream identity, both
  functions at every level, the batch refusal, the snapshot),
  `test-modmed-vector.R` (the vector/scalar identity, parallel mediators,
  rule v1.2 and its switches, the direction source, theta_MM, the
  battery). lavaan added to Suggests.

# dmsa 0.99.13

## `dmsa_mediate()` / `dmsa_modmed()`: heteroscedasticity-consistent pivot (E15), permutation strata (E15), empty reduced model (E16)

* **What changed (E15, pivot).** The Freedman-Lane statistic is now
  studentised with the heteroscedasticity-consistent errors-in-variables
  sandwich by default. The corrected-moments estimator solves the
  estimating equation `sum_i x_i (y_i - x_i'b)/(n-1) + Ofs b = 0`, with
  `Ofs` the moment offsets; its per-person score is `s_i = x_i e_i + Ofs b`
  and `Var(b) = S*xx^-1 [sum_i s_i s_i'] S*xx^-1 / (n-1)^2` (HC0 form;
  Buonaccorsi 2010 ch. 5). The permutation test is the studentised
  Freedman-Lane of DiCiccio & Romano (2017): the meat is recomputed on every
  permuted response, the reduced-model residuals, the corrected `S*xx` and
  the block exchange are unchanged. With no measurement error and
  homoscedastic disturbances it is the ordinary HC0 statistic. `pivot =
  "homoscedastic"` selects the E14 sandwich of 0.99.11-0.99.12 and
  reproduces those results exactly (tested).
* **Why.** Gate v3 found the E14 pivot at .087-.097 on the interaction path
  when the mediator's disturbance variance grows with the moderator - the
  heterogeneity that random slopes in W produce, i.e. the situation in
  which a moderation claim is being made. Under the HC pivot that cell
  holds .053; under homoscedasticity nothing changes (.047 / .057 at the
  null, power .97 / .98 unchanged); under heteroscedasticity the
  interaction path gives up six points of power. GATE_REPORT_hc (prototype)
  and GATE_REPORT_E15 (the shipped engine, full gate v3 grid, the v2/E14
  gate, and n = 70 cells).
* **What changed (E15, strata).** `dmsa_mediate()` and `dmsa_modmed()`
  take `strata =`, a column naming a permutation stratum. With `block`,
  equal-size blocks exchange only within a stratum (a block spanning two
  strata is refused); without `block`, residuals permute within stratum.
  Registered use: the after-wave of the before/after builds
  (GATE_REPORT_delta section 1: stratification alone keeps the test
  honest, the wave dummies remove the bias; both are required).
* **What changed (E16).** An empty reduced model - one predictor and no
  covariates in an equation - is fitted as zero (its residuals are the
  centered response) instead of failing in `.med_reg()`. Plain mediation
  with no covariates now runs; no existing result is touched, because that
  case previously errored. The refusal `dmsa_modmed()` carried in 0.99.12
  is gone.
* **Reading the change.** Estimates are identical between the two pivots;
  only the studentisation differs, so only permutation p-values move.
  Under homoscedastic disturbances the difference is Monte-Carlo noise;
  under heteroscedastic ones the HC p is the honest one.
* Tests: E15 (HC0 identity for observed regressors; the score-sandwich
  formula with error; default; homoscedastic pivot retained and distinct),
  strata (nesting, refusal, within-stratum without blocks, pass-through),
  E16 (empty reduced model runs). The E14 tests now call the homoscedastic
  pivot explicitly.

# dmsa 0.99.12

## New: `dmsa_modmed()` - general moderated mediation over observed and latent roles

* **What it adds.** The general form of `dmsa_mediate()`: each of the
  predictor, mediator, moderator and outcome is EITHER a methylation system
  (a latent tone by the measurement ladder, or a single-indicator composite
  at a fixed reliability - .70 is the registered choice) OR a column of
  `data` (observed, or a questionnaire composite at a fixed reliability).
  The moderator can act on the a-path (Hayes 7), the b-path (14) or both
  (58, dual stage), with optional direct-effect moderation (8 / 15 / 59)
  and an optional curvilinear safeguard (`quadratic = TRUE`: x^2, w^2 in
  the M-model and m^2, x^2, w^2 in the Y-model with their fourth-moment
  offsets). Latent outcomes are allowed (the outcome's error variance
  enters the moment correction). The dual-stage moderation is tested by the
  Holm union over the joint-significance p of the three products a3*b1,
  a1*b3, a3*b3 - the difference of the conditional indirect effect at
  +/- 1 SD of W is zero iff all three are. The full-pipeline family-cluster
  bootstrap holds every system's measurement rung at its full-sample value
  (`boot_rung = "fixed"`) and reports percentile and basic intervals.
* **What did not change.** Nothing in `mediate.R`: `dmsa_modmed()` is built
  on the same `.med_measure()`, `.med_moments()`, `.med_reg()` and
  `.med_fl()` (E14 sandwich pivot), and with every role a system,
  `w_on = "a"` and `quadratic = FALSE` its numbers are exactly those of
  `dmsa_mediate(model = "7")` - tested to 1e-12 on the coefficients and
  bit-for-bit on the permutation p-values. Existing results are unaffected.
* **Why.** The advanced-uses paper needs mixed observed/latent designs
  (a questionnaire moderator on a latent mediator, a methylation outcome,
  a dual-stage model) that `dmsa_mediate()` cannot express; the general
  engine and its tests were derived in MATH_v1_general_modmed and gated
  (GATE_REPORT_v3: level holds in every cell type; an omitted curvilinear
  main effect produces a spurious interaction 55% of the time, which
  `quadratic = TRUE` removes).
* **Engine note (E16, not changed).** Plain mediation (`w = NULL`) with a
  single predictor and no covariates leaves the a-path's reduced model
  empty, which `.med_reg()` does not fit; `dmsa_modmed()` refuses with a
  message asking for a covariate. `dmsa_mediate()` has the same edge.
* Tests: `tests/testthat/test-modmed.R` (52 expectations).

# dmsa 0.99.11

## `dmsa_mediate()`: the permutation pivot is now the errors-in-variables sandwich (engine finding E14)

* **What changed.** `dmsa_mediate()` estimates structural coefficients from
  measurement-error-corrected moments, and tests each path by Freedman-Lane
  permutation. The permutation statistic previously divided the CORRECTED
  coefficient by the plain-OLS residual variance (`Var(y) - b'Sxy`), which
  over-subtracts once the coefficient is disattenuated. The level was exact
  (the same statistic is computed on observed and permuted responses), but
  the null distribution grew heavy tails and power on error-prone
  main-effect paths was lost - the anchor `b[m]` and the `a[x]` path above
  all. The statistic is now the errors-in-variables sandwich,
  `Var(b) = sigma2_obs * S*xx^-1 S_obs S*xx^-1 / (n - k)` with
  `sigma2_obs = Var(y) - 2 b'Sxy + b'S_obs b` (Meijer, Spierdijk & Wansbeek
  2020; Devlieger, Mayer & Rosseel 2016), held fixed across permutations
  like the corrected design moments. Freedman-Lane geometry (raw
  reduced-model residuals, family-block exchange, within-chip permutation)
  is unchanged. For error-free regressors the new pivot IS the classical
  OLS t (tested exactly).
* **Effect.** Identical generator, n = 200, reliability .70 on X, M and W,
  B = 499: power `a[x]` .855 -> .93, `b[m]` .735 -> .82, interaction
  paths unchanged (.76 / .87), Type I on the interaction paths at the
  null .047 / .043, zero degenerate permutations (the old pivot's
  `rssp <= 0` clamp no longer fires). The v2 simulation gate was rerun
  under the new pivot (corner Type I and power surfaces; see
  GATE_REPORT_E14).
* **Scope.** `.med_fl()` is used only by `dmsa_mediate()`. Nothing else in
  the package - `dmsa_triangulate()`, the lenses, maxT/minP,
  `dmsa_moderate()`, `dmsa_change()`, the outcome families - computes a
  corrected coefficient, and none of their numbers move. The `Vy` argument
  of the internal `.med_fl()` is retained for its callers but no longer
  enters the pivot.
* Tests: two new tests in `test-mediate.R` (exact OLS reduction for
  observed regressors; sandwich formula for an error-prone regressor; no
  degenerate permutations and plausible level at a true null).

# dmsa 0.99.10

## `blocks =` replaces `random_effects =`, and there is no default (spec 42)

* **The exchangeability blocks are declared with `blocks = `.** The old
  name `random_effects = ` always meant the blocks - the rows that travel
  together under permutation - not lme4 random effects; it is kept as a
  deprecated alias (a message says so), and giving both is an error.
* **No default blocks (PI ruling, 2026-09-02).** `random_effects = "cID"`
  used to be the default, which hard-errored on any dataset without that
  column. Now, if you declare none, `dmsa_frame()` says live that the
  permutation null treats every row as independent and asks you to make
  sure that is intended (couples, family members and repeated measures
  need `blocks = "<id column>"`); the corrections table records "none
  declared". A block column you DID name but that is missing from the
  data is still a hard error (spec 43). `covariates = "contract"` carries
  the Alpha build's block structure and is adopted, with a message, when
  you give none.

## Unresolved polarity is counted and said; user polarity merges (spec 40, E13)

* **A user polarity table on the bundled Alpha cascade now MERGES.**
  `dmsa_sets("alpha", polarity = <table>)` used to take the table as the
  cascade's entire polarity, so a one-gene override left every other gene
  without a sign, and the report then weighted those genes 0 at the system
  level with nothing said. It now behaves as `dmsa_align(polarity = )` has
  since E4: user rows override gene for gene (grade "user"), the curation
  fills the rest, one message states the counts. A user cascade keeps its
  own table as its polarity, unchanged.
* **What the system tone is built from is stated.** Per system:
  how many testable genes carry a signed polarity, how many are explicitly
  off-axis, and how many are UNRESOLVED (no entry or NA) - those are
  weighted 0 in the system score, never assumed +1, and now the frame
  print says so naming the genes, summary.md states the composition of
  every system tone, `tables/polarity_audit.csv` lists the per-gene status,
  and `dmsa_coverage(level = "system")` gains
  `n_genes_polarity_signed` / `n_genes_polarity_unresolved` /
  `genes_polarity_unresolved`.

## The statistic columns, labelled (spec 26/28)

* summary.md's reporting block ends with a glossary of every p-value column
  in `tables/units.csv`: what it tests, which multiplicity it is corrected
  for, and whether it names a finding. Only the per-lens family-adjusted p
  names (any lens below alpha within its own family; `selected` is the
  verdict); the calibrated union and the family-corrected omnibus confirm
  (the dagger tiers); the raw omnibus and the joint union are synthesis -
  quoted, never named by.

## Report failures are recorded, never swallowed (spec 53)

* **`tables/qc_report.csv` and a "Report QC" section in summary.md.**
  Thirteen places in the report used `tryCatch(..., error = function(e)
  NULL)`: a component that failed simply vanished, and a reader could not
  tell "absent because not applicable" from "absent because it broke".
  Every one now goes through a recording wrapper: the failure is logged
  with its component, outcome, level, unit, what took its place, and the
  error message; the analysis continues exactly as before. The ledger is
  always written when tables are on (a header-only file states that
  nothing failed), listed in summary.md, counted once on the console, and
  carried on the report object as `$qc`. Fallbacks that were already
  announced - the locus panel's bare-panel fallback, a gene model that
  could not be fetched - are recorded too, as "fell back".

## The locus panel draws the tested estimator (spec 52)

* **One estimator for the figure and the test.** The locus panel's per-probe
  bars used to come from a report-side least-squares fit of each probe's
  raw M-value on the outcome, while the gene test pooled per-probe fits of
  the winsorised, standardised probe on the standardised term (with the
  random-intercept transform when a chip was declared random). Two
  estimators, one figure, no label. `dmsa_triangulate()` now hands its own
  per-probe fits over (`attr(., "probe_fits")`: column, b, se, plus the
  residual df and the RI gamma), the panel draws exactly those, the
  greying p is the nominal p of that same fit, and the panel's context line
  says so. The raw-M OLS remains only as a labelled fallback for a panel
  drawn without a test.
* **`tables/locus_probes.csv`** lists every number a locus panel drew - one
  row per probe with b, se, z, nominal p, the estimator and its scale - so
  figure and test can be audited against each other row by row; summary.md
  points to it.

## Regression tests for the lenses, custom polarity and the display rules (spec 62/63/67)

* **The three lenses are now tested separately** (`test-lenses-separate.R`).
  Planted shapes calibrated by simulation pick out a lens by its own
  definition: a sign-cancelling signal is seen by the diffuse lens alone; a
  dense, tiny, consistent signal by the directional lenses and not the
  diffuse one; pure noise by none. At the report level the same shapes pin
  the carrying-lens label in `units.csv` and the "carried by the ... lens"
  sentence in summary.md - the sentence that goes into papers.
* **A user polarity table runs end to end** (`test-polarity-custom.R`):
  flipping every gene's polarity mirrors the system level exactly (opposite
  direction, identical p under the same seed - a RULE ZERO invariant);
  polarity touches only the system level (gene and probe results are
  identical between runs); the frame carries the user's table gene for
  gene; and each `missing_polarity` policy does what it says.
* **The moderation display rules are pinned** (`test-figure-displays.R`):
  a two-level moderator draws one panel with each level's own b and p in
  the legend and no Johnson-Neyman continuum; a two-level outcome is shown
  as fitted probabilities on the log-odds scale; a continuous moderator
  keeps the quantile ladder and the JN band; the overview carries only its
  title. Figures are read back from uncompressed PDF text, so the tests
  assert on what the figure says, not on pixels.
* **summary.md is written as UTF-8 in every locale.** Under a C/POSIX
  locale the badge daggers came out as `<U+2020>` / `<U+2021>`;
  `writeLines()` translates to the native encoding first, even through a
  UTF-8 connection, so the bytes are now written as they are.

# dmsa 0.99.9

## The all-one-direction note tells the user where to check

* **"Every one of N direction calls is -1" no longer distrusts
  cpgdirection.** When the calls carry evidence tiers (resolved by
  cpgdirection), an all-one-direction set is usually genuine curation:
  the old warning is now an informational message that says so and names
  exactly where to verify - `tables/analysis_set.csv` (columns
  `direction`, `direction_tier`, `evidence`) and `frame$map`. Calls with
  NO evidence tiers (a user-supplied map or bare `d` vector) keep an
  actionable warning that names the check
  (`table(map$best_direction)`). Either note is shown once per probe set
  per session, not once per lens battery.

## Binary moderator: one panel; leaner overview; stale frames

* **A two-level moderator draws ONE panel.** The second panel is gone
  entirely (a significant interaction means, by definition, that the two
  slopes differ): each level's own slope and its significance now sit in
  the legend of the single panel ("Husband [1]: b -1.507, p < .001"),
  the caption states that the tested interaction IS the difference
  between the slopes, and the device narrows so the panel keeps its
  proportions.
* **The overview figure carries only its title.** The explanatory
  subtitle is removed; the visual codes (bold, shading, daggers) are
  defined once in summary.md.
* **Stale frame objects are detected and worked around.** A frame built
  by an older dmsa installation inside a long R session carries stale
  stored fields - the pills battery drew a linear display for a 0/1
  outcome exactly this way. dmsa_report() now compares the frame's build
  stamp with the running installation and says plainly to rebuild the
  frame; independently, the moderation figure re-checks two-levelness on
  the data itself, so the binary display appears even from a stale frame.
* The empty brain-bridge line prints "none recorded" instead of
  "ensemble percentiles NA-NA".

## Moderation joins the report; Ensembl fetch repaired

* **Moderation no longer replaces the main analysis.** `moderation = TRUE`
  used to run ONLY the interaction battery - no system, module or gene
  results, no overview or locus figures. Both batteries now run in one
  report (the price is runtime: both permute).
* **The gaussian gate applies only where it is true.** Moderation with
  `frame_role = "outcome"` fits `S ~ outcome x moderator` - the TONE
  SCORE is the response - so a two-level outcome on the right-hand side
  is an ordinary group contrast and is no longer refused. The gate
  remains for `frame_role = "predictor"`, where the outcome is the
  response of a linear model.
* **The moderation figure is binary-aware.** A two-level outcome gets a
  LOGISTIC display: panel 1 draws fitted probability curves per moderator
  level over the jittered 0/1 observations ("P(taking pills [1])"),
  stated as a display of the same interaction - the tested statistic
  stays the permutation-calibrated linear product. A two-level moderator
  gets no Johnson-Neyman continuum (a JN band over sex = 1..2 sweeps
  impossible values): panel 2 shows the two simple slopes with 95% CIs,
  whose difference IS the tested interaction, with a quasi-separation
  warning when the log-odds scale is unstable.
* **Console survivors deduplicate**: a probe that also survives as its
  gene's or module's only member prints once, tagged "same data at N
  levels" instead of three near-identical lines.
* **Ensembl fetch repaired**: lookup/symbol serves json/xml/jsonp ONLY -
  the old request for text/x-gff3 was an unconditional HTTP 400 from
  Ensembl on every machine, which surfaced the moment gene_models =
  "auto" first exercised the path. The lookup now asks for JSON and
  reads its four flat fields with no new dependency; the exon/CDS
  feature call stays GFF3 (overlap/region supports it). Exon models
  should now actually draw on connected machines.
* `palette` is documented: the viridisLite palettes (viridis, magma,
  plasma, inferno, cividis, mako, rocket, turbo) with colour-blindness
  and greyscale notes.

## Covariate/moderator overlap said out loud

* **A column given both as a covariate and as a moderator is flagged
  live.** With `moderation = TRUE` it is dropped from the covariates
  (the moderated model already carries its main effect once) with a
  console message and a corrections note saying how to undo it; with
  `mod = ` given but `moderation = FALSE`, a live message states that NO
  moderation is tested and the column stays a plain covariate - the
  "moderators ignored" fact used to live only in the corrections table,
  pages away from the user who typed `mod = "sex"` meaning to test it.

## Labels that cannot lie, levels that cannot flip

* **`predictor_labels` now names the `predictors = ` column** - the only
  reading a human makes of it. It used to silently label the covariates
  (a legacy of the old orientation vocabulary), which errored the PI's
  own call with "1 label(s) for 6 column(s)". Covariate display names
  have their own argument, `covariate_labels`; using `predictor_labels`
  with `outcomes = ` errors and points to the right spelling.
* **Level labels can be matched BY VALUE.** Every level-label pair -
  `outcome_levels`, `predictor_levels`, and the new `mod_levels` for a
  two-level moderator - is safest as a named vector:
  `c("1" = "Husband", "2" = "Wife")`. Names are matched against the
  variable's actual two levels whatever the coding (1/2, 0/1, "a"/"b"),
  in any order, so a positional flip is impossible; a name that is not
  one of the real levels errors immediately, naming them. Unnamed pairs
  keep the documented sorted-value order (low, high).
* **A two-level moderator gets a two-line slopes panel.** The quantile
  ladder used to draw a third line at an impossible middle value
  (sex = 1.5); a moderator with exactly two observed levels now draws
  one line per level, labelled by `mod_levels` when declared. A length-2
  `mod_label` (one display NAME, not level labels) errors and points to
  `mod_levels` instead of silently keeping only the first entry.

## Predictor_levels

* **Level labels in the user's own vocabulary.** A frame declared with
  `predictors = "pills_past_T1"` labels its levels with
  `predictor_levels = list(pills_past_T1 = c("not taking pills", "taking
  pills"))` - the user who thinks of pills as a predictor is never asked
  to reach for an "outcome" argument. `outcome_levels` stays the spelling
  for `outcomes = `; giving both is refused in plain words, and each
  spelling requires its matching declaration so a mismatch is caught at
  the call, not discovered in the report.

## The family-corrected omnibus and its badge

* **New: `p_omnibus_adj` - the ACAT omnibus, family-corrected.** The
  omnibus already absorbs the cross-LENS multiplicity (ACAT is a valid
  combination under arbitrary dependence; Liu & Xie 2020); this corrects
  it across the FAMILY too, by permutation minP on the same stream with a
  symmetric leave-one-out p convention on both sides. The result is an
  exact family-wise p (measured FWER .050 in the undermining sims,
  dmsa_patch/undermine_omnibus.R) that defends against BOTH
  multiplicities and is most powerful precisely where all three lenses
  agree - the signal class the union test is weakest for. The calibration
  distribution rides along as attr "omnibus_null_min". Additive: nothing
  previously reported changes value. (This corrects an earlier note that
  a family-corrected omnibus was unusable at feasible B: the granularity
  floor affects only ~3G/B of null draws, not all of them, and the sims
  confirm exactness.)
* **The OMNIBUS-CONFIRMED badge (double dagger).** A named unit whose
  family-corrected omnibus clears alpha is marked alongside the union
  badge: dagger = exact-confirmed (strongest single lens), double dagger
  = omnibus-confirmed (all three lenses combined) - two exact
  family-wise claims riding on top of the pre-registered naming, each
  labelled by what it rests on. Gene table gains an "omnibus adj"
  column; the Results template reports omnibus raw AND family-corrected
  with their scopes in words; console print and the overview labels
  carry the badge.
* Undermining recorded: analytic corrections (Bonferroni/Holm) on the
  ACAT omnibus are valid but ~5x conservative under the battery's lens
  dependence and granularity-blind to single-lens signals; the
  permutation calibration is exact and adaptive. On the pills battery:
  PRLR family-corrected omnibus ~.005 and GREB1 ~.012 (both
  omnibus-confirmed), while the analytic route gives .0315 / .0795.

## Labels, precise omnibus wording, exons back by default

* **Binary outcomes get real labels.** `outcome_levels` now also accepts a
  named list - `list(pills_past_T1 = c("not taking pills", "taking
  pills"))` - so every outcome of a multi-outcome frame can label its two
  levels. The report's direction sentence becomes "Taking pills [1] was
  associated with methylation consistent with a HIGHER expression tone of
  GREB1, relative to not taking pills [0]" - declared label first, coded
  value in brackets so the sentence is self-verifying against the data.
  Two-level outcomes are now detected as binary whatever their storage
  type (numeric, logical, factor, or character); a factor-coded 0/1
  column used to fall through to the "Higher X was associated ..."
  continuous wording.
* **The omnibus says exactly what it defends against.** The Results
  template now states, beside every ACAT omnibus value, that it is one
  valid test combining the three lenses (the cross-LENS multiplicity is
  absorbed by the combination) but is NOT corrected across the gene
  family - and that the exact union p beside it is corrected across both.
  A raw .0017 can no longer read as a family-level claim.
* **Exon models are back by default.** `gene_models = "auto"` (the new
  default) fetches the exon/TSS model from Ensembl for NAMED genes only -
  a handful of small announced requests - and degrades to the true
  coordinate axis with a clear message when offline. TRUE forces the
  fetch, FALSE never touches the network, and a prebuilt
  `dmsa_gene_model()` table still works fully offline.

## Exact second-level calibration, measured error, the badge

* **New: second-level Westfall-Young minP calibration of the any-lens
  rule** (Westfall & Young 1993; Ge, Dudoit & Speed 2003), computed from
  the same permutation stream at no extra permutation cost. Every
  triangulation now returns `p_union_exact` - an EXACT family-wise p for
  the claim "some lens carries this unit", adaptive to the measured lens
  dependence (perfectly correlated lenses cost nothing extra; independent
  ones the full threefold - "correction per what is needed") - plus
  `attr(., "union_null_min")`, the calibration distribution itself.
* **The run's own realized family-wise error is measured and stated.**
  `fwer_realized` per family = the share of null draws whose best
  adjusted p anywhere in the family clears alpha - the design-specific
  number behind the MS's ".04-.12" disclosure. summary.md states it per
  gene family; the Methods template cites the construction.
* **The EXACT-CONFIRMED badge (two tiers, PI ruling).** Naming stays the
  pre-registered any-lens rule; a named unit whose exact union p also
  clears alpha is additionally marked exact-confirmed - dagger in the
  gene table and overview figure, "EXACT-CONFIRMED" in the prose and
  console. The printed honesty line beside every named unit is now the
  exact union p (directly comparable to the lens adj it sits beside);
  `p_unit_adj` remains computed and exposed in units.csv. Nothing
  previously reported changes value.
* Undermining recorded (dmsa_patch/undermine_calibrated.R): the exact
  calibrated union's power matches the joint max in every regime (the
  ratio arithmetic .05 x .05/FWER UNDER-corrects - measured threshold
  .020 where the ratio predicts .029 - so the permutation quantile, not
  the ratio, is the only defensible calibration), which is why it is a
  badge and not the naming rule.

## The rule carried through every output

* **Every layer of every report follows the naming rule, and a named gene
  always gets its plot.** The system-level summary line now leads with the
  same naming statistic as every other level (best lens's family-adjusted
  p, star = named; the ACAT value is labelled "omnibus raw" and the joint
  union p rides along). The locus panel of a named gene is guaranteed: if
  the full panel (gene model, genomic coordinates) fails to draw, the
  report falls back to the bare, always-drawable panel instead of
  silently skipping; summary.md points to each named gene's panel by file
  name, and in the (now rare) case a panel could not be drawn at all,
  summary.md says so beside the findings instead of leaving a hole.

## The naming rule, re-ruled

* **A finding is NAMED by the any-lens rule again - the manuscript's
  pre-registered rule.** A gene (or module, or system) is named when any
  single lens's family-adjusted p falls below alpha, with the carrying
  lens reported; the realized family-wise error of this union rule is
  disclosed (the MS benchmark measured .04-.12; ~.09 under the Alpha
  battery's lens dependence) rather than exact. The 0.99.9 interim rule -
  gating on the exact joint any-lens union test (`p_unit_adj`) - was
  re-ruled by the PI after simulation showed every exact-alpha
  alternative (joint max, continuous ACAT, Stouffer, hybrids) costs ~10
  power points in every signal regime; the extra power of the any-lens
  rule IS its disclosed inflation, and rank-p ACAT under maxT carries the
  granularity floor the MS already documents. `p_unit` / `p_unit_adj`
  remain computed, exposed in the units table, and printed beside every
  named unit as the honesty line; they no longer gate. Simulation
  scripts: `undermine_naming.R`, `undermine_naming2.R` (in the patch
  archive).
* **One statistic everywhere.** The overview figure's ordering, bolding
  and near-threshold shading, the summary's gene-results table ordering
  and bolding, the Results/Methods templates, and the console print are
  all keyed to the same naming statistic (the best lens's family-adjusted
  p). The gene table's raw ACAT column is now headed "omnibus (raw)" so
  an uncorrected number cannot read as a verdict, and the Methods
  template states the naming rule with its disclosed error rate.
* **The probe-QC exclusion count is scoped to the frame.** The console
  line and corrections note now lead with how many excluded CpGs map
  into the systems being analysed, with the whole-file count as context
  ("0 CpG(s) mapping into the analysed system(s) excluded ... (486
  across the whole file)"), resolved by one extra flag-mode
  `cpg_gene_pairs()` call over just the excluded CpGs; when their mapping
  cannot be known, the note states the file-wide scope in words.
  cpgdirection's own file-scoped message is muffled in favour of the
  scoped one.

## Report figures and the gene table

* **Overview figures follow the figure rules (spec 29-38).** One panel per
  SYSTEM, in that system's own accent colour, at most 40 units per page
  with `_p2`, `_p3`, ... pagination (a real Alpha battery had put ~70
  genes from 2 systems in one uncoloured 6000-px panel). Within a system,
  units are ordered by the joint any-lens adjusted p; rows at
  `p_unit_adj < 0.20` get a light background shade; survivors stay bold.
  File names carry the system slug and page number only when needed.
* **summary.md now states the gene results.** A `### Gene results` section
  lists, per outcome x system, the top 10 genes by the joint test - CpGs
  tested, concordance, direction, best lens, lens-adjusted, unit-adjusted
  and omnibus p - with survivors in bold and a pointer to the full
  `tables/units.csv`. A summary that jumped from "no survivors" to the
  module table left the reader with no idea what the genes DID.
* **Gene-significant, probe-silent loci are drawn and explained.** The
  report now passes each probe's nominal p (from the same fit as the
  panel's effects) into the locus figure, so CpGs above alpha are greyed;
  when NO probe clears alpha on its own, the panel says so under its
  title and states that the gene-level result pools small,
  direction-consistent shifts across the CpGs. summary.md carries the
  same explanation once, at the top of the gene results.
* **New: `dmsa_save_analysis_set()`.** One row per CpG x gene pair that
  actually entered the analysis - probe, gene, system, module, direction
  call with tier and evidence, co-effect flag, hg38 coordinates when the
  cascade covers the probe, and the frame's direction source, tissue and
  reference - written for the supplement of a paper and for user
  validation. Every report writes it automatically as
  `tables/analysis_set.csv`; the pair ledger keeps being the superset
  that also lists every pair NOT used and why.

## Orientation spellings

* **`outcomes =` and `predictors =` replace `outcome =`.** The phenotype
  columns are now named `outcomes` (which column(s) of `data` the
  methylation is tested against), or - when the question runs the other
  way - `predictors`, which also sets `frame_role = "outcome"`
  automatically, so the reversed orientation is one argument instead of
  two; a contradicting explicit `frame_role = "predictor"` alongside it
  errors in plain words. The former `outcome` argument is OUT of the
  signature but still accepted for backward compatibility (it means
  exactly `outcomes`), so no existing script breaks. Two side effects of
  the compatibility slot: every argument after `data` and `methylation`
  must be given by its full name (abbreviated argument names no longer
  match - arguably a feature in analysis scripts), and an unknown argument
  is an immediate error naming it.

* **`random_effects` now speaks lme4.** State your grouping factors the way
  you think about them. With TWO factors `dmsa_frame()` detects the
  structure and announces it in lme4 notation: NESTED (couples on the same
  chip, `(1 | chip/cID)`) keeps the existing block behaviour exactly, with
  the coarse factor's variance component left to `chip =`; CROSSED
  (partners on DIFFERENT chips, `(1 | cID) + (1 | chip)`) routes the
  finest dependence factor to the exchangeability blocks and the other
  into the `(1 | .)` random intercept - the case that used to stop the run
  with an all-singleton block error. Aliased factors are detected; three
  or more are refused; a conflicting explicit `chip =` or
  `chip_effect = "none"` errors in plain words. This settles spec item 42
  ("random_effects doesn't mean what labs think"): the name now means what
  labs think.
* The all-singleton exchangeability-block error now diagnoses itself: with
  crossed block columns it names which column alone WOULD permute (e.g.
  couples on different chips make `cID x chip_T1` all-singleton while `cID`
  alone is fine) and points a chip/slide/plate/batch-named column at
  `chip = `, where it enters as a random intercept.

## Group selection

* **New: CpG-set selection by name pattern.** Methylation files often carry
  several biological groups in one table (in Alpha's child file,
  `_maternal` and `_paternal` tag the parents' probes; the child's carry no
  tag). `dmsa_frame()` gains `cpgs_include` / `cpgs_exclude` - fixed
  substrings matched against CpG column NAMES, include first, then exclude -
  so analysing the child is one line:
  `cpgs_exclude = c("_maternal", "_paternal")`. Only CpG-site columns are
  candidates (outcomes and covariates are never touched), the drop is
  recorded in the corrections table, patterns that remove everything are a
  hard error, and both the pair path and the bundled path honour the same
  filter. The standalone helper `dmsa_cpg_columns(x, include=, exclude=)`
  applies the identical rule anywhere a column vector is needed.

## Engine review pass — PI-ruled fixes E1-E10

Ten findings from the whole-package adversarial review were ruled on
individually by the PI and fixed. These are the first changes to engine
files in this release cycle; each was approved explicitly.

* **E1 (fit.R)** A vestigial assignment inside the mixed-model loop
  reassigned the focal-term index, so whenever lme4 dropped or reordered a
  fixed effect, the permutation null pooled a DIFFERENT coefficient than the
  observed statistic. Removed; in every run where term orders agreed (the
  usual case) no number changes.
* **E2 (report)** A unit is now NAMED a finding by the joint any-lens test
  (`p_unit_adj` - the max over the three lens statistics put through the
  permutation machinery as one object, so the any-lens rule itself carries
  exact family-wise control). Per-lens survivors remain reported as a
  labelled second tier ("lens-level survivors"). `selected` now keys on
  `p_unit_adj`; prose, print and the how-to-report templates follow.
* **E3 (triangulate.R)** In the single-engine modes (flat/reliability) the
  cross-lens unit statistic compared raw incomparable scales, so the diffuse
  lens could never drive the any-lens test. The three lenses are now
  null-standardised INSIDE the unit statistic only; per-lens tests keep
  their raw-statistic maxT (and its width advantage - the AVP anchor is
  untouched), and the combined engine is unchanged. NOTE for the MS:
  equation 5's second line still describes the earlier rank-p ACAT engine
  fusion; the code uses a continuous standardised-mean fusion because rank
  granularity floors family-wise p at ~family/(B+1).
* **E4 (dmsa_align)** A user polarity table now MERGES with the bundled
  curation, as the documentation always promised: user rows override
  gene-for-gene, the bundled table fills the rest, and a message states the
  arithmetic. It used to silently replace the whole table.
* **E5 (dmsa_lfdr)** The pi0 estimator read the mass between the sample's
  own quartiles - 0.5 by definition - so the EB arm always assumed ~no
  signal. Replaced with robust Efron central matching (iteratively trimmed
  null fit, mass counted in a fixed window). EB results are less
  conservative and now respond to the data; the permutation engines are
  untouched. A follow-up adversarial pass on the fix (six simulation fronts)
  found one genuine anticonservative hole - a HEAVY-TAILED null is misread
  as signal by any Gaussian two-group model (under a pure t5 null the
  selection arm called false units in 80-100% of runs) - and added a
  two-radius symmetric-tail-excess guard: both tails significantly heavier
  than the fitted Gaussian triggers a warning and the documented NULL
  fallback to the frequentist rule. With the guard: t5 null false-selection
  0% at 1000 units (16% of runs select ~0.3 units at 200 - documented
  residual), Gaussian null and one-sided signal unaffected (~3% fallback,
  power intact), symmetric two-sided signal falls back conservatively to
  the frequentist arm. The two-groups model's remaining assumptions
  (minority signal; |z| ~ 2 largely invisible - both CONSERVATIVE failure
  directions) are now stated in the help page.
* **E6 (dmsa_omnibus)** The hand-rolled block shuffle scrambled values
  across families whenever rows were not block-contiguous. It now uses the
  same size-stratified whole-block swap as every other engine; p-values are
  layout-invariant.
* **E7 (dmsa_moderate)** The Johnson-Neyman region was always reported as
  "significant OUTSIDE [l, u]"; when the quadratic opens upward (strong main
  effect, weak interaction) the true region is BETWEEN the roots and is now
  reported as such (plus "everywhere"/"nowhere" for the no-real-roots
  cases). All simple-slope/JN ingredients now come from ONE fitted model
  (the mixed model when it is the engine) instead of splicing lm and lmer
  estimates.
* **E8 (dmsa_mediate)** A custom cascade or cascade CSV path no longer
  silently receives the bundled Alpha direction map; it is refused with the
  two explicit options spelled out.
* **E9 (dmsa_align / dmsa_scores / dmsa_triangulate)** `dmsa_align()` on a
  probe vector now returns rows in the CALLER's order, and the consumers
  verify the positional join on canonical CpG ids: an unambiguous scramble
  of the same probes is repaired, an overlapping mismatch is a hard error,
  and unrelated naming conventions keep the positional contract. A silent
  multiplier-on-wrong-probe run is no longer possible.
* **E10** minP + the combined engine is refused in `dmsa_change()` (it was
  silently running maxT while recording "minP"); the minP branch uses one
  p-value convention on both sides with NA-honest denominators; `pc1` no
  longer returns silent all-NA on an all-zero-multiplier set; the
  `dmsa_model(nulls=)` documentation now describes the actual TRUE/FALSE
  contract; the RNG save/restore blocks now truly restore "no prior state";
  a dead helper was removed; Zhou-manifest `*_beg` coordinates (BED,
  0-based) are converted to the 1-based scale of everything drawn beside
  them; and the moderated permutation in `dmsa_outcome()` now residualises
  the product term on S, W AND the covariates (full Freedman-Lane), so the
  observed and permuted regressors share one geometry.

## Second pass, same day

### The pair path is live (spec 1/3): `direction_source`

* `dmsa_frame()` gains `direction_source = c("auto", "cpgdirection",
  "bundled")`, `reference = "alpha"`, `tissue = "blood"` and
  `replicate_policy = c("collapse_site", "keep_correlated")`. With
  cpgdirection installed (the default `"auto"` finds it), the map the engine
  consumes is built at run time from CpG x gene PAIRS:
  measurements -> `cpg_gene_pairs()` discovery -> reference filter ->
  pair ledger -> engine map. Without it, the bundled per-column snapshot is
  used and the fallback is RECORDED in the corrections table. A user-supplied
  `map` table forces the bundled source. The engine itself is unchanged:
  only where its map comes from is new.
* On the pair path the biological reference (`dmsa_reference`; `"alpha"` =
  the 2026c codebook) defines membership; the reference never selects CpGs.
  One CpG mapping to several genes enters as several pairs - including with
  opposite directions - and technical replicate designs of one CpG are
  averaged into one measurement under `collapse_site` (a missing value in a
  replicate keeps the site visible to `missing_methylation`).
* p_plus for each pair is resolved by `dmsa_align()`'s own documented
  evidence ladder, applied once at map construction - not re-implemented.
* `cpg_map` semantics are now pair-level (spec 13/14): "full" = every usable
  pair; "confidence" = the strict subset resting on tier M (measured) or S1
  (top SMR) evidence.

### Coverage is itemised (spec 15/16/17)

* New `dmsa_coverage(frame, level = c("system", "gene", "pairs"))`: reference
  genes vs testable genes per system, the gene-level correction family size,
  per-gene pair counts with abstention reasons, and the full pair ledger.
* Every report on the pair path writes `tables/cpg_gene_pair_ledger.csv` -
  every discovered pair, used or not, with its reason. A gene that lost all
  of its evidence is visible there even when nothing else mentions it.
* The frame stores `$pair_ledger` and `$measurements`; `print()` names the
  direction source.
* Spec 11: the legacy map builder no longer deduplicates by column (a no-op
  for the bundled one-row-per-column snapshot; a user map may legitimately
  carry pairs).

## The 0.99.9 foundation

### The biological reference and CpG mapping (the 2026c re-architecture)

* **`alpha_reference()` now carries the full 2026c codebook: 1,282 genes in
  187 modules across 30 systems** (previously 549 genes - only those with a
  retained Alpha probe). Biological membership (system > module > gene) is now
  independent of measurement availability; a gene no longer vanishes from its
  system because no probe was retained for it. Counts are derived from the
  bundled table, never asserted. `alpha_reference_check()` verifies the
  invariants.
* New CpG-to-gene mapping layer (`R/cpg_mapping.R`): measurement columns are
  canonicalised, technical replicate designs of one CpG collapse under
  `replicate_policy = "collapse_site"` (under this policy the measurement id
  IS the canonical CpG id), one CpG legitimately maps to MANY genes (pairs are
  preserved, never deduplicated), and direction comes from
  `cpgdirection::cpg_gene_pairs()` resolved at call time - one row per
  CpG x gene pair, pair-local direction, pair-specific abstention. A pair
  ledger records why every pair was or was not used; cpgdirection's own
  `pair_id` is preserved as `cpgd_pair_id` beside dmsa's minted key.
* `cpg_map` default flips **"confidence" -> "full"**: every pair with a usable
  direction call is analysed; "confidence" remains as an opt-in sensitivity.

### Adjusted p-values are NOT comparable to earlier versions

* The analysed universe (all user CpGs x discovered targets, rather than
  Alpha-retained probes) and the `cpg_map` default flip both change the SIZE
  of the level-local families that maxT/minP correct within. Adjusted
  p-values therefore move between versions even where a gene's data did not
  change. The correction machinery itself (maxT/minP, level-local families)
  is unchanged.

### Safety (hard errors where silence used to be)

* `system = "oxytocin"` (a name passed to a TRUE/FALSE level switch) is now a
  hard error pointing to `systems =`; all level and control switches validate
  strictly; `systems = "all"` is explicit.
* A requested exchangeability block that is missing or all-singleton is a hard
  error even under `autofix = TRUE` (it used to silently become unrestricted
  permutation).
* Categorical covariates stay categorical: an explicit factor is honoured even
  with numeric-looking levels; word-level factors no longer dissolve into NA.
  Logical columns become 0/1 with a note (a logical OUTCOME used to dissolve
  into all-NA).
* Sample alignment (`sample_id`, `sample_align`): rows are matched by
  identifier wherever identifiers exist - including NUMERIC id columns -
  with reordering announced; a matrix whose rownames match no column of
  `data` is an error, not a silent positional guess.
* Missing methylation needs a declared policy (`missing_methylation =
  "error"` (default) / `"drop_probes"` / `"common_complete_rows"`); a probe
  can no longer be silently deleted by one missing value. The methylation
  scale is declared (`methylation_scale = "auto"/"beta"/"M"`); declared-beta
  values outside [0, 1] are an error, never clamped.
* The outcome family is per outcome (`outcome_type` accepts a named vector);
  a binary second outcome no longer flips the family of a continuous first
  one. Multinomial is refused with the outcome named.
* `frame_role` is validated for consistency with the shape scan and the
  resolved model orientation of every path is stored and printed.
* `dmsa_report(overwrite = FALSE)` refuses to overwrite an existing report -
  including a summary-only one.
* The system level forwards the chip random intercept exactly as the finer
  levels always did.

### Report

* Spec 27: per-lens survivor flags (`selected_coherence`, `selected_composite`,
  `selected_diffuse`, `n_lenses_hit`, `any_lens_hit`) and the joint-statistic
  `p_unit` / `p_unit_adj` are carried in the units table. The headline
  `selected` flag remains any-lens survival.
* Many small report fixes: multi-outcome level labels, PDF locus panels,
  NA-safe shape printing, the run's alpha used in print, figures recorded
  only when actually written.

# dmsa 0.99.8

* Bioconductor pre-submission polish; no behaviour changes. The R dependency
  rises to 4.6.0 (per BiocCheck), `seq_len()` replaces literal `1:n` in four
  spots, conditional message fragments move out of the signal calls, and the
  pipelines vignette gains chunk labels and a session-info section.

# dmsa 0.99.7

* **New: `dmsa_import()` - one door in from every preprocessing pipeline.**
  Takes the end object of minfi (any `SummarizedExperiment` - a
  `GenomicRatioSet`, `MethylSet` and relatives are recognised structurally,
  no minfi dependency), sesame (`openSesame()` beta matrix), meffil,
  ewastools and ChAMP (beta/M matrices in either orientation), RnBeads and
  methylprep (exports), or a `.csv`/`.rds` path, and returns a
  `dmsa_import` object holding a samples-by-probes beta matrix plus an
  aligned phenotype sheet. `dmsa_frame()` accepts it directly as its first
  argument, so the whole bridge is
  `dmsa_frame(dmsa_import(x, pheno), outcome = ...)`.

  Nothing is guessed silently: probe-axis orientation is read from the
  dimension names or must be declared; M-values are detected by range and
  converted to betas with a note (`dmsa_frame()` converts betas to M-values
  itself, so arriving in M twice would corrupt the frame); and the sample
  sheet is matched to the array names by identifier - literal columns,
  minfi `Basename` paths, or the `Sentrix_ID`/`Sentrix_Position` pair -
  with ambiguity and failure as errors that name what was tried.
  Raw containers (an `RGChannelSet`, a list of sesame `SigDF`s) are refused
  with the preprocessing recipe rather than half-handled.

* **New vignette - "One door in: from any preprocessing pipeline to DMSA".**
  The adapter's rules and per-pipeline recipes, plus a real-data invariance
  check: 94 EPIC arrays - a different physical scan of the canonical
  participants, sharing zero barcodes with the bundled maps' source - through
  minfi noob, minfi raw, meffil, sesame and ewastools, then the same two
  pre-registered questions asked of every result. The positive control
  (immune, inflammation & HLA x immune-cell fraction) hit the permutation
  floor (p = 0.0005 = 1/(B+1), B = 1999) in **all five pipelines**, raw
  unnormalized betas included, with all 8 modules selected everywhere and 17
  genes significant in all five (Jaccard 0.81, gene directions identical for
  92% of shared genes). The theory null (HPA x BSI) stayed null in all five
  (p = 0.41-0.81) with the cascade never descending. The pipeline choice
  moved probe counts (sesame's default masking: 26 usable vs 34-35), never
  the answer.

* **Fixed: a methylation matrix keyed by bare probe ids failed downstream.**
  `dmsa_frame(methylation = M)` with probe-id column names was accepted at
  intake (the map keys on `column` OR `probe`), but the extracted block kept
  the probe ids while every downstream consumer - the per-set column lists,
  the reporters, the aligner - keys on the map's `column` names, so the
  B = 49 pilot indexed out of bounds and the frame could not run. The
  extracted block is now normalised to the canonical column names at
  extraction. (Alpha builds carry methylation inside `data` under `column`
  names and never hit this path.)

# dmsa 0.99.6

* **Permutation groupings now accept a column name - and the old failure mode
  is refused, not silently endured.** `block`, `ri_group` (and the new `id`)
  in `dmsa_triangulate()`, and `block` in `dmsa_model()`, `dmsa_gate()`,
  `dmsa_gate2()`, `dmsa_change()` and `dmsa_block_index()`, were documented as
  per-row label vectors. Passing a column NAME - the natural thing to do - did
  not error: a length-1 string recycled into ONE block holding every row, the
  only within-block permutation was the identity, every permuted statistic
  equalled the observed one, and **p = 1, always, silently**. Total power
  loss. (`dmsa_frame()` was unaffected - it resolves its block columns
  internally, which is why published frame-based results are sound.)

  A column name (or several, combined by interaction, matching
  `dmsa_frame()`'s semantics) now resolves against `data`; a vector is
  validated against `nrow(data)`; and the degenerate all-rows-in-one-block
  case errors with the reason instead of returning p = 1.

* **The design now declares the data shape, and the declaration is
  verified.** The engines cannot tell wide from long by looking - two rows
  per person and two people look identical without an identifier - and the
  cost of guessing wrong is not cosmetic: with two-wave data (ICC ~ .6)
  treated as independent rows, the measured true-null rejection rate is
  **.175 at nominal .05**, roughly 3.5x.

  `dmsa_design()` gains `id` (participant column) and `format`
  (`"wide"`/`"long"`). Long requires an id - rows cannot be kept together
  without knowing whose they are; wide does not, and is verified when an id
  exists, taken on trust (with a note) when not. `dmsa_check_design()`
  verifies the declaration against the rows: declared wide with repeating
  ids, declared long with none, and a repeating id that sits in neither
  `random` nor `exchangeable` are all named problems.

  `dmsa_triangulate()` gains `id =`: when any participant contributes more
  than one row, a missing `block` - or a block that lets an id span two
  blocks - is an immediate error instead of a silent 3.5x inflation.

  `alpha_design()` declares the shape for all four Alpha builds (1: wide,
  id = "ID"; 2 and 3: long, id = "ID"; 4: wide - the child build has no
  per-person column, so its wide declaration is taken on trust).

* Regression tests: name-vs-vector equivalence, the degenerate-block error,
  the id guards, and the full shape-declaration matrix.

# dmsa 0.99.5

* **`dmsa_align()` now refuses a probe whose confidence is missing, instead of
  losing it silently or assuming it certain.** A
  `cpgdirection::cpg_expression_direction(tissue = "all")` result fills
  `best_confidence` only for its catalogue layers. SMR rows arrived without
  one, were given an `NA` weight, passed the usability check, and then
  disappeared inside `dmsa_triangulate()` - no error, no warning, nothing in
  the output. On one Alpha system that was 940 of 1,548 probes.

  P(d = +1) is now recovered in order from `p_plus`, `probability_plus1`,
  `best_confidence` / `confidence`, and then - for SMR rows - **the published
  accuracy of the row's own SMR tier** (S1 0.96, S2 0.85, S3 0.705), never the
  universal/distance probability, which would substitute a 0.60-0.65 prior for
  a 0.95-0.97 validated causal call.

  A table with no confidence *column* at all is still treated as certain, which
  is the documented behaviour for a plain `d`-only table. A table that has
  confidence columns but leaves a row empty now marks that row **unusable**
  with reason `no_confidence` and warns - filling it with certainty would give
  the least-supported rows the most weight.

* **`dmsa_align()` warns when a direction layer is one-way by construction.**
  cpgdirection's distance curves require unanimity across three tissues and the
  blood curve never exceeds 0.449, so `distance_only` cannot return `+1` for
  any CpG at any distance: every direction it contributes is `-1`. On the Alpha
  panel that is 6,297 calls, 0 of them positive. Aligning a set on a block of
  probes that all point one way by construction manufactures coherence that is
  not in the data, so `distance_only` and `targeted_last_resort` are now named
  in a warning. A generic backstop also fires when every call in a set of 50 or
  more shares one sign, whatever the source.

* Four regression tests: the recovery order including SMR tiers, the
  column-versus-value distinction, the one-way guard, and the invariant that
  "usable" and "reaches the test" are the same set of probes.

# dmsa 0.99.4

* **Critical fix: `alpha_polarity()` read the wrong table.** It returned a
  115-gene, six-system draft while `dmsa_polarity()` returned the audited
  1,234-gene 2026c table covering all 30 systems. Because a gene with no
  `w_g` entry stops an alignment by design, `dmsa_align(level = "system",
  polarity = "alpha")` **failed on most of the panel** - including the
  oxytocin/vasopressin system, where 34 calcium-handling genes had no entry.
  System level is the level this method is meant to be used at, so this
  disabled the primary use case rather than degrading it.

  `alpha_polarity()` now reads `alpha_polarity_2026c.csv`, falling back to the
  legacy draft only if the audited file is absent. Its `"status"` attribute is
  `"audited"` accordingly, `system_id` is coerced to character so a numeric
  `system_id =` matches, and an empty polarity subset now errors with the id
  rather than proceeding. With the fix, all 30 systems align: a median of 57
  usable probes each, 23 systems with 30 or more.

* `alpha_polarity()` gained documentation of what the two bundled tables are
  and why the audited one matters, plus regression tests covering panel-wide
  coverage and a system outside the legacy six.

# dmsa 0.99.3

* **A second direction layer ships: nasal epithelium.** `dmsa_align()` and
  `dmsa_directions()` gain `tissue = c("blood", "epithelium")`. Direction is
  tissue-specific, and the studies that will use this package are not one
  tissue: adult saliva by passive drool runs about 82% immune and takes blood
  calls, while neonatal and infant buccal swabs are .89 to 1.00 epithelial and
  take epithelial ones.

  Shipping blood alone would have made the wrong choice silent. Of the 126,112
  probe-gene pairs both layers call, 31% disagree in sign - and the
  disagreement is not a confidence artefact: among pairs both layers call at
  tier A it is 32%. The layers are also complementary rather than nested, each
  calling tens of thousands of probes the other declines to call, so the wrong
  layer does not merely mis-sign a sample - it can have nothing at all to say
  about the genes in the set. A buccal sample scored against blood calls would
  have reported full coverage while drawing every direction from the wrong
  tissue.

  The epithelial layer is 688,843 pairs over 277,830 probes and 4,447 genes,
  in 1.6 MB, on the same schema and versioning as blood. Which layer produced
  an alignment continues to travel with the result in the `"direction_map"`
  attribute.

# dmsa 0.99.2

* **Direction calls now ship with the package.** `dmsa_align()` accepts a plain
  character vector of CpG identifiers and takes both the direction calls and
  the probe-to-gene mapping from a bundled map, so a user who has just finished
  normalising can align without obtaining a direction resource first. That
  round trip was the largest barrier between a normalised matrix and a DMSA
  result, and it is now gone for the common case.

  The bundled layer is blood, restricted to probe-gene pairs carrying an actual
  call: 488,204 pairs over 237,761 probes and 5,376 genes, in 1.3 MB. About
  seven in ten rows of the full table abstain, and an abstention cannot be
  aligned to, which is why the reduction costs nothing DMSA could have used.
  Other tissues, the brain bridge, the SMR layer and the full 3.1M-pair table
  remain in `cpgdirection` (\doi{10.5281/zenodo.22024185}), which stays the
  complete resource. Supplying your own table works exactly as before.

* **`dmsa_directions()`** exposes the bundled map directly, for checking what
  fraction of a probe set carries a call before committing to an analysis.

* **Every alignment records which direction map produced it**, in the
  `"direction_map"` attribute. Annotation drift silently changes enrichment
  results, and a direction map is an annotation; a DMSA finding that moved
  because the map moved would be the failure this package exists to name. The
  map is versioned and the version travels with the result.

# dmsa 0.99.1

* **`dmsa_pdance()` is now part of the package.** The p-dance test - does a
  finding's p-value move when the analyst changes which *other* probes are in
  the analysis? - existed only as an analysis script behind the manuscript's
  Figure 3, while the README advertised it as a function. It is now exported,
  documented and tested, in both forms: `"dropout"` removes a random fraction
  of the non-focal probes, `"addition"` adds probes from a reservoir outside
  the set. The two move competitive and whole-set engines in opposite
  directions, which is why running only one understates the problem.

  It is engine-agnostic by construction: an engine is any function from a set
  of probes to a p-value, so a competitor is tested by exactly the procedure
  that tests DMSA. The focal unit's own probes are never perturbed, so any
  movement is attributable to the analysis set rather than to the evidence.

# dmsa 0.99.0

Bioconductor submission version. No change to any estimate, permutation, or
correction; every number in the manuscript reproduces from v1.20.2, archived at
[10.5281/zenodo.22042907](https://doi.org/10.5281/zenodo.22042907).

Bioconductor requires new packages to enter review at `0.99.z` and to release at
`1.0.0`, so the version number resets here. The GitHub and Zenodo lineage
continues at `1.20.x`; the concept DOI
[10.5281/zenodo.22023957](https://doi.org/10.5281/zenodo.22023957) resolves to
whichever is current.

* **`cpgdirection` is no longer a declared dependency.** Bioconductor cannot
  install packages from GitHub, so it forbids declaring them even in `Suggests`.
  The dependency was already optional in fact — every call site guarded it and
  either named an alternative or degraded — so it is now resolved by name at
  call time through a documented `.cpgd()` helper. Behaviour is unchanged when
  the package is installed, and the error messages now name the repository and
  the DOI when it is not.
* **`biocViews` added**, and `VignetteBuilder: knitr`.
* **A vignette**, `vignette("dmsa")`: a worked analysis in which a
  direction-blind test reads null and DMSA recovers the signal from the same
  simulated data, with the cancellation constructed explicitly.
* `NEWS.md` now ships inside the package rather than being excluded from the
  build, and a stray log file was removed from the package root.

# dmsa 1.20.1

Test-suite fixes only. No change to any result.

* `test-frame.R` asserted the literal word "binary" in a design note. Rewording
  that note in 1.19.0 - necessary, because "switched to logistic" wrongly
  implied a logistic model is fitted - therefore failed a test whose behaviour
  had not changed. The test now asserts the claim (a note exists for the
  outcome, it identifies a two-level outcome, and it describes a contrast)
  rather than the prose.
* Three legitimate warnings in the reference tests were surfacing as uncaptured
  warnings in the suite. They are now captured with `expect_warning()`, so a
  clean run reports no warnings.

# dmsa 1.20.0

Three review items, all verified on real data.

* **A binary outcome now gets its response curve.** The shape scan used to skip
  a non-gaussian outcome entirely. The reasoning was half right: putting a 0/1
  outcome on the left of a least-squares fit is a linear probability model, and
  the Sasabuchi and Fieller statistics are undefined on it. But a binary outcome
  IS non-linear - probability is bounded at 0 and 1 - and on the LOGIT scale the
  quadratic term is perfectly well defined. Curvature is now tested by
  likelihood-ratio against the linear logistic model and the fitted response
  curve is drawn on the probability scale, with observed proportions in
  equal-count bins carrying Wilson intervals, the data as a rug, and the same
  solid-inside-90% / dashed-beyond extrapolation rule as the other shape
  figures. On the Alpha parent build this surfaced curvature at **p = .0062**
  that the old skip reported as nothing at all.

* **A duplicate unit no longer survives on an uncorrected family of one.** A
  gene with a single direction-called probe creates a PROBE family of one, where
  maxT has nothing to correct against. The identical measurement therefore came
  out at the raw p (.028) while its gene-level twin, corrected across a 6-gene
  family, failed (.108) - and because the probe was flagged a duplicate it was
  also excluded from the figures, so the run reported a survivor it never drew.
  Duplicates now inherit the adjusted p of the unit they duplicate, which is the
  multiplicity actually incurred; each unit's own value is kept in a
  `*_adj_own` column so the inheritance is auditable.

* **Legend collision in the moderation panels.** The legend was drawn at
  "topleft" over a full-height y range and sat on the scatter and the slope
  lines. Headroom is now reserved and the legend has an opaque backing.

* The overview figure title used the raw outcome column name rather than the
  supplied label, and wrapped by character count rather than measured width.

# dmsa 1.19.0

Five review comments from the Set 1 rev-3 output, all implemented.

* **Categorical outcomes were described as if they were continuous.** One
  sentence template said "Higher X was associated with ...", which produced
  **"Higher sex was associated with ..."** and **"Higher time was associated
  with ..."** in shipped output. Three cases are now told apart, and by the
  DATA rather than the column name - a two-level outcome that varies *within*
  the finest permutation block of a before/after build is a within-person
  change, one that does not is a group contrast:
  - continuous: "Higher X was associated with ..." (unchanged)
  - binary: "Relative to female, male showed methylation consistent with ..."
  - wave: "Within person, from T1 to T4, the methylation of MEG3 shifted
    toward a LOWER expression tone."
  New `outcome_levels = c("T1", "T4")` names the two levels; without it they
  are named by their coded values, which is ugly but never ambiguous.
* The design note `binary outcome under gaussian -> outcome_type switched to
  'logistic'` was misleading - methylation is the response and the outcome is a
  predictor of it, so **no logistic model is fitted**. It now says what actually
  happens.
* **Labels reached exactly one place in the summary.** With
  `outcome_label` set, the Results paragraph honoured it and the report header,
  the survivor list, the system line and the Methods paragraph all still printed
  the raw column name. All of them now use the label. Tables gain a
  `unit_label` column and keep `unit` as the join key; design notes keep raw
  names, because they refer to columns in the data file.
* **Figure text was cut off.** A module-level three-way title measured **15.25
  in of text into 6.6 in of usable width - a 2.31x overflow**, centred, so it
  was cut at both ends. Three separate causes: the p-values were concatenated
  into the title (they are now a subtitle); wrapping used a fixed CHARACTER
  count rather than measured width; and `strwidth()` measures against the
  panel's `par("cex")` while `mtext(outer = TRUE)` draws at device scale, which
  under-measured every outer title by 1.2x in a 2x2 layout. Titles are now
  fitted before the layout is set, the outer margin is sized to the number of
  lines they need, and long labels shrink rather than clip. Also fixed a
  non-terminating loop in the truncation fallback.
* **A single CpG is now named `GENE's cgXXXXXXX`** in figure titles, figure
  file names, axis labels, the summary prose and the tables. A bare probe id
  does not say which gene's level-local family it was corrected inside, and the
  family IS the gene.
* **Moderated non-linear models gain a surface figure** - 3D perspective plus a
  filled contour with the observed data drawn on it. Only for a moderated
  non-linear model: a plain shape scan has one predictor and therefore no second
  axis. Grid cells with no observation within one bandwidth are drawn grey (3D)
  or blank (contour), and the figure prints what fraction of the plotted surface
  is extrapolation - a quadratic-by-moderator surface extrapolates violently at
  the corners, and a 3D surface makes that look like data.

# dmsa 1.18.0

**Chip now enters as a random intercept, as the Project Alpha covariate
contract specifies.**

`tutor/covariate_sets.csv` states `1 parent T1, chip_T1, random, (1 | chip_T1)`.
Through 1.17.0 dmsa entered chip as a FIXED factor - 61 levels, 60 df on the
parent build, with 7 singleton chips fitted perfectly. That was an undocumented
deviation from the pre-specified contract.

Why chip is needed at all (parent build, n = 396):

* chip is **not** associated with any outcome (F ~ 1, p = .25-.73 across five
  outcomes), so it is a precision variable, not a confounder - omitting it
  cannot bias an effect estimate;
* but it is a large technical-variance component in methylation: **61%** of
  probes show chip signal beyond the control-probe SVs at p < .05 (5% expected)
  and **30%** at p < .001 (0.1% expected);
* the control-probe SVs do not substitute for it (ctrlSV3 R2 = .32,
  ctrlSV5 R2 = .36 against chip);
* families **span** 2-3 chips, so chip is crossed with the cID permutation
  block and blocking on cID does not absorb it.

Implementation: for a one-way random intercept the GLS solution is exactly OLS
on quasi-demeaned data (Balestra-Nerlove), with the variance ratio estimated by
REML. No mixed-model package is required - the REML step is a 1-D optimisation
in `R/ranef.R`. Validated against `lme4::lmer(REML = TRUE)` on the Alpha parent
build: **coefficients and t-statistics agree to 5e-13**, variance components to
4e-06, including the sigma^2_u = 0 boundary.

Because the transform is a linear map applied once, before any permutation,
Freedman-Lane block permutation is unaffected. Validity does not depend on the
variance ratio being exactly right: the null is generated under the same
transform, so the variance ratio affects power, not type-I error.

What it costs and what it buys: the random intercept spends a median of **~16
effective df** instead of 60, and **0** on the ~20% of probes where the chip
variance is estimated as nil. On the two deposit genes it is appropriately more
conservative than the fixed factor, and both still survive:

| gene | random (contract) | fixed (1.17.0) | chip omitted |
|---|---|---|---|
| PLAGL1 | adj .0260 survives | adj .0100 survives | adj .1360 |
| AVP | adj .0140 survives | adj .0080 survives | adj .0140 survives |

Omitting chip loses PLAGL1 altogether.

* New argument `chip_effect = c("random", "fixed", "none")`, default `"random"`.
  `"fixed"` reproduces <= 1.17.0 behaviour and is ~5x faster, which is useful
  for structural test runs; the Design notes now record which was used and flag
  `"fixed"` as a deviation from the contract.
* Fixed a silent bug introduced with this change: under `"random"` the chip
  column was dropped from the retained data, which would have left the grouping
  vector NULL and the random intercept never applied.

# dmsa 1.17.0

Five defects found by auditing a 29-model battery run on all four Project Alpha
builds. All five were silent: every one produced a report that compiled cleanly
and said something false or incomplete.

* **Surviving modules were never named in the report.** `hits.csv` carried them,
  the summary and the "How to report this" Results block did not. In one run
  four modules survived, the strongest at adjusted p = 0.002 - better than two
  of the three genes that *were* reported - and a reader of `summary.md` alone
  would never have known. Module survivors are now listed in the summary and
  get a Results paragraph, with the caveat that a module pools several genes.
  `"nothing survived"` is also no longer printed when a module did.

* **The significance star at system level marked the wrong quantity.** The star
  came from `selected`, which is driven by the best lens's *adjusted* p, but it
  was printed beside the ACAT *omnibus* p. That produced `omnibus 0.0658 *` and
  `omnibus 0.0820 *` - a significance mark on a p above alpha. The line now
  reads `omnibus 0.0820 (best diffuse adj 0.0220 *)`, so each star sits on the
  number it refers to.

* **The main-effect Methods block never stated the covariates.** The moderated
  block always did. A reader of a linear report could not tell what had been
  adjusted for - including in runs where a mistyped covariate was silently
  dropped, or where every covariate passed was unknown and the model was
  effectively unadjusted. The Methods text now names the covariates, the block
  structure, and warns explicitly when requested covariates were dropped.

* **"Coverage: 0 gene(s) ... untestable" in three of the four builds.** The
  untestable-gene count re-read the direction map from disk, so it still held
  parent-T1 column names and matched nothing in the child and before-after
  builds. Those reports claimed complete coverage while the parent build, on
  the same two systems and the same 18 genes, correctly reported 36 untestable.
  The raw map is now translated to the build before the intersection.

* **The console under-reported the table count by one** whenever `probes.csv`
  was written; it was never added to the tally.

# dmsa 1.16.0 (2026-08-19)

Everything in this release came out of running the 25-model Set 1 battery on all
four Project Alpha 2026d builds. Five of the twenty-five models failed and two
more passed while testing something other than what they claimed to.

## Three of the four builds could not be analysed at all (critical)

* **The direction map names probes in the parent-T1 form only.** The same CpG is
  `cg02187522_TC21_T1_AVP` in the parent build, `cg02187522_TC21_AVP` in the two
  before/after builds, and `cg02187522_TC21_T4_AVP` (plus `_maternal_T1_` and
  `_paternal_T1_` sets) in the child build. `dmsa_frame()` intersected the map's
  column names with the data's and found **zero** matches for the child, T1-T4
  and Oct-7 builds, stopping with "no mapped methylation columns found". The
  selection cascade lists every variant explicitly (`col_parent_T1`, `col_long`,
  `col_child_T4`, `col_child_maternal_T1`, `col_child_paternal_T1`), so the frame
  now translates the map to whichever form the data uses, prefers the variant
  belonging to the detected build when several match (the child workbook carries
  three), and records the translation in the corrections table. Passing
  `methylation` as an explicit vector of column names overrides the choice - that
  is how to analyse the maternal or paternal array inside the child workbook.

## Typing the contract out by hand gave a different model (serious)

* **`sex_c` was derived only under `covariates = "contract"`.** A caller who
  passed the same six columns explicitly got `sex_c: covariate absent from data
  -> dropped`, so the model was silently **not adjusted for sex** and returned
  different p-values from the identical-looking contract call. `sex_c` is a
  derived column, not a column of any workbook, and is now built whenever it is
  asked for by either route.
* **A missing covariate now raises a `warning()`**, not just a line in the
  corrections table. A mistyped covariate name is indistinguishable from an
  intended one, and dropping it quietly produces an unadjusted model that still
  reads as adjusted in the write-up.

## Moderated curvature was reported in prose and never drawn

* `moderation$selected` was keyed on the linear product's adjusted p alone. A run
  whose only survivors were **moderated curvature** (S^2 x moderator) named nine
  of them in the summary and produced **no figure at all**, and any user
  filtering `moderation$selected` saw none of them. The table now carries
  `curv_selected` and `any_selected`, and the figure loop draws both arms,
  skipping units flagged as repeating the same data at another level.
* The figure title now leads with the arm that actually carried the unit. A panel
  drawn for a curvature survivor used to be headed "family-adjusted p = 1.0000".
* `selected` is now guarded with `is.finite()`; an NA adjusted p produced NA rows.

## `chip` is an argument

* `dmsa_frame(chip = FALSE)` leaves the array batch factor out; a character
  string names the column to use; `TRUE` (default) keeps the previous
  auto-detection. On the parent build `chip_T1` has 61 levels across ~400 arrays
  with 7 singleton chips, so entering it spends 60 degrees of freedom on a
  nuisance factor that the control-probe surrogate variables already partly
  absorb. That is a judgement to be made and reported, so it had to be sayable.
  The before/after contracts forbid chip and still override the argument.

## Shared probes: audited, and the hazard is structurally absent

* 229 probes in the Alpha panel sit inside two panel genes and appear as a column
  under each name (IGF2+INS 31, GABRA5+GABRB3 25, HLA-DPA1+HLA-DPB1 23,
  PEG3+ZIM2 21, and others). Two such genes in one level-local family are not two
  independent units. `.frame_shared_probes()` now reads the cascade's explicit
  `shared_with` column rather than trying to infer sharing from duplicated probe
  names - which cannot work, because an Alpha column name embeds its gene, so one
  probe under two genes is two differently named columns.
* Auditing all 494 documented (probe, gene, partner) rows against the direction
  map: **not one shared probe survives into the map under both names** (62 under
  exactly one, 432 under neither). No DMSA analysis of this panel can
  double-count a shared probe. The check is retained as a guard for future maps
  and reports into `frame$shared_probes`.

## Two failures that gave uninformative errors

* A selected system left with no direction-called probe reached the aligner with
  a zero-row table and died inside it on `replacement has 1 row, data has 0`,
  which says nothing about what went wrong. It now stops with the system name,
  how many probes remain overall, and which systems they cover.
* When **every** mapped probe column has a missing value on the complete-case
  rows the frame now says so, and names the usual cause: the methylation set and
  the covariates were measured on different people. That is exactly the case in
  the Project Alpha child workbook, which carries the mother's and father's T1
  arrays but no maternal or paternal cell-composition or control-probe
  covariates - so those arrays cannot be analysed with the child T4 contract,
  and the frame now says that instead of failing obscurely.

## Smaller

* Moderation figure captions are wrapped to their panel (at most two lines).
  The four-line model statement on the left collided with the legend; a
  single unwrapped line then ran into the neighbouring panel's caption in the
  three-way layout and the Johnson-Neyman note ran off the panel edge. The full
  model statement belongs in the summary.
* Captions name the moderator by its display label throughout; one line still
  used the raw column name.
* The shape-scan summary said figures are drawn "for units whose adjusted
  quadratic p clears alpha". They are drawn for any surviving departure arm -
  quadratic, threshold or exponential. The sentence now says what the code does.
* The "no figures" message for a moderation run names both arms when
  `type != "linear"`.
* `dmsa_frame(chip = "pool")`, see above.
* "1 direction-called probes" now agrees in number.

# dmsa 1.15.0 (2026-08-18)

## `covariates = "contract"` is build-aware (important)

* **The contract was hard-coded to the parent T1 build.** Project Alpha ships
  four builds, each documenting a different covariate set in its own
  `Covariates` sheet. `.frame_contract()` returned the parent T1 set for all of
  them, so on the other three every contract column except `sex` was absent,
  silently dropped as "covariate absent from data", and the model ran on **sex
  alone**. On the two before/after builds that omits **`time`** - the design
  itself. The contract is now selected from the columns present:
  parent T1, child T4, or before-after. When no documented contract matches, it
  **errors** instead of degrading.

* **Repeated-measures builds get the right blocks.** The before/after builds
  list ID and cID as random effects. Permuting within cID alone would let two
  timepoints swap between members of a twin pair, breaking the within-person
  pairing; when the caller leaves `random_effects` at its default, those builds
  now block on `c("ID", "cID")` and say so.

* **chip is not entered where the contract forbids it.** The before/after builds
  list chip as NEVER include; it is now omitted there rather than added as a
  fixed factor.

* **The build and its contract are printed** in the frame's corrections table,
  so the model's covariate set is auditable from the output.

## Overview figure fixes

* **"Bold unit = survives" was wrong in every overview figure.** `axis()` takes
  a single `font`, not a vector, so `ifelse(selected, 2, 1)` applied the first
  element to every label: all units bold when the top row survived, none
  otherwise. Drawn with `mtext()` now, which vectorises.
* **A two-line wrapped title was clipped** at the top of the device.

# dmsa 1.14.0 (2026-08-18)

## The summary is what people will read, so it now says what it means

* **A family of one is marked as uncorrected.** At probe level a gene with a
  single direction-called CpG forms a family of one, where maxT has nothing to
  correct against and the "family-adjusted" p equals the raw permutation p. That
  number was being reported alongside genuinely corrected ones - and because it
  is the smallest, it looked like the strongest result in the run. The summary,
  the "How to report this" sentence and `print()` now all say
  **FAMILY OF 1, so uncorrected**.

* **The same data at three level names is no longer counted three times.** A
  module whose only direction-called gene is X, that gene, and (if it has one
  probe) that probe produce identical statistics. They were listed as separate
  survivors, multiplying the apparent number of findings: "6 survive" was really
  4 distinct results, and "2 survive" was 1. Duplicates are now detected by
  identical outcome, b and t, labelled `= <coarser unit>`, and a distinct count
  is given.

* **The Methods paragraph covers the moderated non-linear model.** It described
  only the two-way product, so a reader copying it would not know S^2 had ever
  been fitted. It now states the second model and its tested term.

* **The curvature paragraph interprets its sign** - positive means the curve
  becomes more convex as the moderator rises, negative more concave - instead of
  listing coefficients with no reading.

* **Curvature panels no longer say "Johnson-Neyman: slope".** The no-boundary
  branch kept the slope wording over a curvature axis.

* **Labels reach the rest of the report**: the moderation and curvature
  paragraphs, panel subtitles and the J-N/curvature axis were still printing raw
  column names while the titles used the labels.

* Singular "1 direction-called CpG".

# dmsa 1.13.0 (2026-08-18)

## Outcome types, labels, and the non-linear moderation figure

* **A nominal outcome is now refused instead of silently scored (bug).** The
  main path regresses methylation ON the outcome and tests a single 1-df
  coefficient. With `outcome_type = "multinomial"` a k-level nominal variable
  was coerced to numeric and analysed as an ordered, equally-spaced score -
  "delivery mode 5 is five times delivery mode 1". `dmsa_frame()` now errors and
  points at `dmsa_outcome()` (which implements a rank-based scheme for
  multinomial and ordinal outcomes) or a binary recode. The gaussian-with-few-
  levels note now says the variable is being treated as ordered.

* **`outcome_type` no longer implies a link function that does not exist.** The
  frame printed "model: logistic linear" while fitting OLS. It now states that
  the outcome enters as a PREDICTOR of methylation, so a binary outcome is a
  group contrast on a continuous response and no link function is needed or
  used. That is correct for the main path - and the shape scan, which puts the
  outcome on the LEFT, now refuses a non-gaussian outcome rather than fitting a
  linear probability model under Sasabuchi/Fieller machinery.

* **A covariate derived from an outcome is excluded up front (bug).** With
  `outcome = "sex"` the contract built `sex_c` from `sex`, then found the design
  rank-deficient, then dropped it - three messages leaving the user to work out
  that the covariate WAS the outcome. It is now excluded before fitting, with
  that stated as the reason. (The final model was correct either way.)

* **The non-linear moderation figure draws the non-linear model.** With
  `type != "linear"` the panel drew straight simple slopes while the table
  reported moderated curvature. It now fits and draws the curved model, the
  right-hand panel shows CURVATURE as a function of the moderator (J-N is about
  a slope, and with curvature the slope depends on the score too), and the title
  carries the moderated-curvature adjusted p.

* **Three-way panels name their `mod2` level.** The row label lived only in the
  outer note, so each panel was ambiguous on its own.

* **Display labels.** `dmsa_frame()` gains `outcome_label`, `predictor_labels`,
  `mod_label` and `mod2_label` - one per column, in order, or a named vector.
  They appear in figures and in the written report in place of raw column names.
  Purely cosmetic: no label reaches the model, the permutation or the correction.

# dmsa 1.12.0 (2026-08-18)

## Moderated non-linear is now a real test, and the shape prose reaches the report

* **`type` was silently ignored under `moderation = TRUE` (bug).** The moderation
  branch of `dmsa_report()` returns before the shape scan runs, so asking for a
  moderated non-linear model produced a moderated LINEAR one with no warning.
  There is now an actual moderated non-linear test: whether the CURVATURE of the
  tone score depends on the moderator, tested as the moderator by squared-score
  product with S, S^2, the moderator and S x moderator all present for
  marginality, maxT-corrected within the same level-local family. On the Alpha
  HPA/SNS panel 6 of 106 units survive it.

* **"How to report this" now covers the shape scan.** It was main-effect only, so
  a run with `type = "non-linear"` handed the user a Methods paragraph
  describing a purely linear analysis - which is what they would then publish.
  There is now a Results-shape paragraph stating the incremental test and its
  adjusted p, a Methods-shape paragraph, and a reference list (Fieller 1943;
  Ganzach 1997; Lind & Mehlum 2010; Morris et al. 2023; Sasabuchi 1980;
  Simonsohn 2018). The linear Methods paragraph now says outright that the
  unit-level tests are linear in the aligned score.

* **The shape table ranked by the quadratic arm only (bug).** Units surviving on
  the threshold or exponential arm were pushed off the bottom of the truncated
  table while still getting a figure, so the table and the figures disagreed.
  Ranking is now by the best adjusted departure p across arms, and any surviving
  unit is always shown.

* **A family of one is marked as such.** At probe level a gene with a single
  direction-called probe forms a family of one, where maxT cannot correct
  anything and the "adjusted" p equals the raw p. The table now says so on the
  row instead of presenting an uncorrected number as corrected.

# dmsa 1.11.0 (2026-08-18)

## The scan now covers the whole declared hierarchy, and the exponential arm is real

* **Probe level added, so `system > module > gene > probe` is true.** `.rp_units()`
  had no probe branch, so the shape scan and the moderated path both stopped
  silently at the gene level while the declared design said otherwise. Only the
  main-effect path ever reached probes. There is now a probe branch, with
  family = the probes of one gene, matching the family the main-effect path
  corrects within. Both the shape scan and moderation use it. On the Alpha
  HPA/SNS panel the scan goes from 22 units to 53: 2 systems, 9 modules,
  11 genes, 31 probes.

* **Every departure arm is maxT-corrected, not just the quadratic.** An
  `type = "exponential"` run reported its own arm raw while the quadratic beside
  it in the same table was corrected. Each arm now gets Westfall-Young maxT
  inside its level-local family on its own permutation stream.

* **The figure gate no longer keys on the quadratic alone.** A run whose
  exponential arm survived drew no figure, because the gate asked only whether
  the quadratic had survived. It now fires on any surviving departure arm, and
  the exponential fit is drawn and carries its p-values in the legend.

* **The panel carries numbers; the summary carries the argument.** The figure
  header had grown to seven lines of prose that collided with the legend. It is
  now three short lines - unit, shape verdict with the Fieller interval, and the
  solid/dashed convention - with every arm's raw and adjusted p in a single
  legend. The reasoning about marginality, Sasabuchi and Fieller stays in
  `summary.md` where it can be read properly.

* **A caution on the exponential arm.** exp(score) is a monotone transform that
  correlates strongly with the linear score, so the incremental test has low
  power and a null there is weak evidence. It asks whether the association
  accelerates, not whether it is non-monotone - only the shape test speaks to
  that. Stated in the summary whenever the arm is present.

# dmsa 1.10.0 (2026-08-18)

## The shape scan runs at every level, and is corrected like everything else

* **Scanned at system, module and gene, not only system.** An aligned tone score
  exists at every level of the cascade, so there was never a reason to ask about
  functional form at one level only. `type = "non-linear"` now scans each
  declared level.

* **Westfall-Young maxT inside level-local families.** The scan previously
  reported raw p-values and was corrected nowhere - the one part of DMSA sitting
  outside the error control the rest of the package is built on. That is
  untenable once the scan covers dozens of units. The quadratic departure term
  is now maxT-corrected within its level-local family on a shared permutation
  stream, exactly like every other DMSA statistic. It changes conclusions: in the
  Alpha HPA/SNS panel, ADRA2A's quadratic term is p = .032 raw and p = .060
  adjusted, so uncorrected it reads as a finding and corrected it does not.

* **Fieller interval for the extremum.** Lind & Mehlum present containment of a
  confidence interval for the turning point within the observed range as the
  test equivalent to the slope-sign test, and warn that the delta method is
  severely biased in finite samples. The scan now reports a Fieller 95% interval
  and whether it is contained. An unbounded interval - which is common, and
  which the delta method would hide behind a tidy-looking standard error - is
  reported as unbounded, because it means the turning point is not determined.

* **Figures are drawn for units that survive their family**, not for every
  scanned unit, and the run says so when nothing survives.

* **Documented against the source literature.** The Sasabuchi hypotheses, the
  max-of-one-sided-p convention, the interior-interval suggestion and the
  Fieller recommendation were each checked against Lind & Mehlum rather than
  reconstructed from memory.

# dmsa 1.9.0 (2026-08-18)

## The shape scan now respects marginality, and will not call a curve a U

Reported shape numbers change in this release. Anything taken from a 1.8.9
shape table came from the wrong specification and should be re-run.

* **Departure terms are now fitted WITH the linear term (breaking).** The scan
  substituted the transformed score *for* the linear score, so the quadratic arm
  fitted `outcome ~ score^2 + covariates` with no linear term. Because the score
  is standardised, that imposes a symmetric curve whose turning point is pinned
  at the mean, and a significant coefficient there is not evidence of curvature
  (the marginality principle; Morris et al. 2023). Each departure term is now
  entered alongside the linear term, so the reported p is the 1-df incremental
  test of curvature over and above the linear component - equivalent to
  comparing the linear and quadratic models directly. The threshold arm is fixed
  the same way.

* **The omnibus no longer pools the linear arm.** ACAT over {linear, quadratic,
  threshold} answered "is there any association", which a significant linear
  term alone can drive and which therefore cannot support a claim about shape.
  The omnibus now combines the departure terms only and is reported as such,
  with the linear p in its own column.

* **A U-shape claim now requires a U-shape test.** A significant quadratic term
  is consistent with a relationship that rises throughout but with increasing
  steepness, so it does not establish a U or an inverted U (Lind & Mehlum 2010;
  Simonsohn 2018). The report adds the Sasabuchi slope-sign test - the fitted
  slope must take opposite signs at the two ends of the observed range, with the
  p-value taken as the LARGER of the two one-sided p-values - and declares a
  shape only when that test clears alpha *and* the turning point falls inside the
  observed range. Otherwise the summary says to report a curvilinear component
  and stop.

* **A trimmed shape test, because the extremes are the sparse points.** Lind &
  Mehlum evaluate the slope at the minimum and maximum of the observed range,
  which in a skewed score are exactly where the data thin out. The scan repeats
  the test on the central 90% and flags any shape that passes on the full range
  but not on the trimmed one, so a U resting on a handful of tail cases is
  visible as such.

* **The figure draws the models its p-values came from.** 1.8.9 drew a
  hierarchical quadratic curve while labelling it with the squared-term-only
  p-value. Every fit shown is now the same hierarchical model that produced the
  p beside it, the quadratic's turning point is marked when it lies inside the
  data, and fits stay solid across the central 90% and dashed beyond.

* **Two cautions written into the summary.** The Sasabuchi test uses parametric
  standard errors rather than the block permutation, so it is anticonservative
  under clustering; and curvature in a multivariable model can be induced by an
  omitted product term (Ganzach 1997), so a curvilinear result is worth
  re-checking with the relevant interaction present.

# dmsa 1.8.9 (2026-08-18)

## Bug fix: `type = "non-linear"` computed a result and threw it away

* **The shape scan was unreachable.** With `type = "non-linear"` or
  `"exponential"`, `dmsa_report()` ran an ACAT omnibus over linear, quadratic
  and threshold arms of each system's aligned score - and then discarded it. The
  rows were assembled into a local list that nothing ever read. Consequences: a
  non-linear run returned a results table **byte-identical** to the linear run
  (verified: max absolute difference in `p_omnibus` = 0), the returned object had
  no shape component, and `summary.md` never mentioned that a scan had happened.
  A user asking for a non-linear analysis got a linear one with no warning. The
  scan now reaches `summary.md`, `tables/shape.*`, `print()` and
  `report$shape`, with the per-arm p-values that `dmsa_model_omnibus()` was
  already computing and `.report_shape()` was dropping.

* **It says what it is, and what it is not.** Each arm is fitted as the ONLY
  score term, so the scan asks which functional form of the whole-system score
  tracks the outcome - it is **not** a test of whether curvature adds anything
  to a linear fit, and it does not re-specify any unit-level test. The summary
  states that the system, module, gene and probe p-values remain the linear
  ones, so the scan cannot be mistaken for a non-linear model of the units.

* **A figure per scanned system+outcome.** The aligned score against the
  covariate-adjusted outcome, with the linear, quadratic and threshold arms
  drawn over the points and each arm's own p-value in the legend. Fits are
  **solid across the central 90% of the score and dashed beyond it**, with the
  count of participants outside and a density strip along the axis: a quadratic's
  tails are set by whatever few observations sit out there, and curvature carried
  by a handful of cases should not be able to look like the finding.

# dmsa 1.8.8 (2026-08-18)

## A surviving moderation now gets a figure, and the model gets stated

* **Simple slopes with partial residuals, plus Johnson-Neyman.** 1.8.7 said a
  moderation run draws nothing. That was the wrong fix for the right
  observation: an interaction that survives correction is precisely the result
  a reader cannot picture from a coefficient. Each surviving unit now gets
  `moderation_<unit>_<outcome>.png` - the tone-score slope at three levels of
  the moderator drawn over the partial residuals it was actually fitted on, and
  a Johnson-Neyman panel showing over what range of the moderator the slope is
  distinguishable from zero. A three-way run draws the pair twice, with the
  second moderator held at -1 SD and +1 SD.

* **The J-N panel says how many participants are in the significant region.**
  A boundary near the top of a skewed moderator can be perfectly real and still
  describe almost nobody. The subtitle gives the count and the percentage, and a
  density strip shows where the participants actually sit.

* **The simple slopes are drawn at observed quantiles, not at +-1 SD.** The
  +-1 SD convention puts a line at an impossible moderator value whenever the
  moderator is bounded and skewed - an adversity count has a floor at zero, and
  mean minus one SD is below it. The panel now uses the 10th, 50th and 90th
  percentiles, falling back down a ladder of quantiles when ties in a
  zero-inflated moderator collapse two of them onto the same value.

* **The Methods paragraph now enumerates the model instead of gesturing at it.**
  "All lower-order terms present" is not auditable. It now names every term in
  the fitted model, lists the covariates, states that with several outcomes each
  is additionally adjusted for the others, and says plainly that the cluster
  variable enters as a permutation block rather than as a random intercept - no
  mixed model is fitted, so there is no variance component to report.

# dmsa 1.8.7 (2026-08-18)

## Moderation runs get the same reporting support as main-effect runs

* **`summary.md` now writes a "How to report this" section for mDMSA too.**
  1.8.6 added the section but returned early when there was no main-effect
  table, so the branch a user is *most* likely to need help wording - a
  slope-on-a-slope, reported from one collapsed score, under one lens - was the
  one branch that produced nothing. The moderated block states the three things
  the main-effect wording never has to: that the unit is collapsed into a single
  subject-level tone score BEFORE the interaction is fitted, that only the
  composite lens is reported because coherence and diffuse do not carry through
  a product, and what the coefficient means in words. With two moderators it
  says so correctly - the coefficient is the three-way term, so it reports how
  much the first moderator's moderation depends on the second, rather than
  mislabelling it as a change in the simple slope.

* **A moderation run with `plots = TRUE` says why it drew nothing.** There is no
  per-probe effect and no three-lens panel once the unit is a single score. The
  run now names that, and names the switch that gets the figures back.

# dmsa 1.8.6 (2026-08-18)

## Reporting the result, not just computing it

* **A dedicated module-level figure.** `modules_<outcome>.png` plots every
  module in the selected systems across all three lenses on one family-adjusted
  scale. Modules are **blocked and named by system, with a rule between blocks**,
  because that boundary is the panel's actual claim: a module is corrected
  against its own system's modules and no others. Each row carries the module's
  literature evidence tier and direction-called probe count in the right margin -
  Moderate rendered in red, so a result resting on a weaker module definition
  cannot be read without noticing. The module level was the one tier of the
  cascade with tables but no picture.

* **`summary.md` gains a "How to report this" section.** DMSA asks a reader to
  hold three unfamiliar things at once - which **lens** carried the finding,
  which **family** it was corrected inside, and a direction that is a claim
  about **expression tone** rather than about methylation level - and a user
  handed `adj .0120, coherence, d = -1` has no way to know how to write that
  up. The report now writes the sentences, with this run's own numbers in them:
  a Results paragraph per surviving unit, a Methods paragraph carrying the
  actual settings, and the direction stated in words. A run with no survivor
  gets the calibrated-null wording instead, which is the harder one to phrase.

## Fixes

* **Every table is also written as .docx**, and if `gt` is missing the run says
  so once, up front, instead of leaving you to notice the absent Word file when
  you go to paste it. csv to compute on, html to read, docx to paste into a
  manuscript - a results table retyped by hand is a table with typos in it.
* **The overview figure's title no longer runs off the right edge.** Outcome
  names are the user's own column names, and "Attachment_Anxiety_General_T1"
  already overran the device. The title wraps to the panel width.
* **A locus panel drawn without a gene model says so on the console**, naming
  the argument that fixes it. The caption always said it, but captions get
  skimmed, and a panel with no exons is half the point of the panel.

# dmsa 1.8.5 (2026-08-18)

* `dmsa_report()` now says where it wrote, every run, with the path resolved:
  `outdir` stays relative by default - the working directory is where every
  other R tool puts results, and changing that would be the more confusing
  choice - but "dmsa_output" is only findable if you already know what it is
  relative to. The message appears whether or not the result is printed, so
  assigning `rep <- dmsa_report(fr)` no longer ends a 35-minute run in silence.
  `print()` resolves the path too.

# dmsa 1.8.4 (2026-08-18)

* The progress bar advances through the probe level and the figure/table
  phases instead of sitting at 71% and jumping to 100% at the end. The two
  trailing phases were counted in the total but never ticked.

# dmsa 1.8.3 (2026-08-18)

## Two conveniences for long runs

A report is minutes of permutation with nothing on screen - the parent cohort
at B = 1999 across 19 systems takes 35 minutes - which leaves no way to tell
"working" from "hung".

* **A 0-100% progress bar.** `dmsa_frame(progress = )`, or override at run time
  with `dmsa_report(frame, progress = TRUE)`. It counts triangulation calls, the
  unit of work that actually costs time, so the bar is proportional rather than
  a decorative sweep.
* **A completion sound.** `dmsa_frame(beep = )`. `TRUE` uses beepr sound 8, a
  number picks another, `FALSE` is silent.

Both default to `interactive()`. A progress bar writes carriage returns, which
turn a piped log into soup, and a package that makes noise during `R CMD check`
or a testthat run is one people learn to mute - so both are live at the console
and off under `Rscript`, check and tests. beepr is a **Suggests**: without it
the sound is silently a no-op, never load-bearing.

# dmsa 1.8.2 (2026-08-18)

## Regression fix: maxT had lost its information advantage under a single engine

Caught by running the shipped package against the pre-registered proof of
concept and comparing every number to the deposited results.

A change on 15 August 2026 standardised each lens statistic by its own
permutation null before the family-wise step. That was introduced to fix a real
problem in the **combined** engine - a rank-p fusion floors the family-wise p
near (family size)/(B+1), which silently zeroes power for large families at
modest B - but it was applied to the **single-engine** path too, where it does
not belong.

Standardising divides out exactly what maxT's advantage is made of. A unit whose
probes all point the same way generates the **widest** permutation null in its
family and is therefore charged the least; rescale every unit to unit spread and
that width buys nothing, leaving minP-like behaviour under a maxT label.

The cost was the pre-registered anchor. AVP's family-adjusted coherence moved
from **.0120** - the reported value, and the one the p-dancing figures quote at
every perturbation - to **.1135**, which is the documented *combined*-engine
value of .125. PLAGL1, CAT and the exploratory sweep all shifted with it, in
both directions, and a gene that is not a finding (ADCY6) appeared to clear .05.

Standardisation now happens only where two engines are actually fused. With one
engine, maxT sees the raw statistic again.

**Every number in the pre-registered proof of concept now reproduces exactly**
from the shipped package at `weighting = "flat"`, `cpg_map = "full"`, B = 1999:
AVP coherence .0120/.0120, composite .0070/.1285, diffuse .0080/.0505; PLAGL1
coherence .0055/.0425, composite .0025/.0080, diffuse .0025/.1720. All twelve
raw and adjusted values match the deposit.

Three tests now pin the property: the concordant unit's toll stays small under
one engine, maxT stays strictly kinder than minP to it, and the two-engine
fusion keeps its standardisation.

# dmsa 1.8.1 (2026-08-18)

## The locus panel was mixing genome builds, silently

Found by running the whole pipeline on the real parent cohort rather than on
test fixtures - which is the only way this class of defect surfaces, because
neither layer is wrong on its own.

* **`dmsa_report()` took probe coordinates from `dmsa_probe_coords()`, which
  reads the companion package's EPIC manifest in hg19, and drew them against an
  hg38 gene model.** AVP moves 19,353 bp between the two builds, so its eight
  promoter CpGs were drawn 19 kb from their own gene, in the middle of the
  neighbouring intergenic space, with nothing reported. Coordinates now come
  from the bundled cascade, looked up by probe id: it is hg38, it is the same
  table the gene models are aligned to, and it is always present. The manifest
  lookup remains as the fallback only.
* **A genome-build guard.** `dmsa_plot_locus()` now measures the distance from
  the probes to the gene model and refuses to pass it off as a map when they
  cannot belong together: if every probe is further from the gene than the gene
  is long (floor 5 kb), it messages naming the build and prints CHECK GENOME
  BUILD on the figure itself. The rule is distance, not overlap, because real
  promoter probes legitimately sit outside the gene span - all eight of AVP's
  do - and a naive containment test would cry wolf on exactly the figure this
  feature exists to draw.

554 unit tests.

# dmsa 1.8.0 (2026-08-18)

## The locus panel draws the gene, not a bar shaped like one

Strip B used to draw a gene as a featureless grey rectangle, and - when a probe
carried no coordinate - to space probes evenly along it. Both told the reader
nothing while looking like they told them something. A CpG in the first exon
and a CpG 40 kb into an intron are different claims about regulation, and which
one it is decides the probe's direction call. The panel now shows it.

* `dmsa_gene_model()` returns real exon, UTR and strand structure on genomic
  coordinates, from any of four sources in one shape: **Ensembl REST** (no
  install, no annotation package, cached to disk), a **GTF/GFF3** you already
  have, a **TxDb** plus OrgDb, or an **EnsDb**. Every row records which.
* `dmsa_plot_locus(gene_model = ...)` draws it: coding exons as full-height
  boxes, UTRs half-height, introns as a line with chevrons running 5' to 3',
  and the TSS marked. These are the conventions every genome browser uses, so
  a reader already knows how to read the strip. `transcripts = "all"` stacks
  every isoform; the default draws the annotation's own canonical transcript
  and names it under the axis.
* `dmsa_plot_locus_gviz()` hands the same probes and the same model to **Gviz**
  for what Gviz gives and this does not: cytoband ideograms, transcript
  collapsing, and any other track stacked in the same coordinate frame.
* `dmsa_frame(gene_models = TRUE)` turns the models on for every locus panel in
  a report. Off by default, because a report should not make network calls
  without being asked; pass a prebuilt table instead to stay offline.
* `dmsa_gene_model_check()` validates a model before it is drawn, including the
  failure this is really guarding against: a region query keeps the NEIGHBOUR
  gene's exons and draws them as if they were this gene's.

Gviz, GenomicRanges, IRanges, GenomicFeatures, ensembldb, AnnotationDbi and
AnnotationFilter are **Suggests**. The Ensembl and GFF sources need none of
them - nothing beyond base R - so the figure is not reachable only by users who
can install a gigabyte of annotation.

## The evenly-spaced fallback is no longer a silent default

It looks like a map and is not one: two probes drawn adjacent may be 40 kb
apart. When the panel has no choice it now prints **NOT TO SCALE** in bold under
the axis, where the caveat survives being cropped into a slide, and messages
once on the console naming how to fix it. The bundled Alpha cascade carries an
hg38 position for every probe, and `dmsa_report()` now reads coordinates
straight from it, so the default path is always to scale.

## Fix

* The minus strand. On a gene that reads right to left the 5' UTR is the piece
  with the HIGHER coordinate, and the TSS is at the high end. Getting this
  backwards puts the promoter at the wrong end of half of all genes; it is
  pinned by tests on both strands.

# dmsa 1.7.2 (2026-08-17)

## `dmsa_sets_check("alpha")` failed the package's own data

The data sanity gate reported the bundled cascade as **not usable** on two
counts, and both were defects in the check rather than in the data. A validator
that cries wolf on the shipped set teaches users to ignore it, so:

* **Uniqueness is keyed on `probe_id`, not `cpg`.** EPIC v2 ships replicate
  designs for the same site - `cg09834343_BC21` and `_BC22` are two columns
  measuring one CpG - so 238 legitimate rows were being counted as duplicates.
  115 CpGs in the bundled set carry more than one design; they are now reported
  as a note, since replicate designs are correlated measurements of one site and
  that matters for a probe-level family.
* **Module evidence is read from the companion audit table.** It is shipped
  beside the cascade rather than inside it, so validating the raw CSV found no
  `evidence_strength` at all. `"alpha"` is now validated as it is actually
  loaded, and a cascade genuinely carrying no evidence gets a **WARN** rather
  than a FAIL - absent evidence changes how a result reads, it does not make the
  set unusable.

## Fixes found by `R CMD check --run-donttest`

* `dmsa_systems("stress|hpa")` - the form used in this package's own
  documentation - errored with "cascade file not found". A bare string that is
  neither `"alpha"` nor an existing file is now read as a search pattern, which
  is what the example always meant.
* `dmsa_polarity_fetch()`'s example moved to `\dontrun{}` and now writes to
  `tempdir()`. It made live HTTP calls during the check: 35 seconds of timeouts
  on a machine without network, and a stray CSV left in the check directory.
* `POLARITY_AUDIT.md`, `SELECTION_SETS.md` and `NEWS.md` added to
  `.Rbuildignore`; they are repository documents, not package files.

494 unit tests.

# dmsa 1.7.1 (2026-08-17)

* `dmsa_polarity()`, `dmsa_polarity_review()` and `dmsa_evidence()` accept a
  `dmsa_sets`, `dmsa_selection` or `dmsa_frame` object, not only a file path.
  The print banner recommends `dmsa_polarity_review()`, so the object already in
  the user's hand is the obvious argument; passing one used to fail with
  "cannot coerce class dmsa_selection to a data.frame".
* The module-label check in `.pol_toward_system()` uses word boundaries. Without
  them "sex de**termination**" matched the brake vocabulary and every gonadal
  module was flagged.

# dmsa 1.7.0 (2026-08-17)

## Gene-to-system polarity becomes auditable data

`w_g` is the multiplier that decides which way a system score points, so a
single wrong sign does not degrade a result - it inverts it. Until now those
signs lived inside a curated CSV with no provenance. They are now a first-class
table in which every row states where its sign came from and what would change
it.

* `dmsa_polarity()` loads the table: 1,234 genes across 30 systems, each row
  carrying `w_g`, a `role`, a `confidence`, a `w_g_source` and, for 451 rows, a
  citation. Sources: 112 curated, 121 OmniPath, 46 GtoPdb, 23 TRRUST, 1 UniProt,
  228 literature, 676 rule-based, 27 unresolved. Signs: 583 at +1, 291 at -1,
  360 held deliberately off-axis at 0.
* `dmsa_polarity_review()` surfaces the 190 rows a human still has to decide -
  low-confidence signs, unresolved genes, framing choices and the one database
  contradiction left standing - and nothing else. `which = "disagreement"` and
  `which = "unresolved"` narrow it.
* `dmsa_polarity_check()` runs the invariants: an anchor may never lower its own
  system's tone, no system may be entirely one-sided, no sign may point at its
  module's sub-process rather than the system, and no row may carry a sign
  without a source. It reports the four known, real exceptions in the Alpha
  panel rather than passing them silently.
* `dmsa_polarity_fetch()` drafts signs for a user's own cascade from GO,
  OmniPath, TRRUST and SIGNOR. When a resource is unreachable it returns
  `unresolved` at `w_g = 0`, never a guess, and flags every row for review.
* `dmsa_polarity_sources()` states each resource's caveat, including the one
  that matters most: CollecTRI and DoRothEA both assign an **activating** mode
  by default when the mode is unknown, so OmniPath's `consensus_stimulation` at
  `curation_effort = 1` is a default rather than a finding.

## The rule the audit turned on

A gene is signed against **the system's declared tone**, not against its
immediate molecular partner. FKBP5 restrains GR and GR mediates HPA negative
feedback, so FKBP5 is **+1** to the HPA system - a brake on a brake. The role
vocabulary therefore carries `brake-of-brake` and `feedback-enabler` as
first-class labels, and separates roles that commit to a direction (`driver`,
`brake`, ...) from roles whose sign is inherited from what they serve
(`accessory`, `transducer`, `capacity`, `readout`).

A database may raise the evidence grade of a sign it agrees with. It may never
flip a sign on its own. All 38 database contradictions were adjudicated one at a
time against the literature; four proved to be defects at source. `w_g` for the
Project Alpha panel is unchanged from the curated values in every one of the
112 rows the PI set by hand. `POLARITY_AUDIT.md` records the method, the
rejected resources and the escalations.

* The polarity table travels with `dmsa_sets()` and narrows with
  `dmsa_select()`, so a selection reports the coverage of its own genes.
* `print()` on a set, a selection or a frame names the unsigned and flagged
  counts instead of hiding them.

# dmsa 1.6.0

## The selection cascade becomes a first-class object

`system > module > gene > probe` now ships with the package as a curated
selection set, and it is the object `dmsa_frame()` resolves `systems =` against.
The point is not convenience: because families are declared in this file and
never harvested from whatever the analyst happened to load, the family a unit is
corrected within is a property of the map rather than of the session. That is the
mechanism behind DMSA's invariance to set selection, and it is now visible.

* `dmsa_sets()` loads a selection cascade - the bundled, module-audited Project
  Alpha 2026c set (30 systems, 188 modules, 1,234 genes, 16,823 CpGs) or any CSV
  in the same shape.
* `dmsa_systems()` prints the systems with the **short name** for each, which is
  what `systems =` now takes: `systems = c("hpa", "oxytocin")`. Ids, full names
  and unique prefixes resolve too; ambiguity errors with the candidates named
  rather than picking one.
* **Everything below the chosen system defaults to full.** `"hpa"` takes all 8
  modules, 49 genes and 538 CpGs the set assigns to it. `dmsa_select()` narrows a
  level when a narrower question was declared, and says so in its print output,
  because narrowing changes the family and therefore the toll.
* `dmsa_select(columns = "auto", data = dat)` keys the cascade to the caller's
  matrix through the `col_*` columns instead of requiring renamed columns.
* `dmsa_frame()` gains `sets =` (default `"alpha"`). The module level now takes
  its membership from the cascade rather than the legacy `modules_alpha.csv`, so
  the family a unit is corrected within is the family the user selected from.
  `sets = NULL` restores the legacy path.

## Module evidence travels with the result

Every module label in the bundled set was checked against the literature. That
check is now part of the output, on the principle that a module-level finding is
only as good as the module's definition.

* `dmsa_evidence()` returns the evidence tier (High / Moderate), the audit status
  (`supported_as_named`, `renamed_*`, `heterogeneous_locked_membership`,
  `measurement_defined_module`), the citation keys and the evidence note, for a
  set, a selection or a frame. `which = "flagged"` reduces it to the modules that
  need a caveat when reported.
* The counts print with `dmsa_sets()`, `dmsa_select()` and `dmsa_frame()`, and
  `summary.md` from `dmsa_report()` gains a **Module evidence** table covering the
  modules the run actually reported on.

## Bringing your own panel

* `dmsa_sets_template()` writes the CSV skeleton plus the full instructions;
  `dmsa_sets_check()` validates a candidate file against the five structural
  rules before it is trusted (one system per module, one module per gene, unique
  short names, constant module metadata, no duplicate rows) and reports CpGs
  annotated to more than one gene rather than rejecting them.
* `SELECTION_SETS.md` documents the schema, the resolution rules and the
  separation between the selection set, the direction map and the polarity table.

## Fixes

* `module_id` is read as text throughout. Parsed as a number, `24.10` collapses
  onto `24.1` and two different immune modules become one - a silent, total
  corruption of the family structure.
* `systems =` no longer errors when a system id is character rather than numeric.
* `data.table` moved from an undeclared `::` call to `Imports`.
* Documented the previously undocumented `weighting`, `w_floor`, `gate`,
  `sparse_reach`, `family_correction` and `nulls` arguments.

# dmsa 1.5.0 (2026-08-15)

* ENGINE DEFAULT CONFIRMED: `combined`. The evidence points two ways and both
  directions are on the record. On REAL data, pre-registered application to the
  PTSD and children cohorts showed neither flat nor reliability weighting
  dominates - reliability gains on concentrated single-gene signals (FKBP5,
  SGK1) and loses on diffuse system and module signals - so hedging the engine
  choice keeps every real finding. In SIMULATION, however, flat weighting needs
  FEWER participants than combined in all fourteen sample-size cells, by 3% at
  a twelve-gene family (326 vs 337) and by 62% on a three-system family (472 vs
  764). The default is `combined` because the shape of a real signal is not
  known in advance and the premium for that insurance is roughly 10% more
  participants in a typical cell; `weighting = "flat"` is a first-class option
  for anyone who knows their signal is coherent, and remains the exact
  reproducibility reference for the published anchors.

* Documentation rendered: `man/` now ships 48 Rd files, so `?dmsa_frame` and
  the rest resolve. The package previously carried roxygen comments with no
  generated documentation.

* `inst/scripts/artifact_probe_mask.R` - checks EPIC probes carrying a finding
  against five published artifact lists (cross-reactive, cross-hybridising,
  SNP at target, INDEL at target, 1000G variant in probe body), with the
  provenance of every list documented in the header.

* Two tests in `test-design.R` now qualify an internal call as `dmsa:::`, so
  the suite runs outside the package namespace. 343 tests, 0 failures.

* FIXED: the combined engine's rank-p fusion had a DISCRETENESS FLOOR. Rank
  p-values are granular at 1/(B+1); under maxT across K units, ties at the
  extreme floored the family-wise p at about K/(B+1). At K = 42, B = 399 that
  is .105, so power was exactly 0.000 regardless of effect size (found on an
  independent K = 42 sweep, 15 Aug 2026; the defect was invisible in our own
  tests, which never used a family larger than five). Engines are now fused on
  a CONTINUOUS scale: each lens/engine statistic is standardised by its own
  permutation null (robust centre/scale) before combination, so no ties and no
  floor. Verified at K = 42, B = 399.

* NEW: `p_unit` / `p_unit_adj` - a unit-level test on ONE joint null across
  lenses (and engines). The per-lens adjusted p-values correct across UNITS
  only, so the previously used rule - "name a unit if the minimum of the three
  lens-adjusted p's is below .05" - is not family-wise controlled across
  lenses and inflates the error rate roughly threefold (measured: .06 where
  .05 was intended, and up to .13 in some grids). `p_unit_adj` takes the
  maximum standardised lens statistic per unit and applies maxT over units:
  one rule, valid across lenses and units. Measured family-wise error .030
  (flat) and .040 (combined) at K = 12 against the old rule's .060. USE
  `p_unit_adj` FOR NAMING; the per-lens columns remain for reporting which
  lens carried a finding.

# dmsa 1.4.0 (2026-08-15)

* THREE ENGINES; the COMBINED engine is now the DEFAULT. Pre-registered
  testing on the PTSD and children cohorts showed neither flat nor
  reliability dominates on real data: reliability gains on concentrated
  single-gene signals (FKBP5, SGK1) but loses on diffuse system/module
  signals (the GR-module PTSD lead, the child HPA-system finding). The
  combined engine fuses the two per lens by an ACAT (Cauchy) combination of
  their rank-p's on ONE shared permutation stream, so the family-wise maxT is
  an exact JOINT NULL - not a post-hoc ACAT of two adjusted p's. It keeps
  every real signal (validated: GR module .0445 flat / .118 reliability /
  .064 combined; FKBP5 .681 / .189 / .368) for the small toll of hedging the
  engine choice, in the same spirit as the three-lens triangulation. Both
  per-permutation regressions are shared across engines (computed once), so
  combined costs ~2x lens work, not 2x total. `weighting = c("combined",
  "reliability", "flat")` throughout; flat remains the exact reproducibility
  reference. Calibration of the combined omnibus verified (outcome-free
  weights + pooled-rank fusion keep type I at alpha). Frame print names the
  active engine. +2 assertions (combined calibration + both-archetype
  recovery).

# dmsa 1.3.0 (2026-08-15)

* TWO ENGINES, reliability weighting now the DEFAULT. Each aligned probe is
  weighted by its item-rest (centrality) correlation with the rest of its
  unit - `dmsa_relweights()`, computed from methylation ALONE so the
  permutation null stays exact and weighting changes power, not type-I error.
  In simulation it gains power on units whose probe quality is heterogeneous
  and ties the flat engine when it is uniform; it therefore ships as the
  default. `weighting = "flat"` restores the original equal-weight engine and
  is the exact reproducibility reference. Threaded through `dmsa_frame()`,
  `dmsa_report()`, `dmsa_triangulate()`, `dmsa_scores()` and `dmsa_change()`
  via `weighting` and `w_floor` (leading-eigenvalue coherence floor, default
  1.5). Guards: units with < 3 usable probes, no shared axis, or all-zero
  clamped weights fall back to flat. `dmsa_test()` gains a `w` argument.
  Frame print reports the active engine. +7 assertions (test-weight.R):
  calibration, flat reproducibility, power gain on heterogeneous units,
  guard fallbacks.

# dmsa 1.2.0 (2026-08-14)

* `dmsa_change()` - longitudinal mDMSA on the 2-wave identity: time x S x E
  == Delta-S x E. Per-unit aligned scores at each wave from the SAME map,
  Delta-score x exposure through the composite machinery (Freedman-Lane,
  family blocks), maxT/minP across the level-local family on one shared
  permutation stream; baseline score as regression-to-the-mean guard;
  `mod2` adds the second moderator (Delta-S x E x mod2). Validated on the
  Oct-7 cohort (n = 172) before shipping, per the pre-set condition.
* The frame layout `moderation`/`mod`/`mod2` now carries the three-way end
  to end (planted-interaction test locked in the suite); in a longitudinal
  frame Delta IS time, so `mod`/`mod2` are the two remaining factors.
* Frame fixes from the real child-data run: compositional check no longer
  trips on zero-variance aliased columns; `chip_f` autofix accepts any
  `chip`/`chip_T*` column; factor covariates and chip are droplevels-ed
  after complete-case subsetting (empty levels made designs singular).

# dmsa 1.1.2 (2026-08-14)

* A completed moderation run with zero survivors now SAYS so: print and
  summary report units tested, families, B, and the strongest interaction
  (unit, b, t, raw and adjusted p) - a calibrated null no longer looks like
  a run that did nothing. Sub-minute runtimes print in seconds.

# dmsa 1.1.1 (2026-08-14)

* Fix: the untestable-gene count in `frame$untestable` and the report Design
  notes was computed across ALL systems in the data instead of the chosen
  `systems` subset (spotted by comparing the Mac and cloud anchor runs -
  the statistics were never affected, only the coverage sentence).

# dmsa 1.1.0 (2026-08-14)

The two-call user interface.

* `dmsa_frame()` - declare everything once: data, levels (system / module /
  gene / probe), outcomes, covariates ("contract" default), random-effect
  blocks, moderation, map, styling. Runs the TEST DRIVE automatically
  (autofix = TRUE): numeric coercion with per-column loss report, QR rank
  audit dropping aliased columns by name (compositional groups flagged),
  block audit with equal-size permutation strata, confounded-block flags,
  hard error only when fewer than ~2^10 distinct permutations exist, binary
  outcome under gaussian autofixed to logistic, B = 49 pilot with a printed
  ETA. Every action lands in `frame$corrections`, printed by `print(frame)`
  and repeated in the report's Design notes.
* `dmsa_report(frame)` - runs `dmsa_triangulate()` at each requested level
  with correction ONLY inside level-local families (systems named; modules of
  the system; genes of the system; probes of surviving genes), writes
  figures/ tables/ summary.md under `frame$outdir`, and returns the paths.
* Moderation layout: `moderation = FALSE` (default), `mod = ""`, `mod2 = ""`.
  `moderation = TRUE` + `mod` tests frame x mod through the composite lens
  (the one lens that carries to products), maxT across units on a shared
  permutation stream; `mod2` declares frame x mod x mod2 - the 3-way engine
  ships in 1.2.0 after tonight's longitudinal test drive, and 1.1.0 raises a
  clear error if `mod2` is set.
* `frame_role = "predictor" | "outcome"` - methylation as predictor (default)
  or as the outcome of the named columns; the composite lens flips the score
  to the left-hand side, narration switches to expression-tone language.
* `cpg_map = "confidence" | "full"` - both alignments are always built; any
  gene or system whose pooled direction flips between maps is analysed under
  the chosen map but warned about, stored in `frame$map_conflicts`, and
  repeated in the report.
* Brain-bridge post-hoc, always on: surviving units' probes go through
  `cpgdirection::cpg_brain_bridge()`; hit rule fixed as `bridge_usable |
  has_t2_direction | !is.na(bridge_grade)`; `ensemble_pctile` is context.
  Errors never kill a report.
* `dmsa_triangulate()` now returns NA (not a p = 1 placeholder) for units
  with zero usable direction-called probes; `dmsa_frame()` tracks them in
  `frame$untestable` and the report counts them per family.
* `dmsa_model(nulls = TRUE)` returns its permutation null (`null_t`) for
  family maxT across separately fitted models.
* The former figure helper `dmsa_report()` is now `dmsa_report_panels()`.
* gt moved into Suggests: tables use gt when installed, otherwise a
  self-contained HTML fallback (CSVs are always written).

Anchor (real Alpha data, systems 1 + 21, B = 1999, seed 1, cpg_map = "full"):
AVP coherence adj .0120, PLAGL1 composite adj .0080 - identical to the
locked triangulation run. New from the always-on bridge: two PLAGL1 probes
carry a formal brain bridge (cg08263357, ensemble percentile 94; cg11532302).
