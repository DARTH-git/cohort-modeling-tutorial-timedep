# cohort-modeling-tutorial-timedep 0.3.0 (2026-11-10)

This release fixes a bug in the Simpson's 1/3rd within-cycle correction (WCC).
**It changes the numerical results reported in the published article**, but
not its conclusions. It also fixes the probabilistic sensitivity analysis
(PSA) of the state-residence model, which stopped with an error.

## Bug fixes

* `gen_wcc(method = "Simpson1/3")` (`R/Functions.R`) assigned the 4/3 and 2/3
  weights by position in the 1-based R vector instead of by time
  t = 0, 1, ..., n. It returned `1/3, 2/3, 4/3, 2/3, ..., 1/3` instead of
  `1/3, 4/3, 2/3, 4/3, ..., 1/3`, so the weights summed to n - 2/3 (even n)
  or n - 1/3 (odd n) instead of n. Weights are now indexed by time. This is
  the same fix as in
  [cohort-modeling-tutorial-intro v0.3.0](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/releases/tag/v0.3.0)
  and [darthtools v0.4.1](https://github.com/DARTH-git/darthtools/releases/tag/v0.4.1).

  The same error is in the `within-cycle-vector` code chunk of the code
  walkthrough appendix (Supplementary Material), and in the introductory
  tutorial that the article refers to for the WCC. The files in `manuscript/`
  are kept unchanged as the archival source of the article (see
  `manuscript/README.md`).

  Impact on the base-case results (both models use 75 annual cycles):

  | Model | Strategy | QALYs (old → new) | Costs (old → new) | ICER (old → new) |
  |---|---|---|---|---|
  | Simulation-time | SoC | 18.879 → 19.199 | $116,374 → $117,126 | — |
  | | B | 20.199 → 20.519 | $202,536 → $203,707 | $65,288 → $65,611 |
  | | AB | 21.097 → 21.423 | $296,300 → $297,850 | $104,461 → $104,125 |
  | | A | 19.636 → 19.962 | $218,789 → $219,923 | Dominated |
  | State-residence | SoC | 20.005 → 20.325 | $103,112 → $103,872 | — |
  | | B | 21.104 → 21.423 | $186,262 → $187,436 | $75,678 → $76,060 |
  | | AB | 22.098 → 22.424 | $274,408 → $275,960 | $88,668 → $88,477 |
  | | A | 20.882 → 21.208 | $198,058 → $199,197 | Dominated |

  QALYs increase by about 0.32 per strategy and ICERs change by less than
  0.6%. The efficient frontier (SoC, B, AB) and the dominance of Strategy A
  are unchanged.

* The PSA of the state-residence model (`analysis/cSTM_time_dep_state_residence.R`)
  stopped at parameter set 220 with "This is not a valid transition array".
  `check_sum_of_transition_array()` compared the sum of each cycle's
  transition probabilities with the number of states using exact equality,
  so rounding error (about 1e-14) in a valid array made it fail. It now uses
  a numerical tolerance (`all.equal()`), as in darthtools.
  The script then stopped when combining the PSA figures, because it did not
  load `patchwork`; it now does.

* `calculate_ce_out()` of the simulation-time model
  (`R/Functions_cSTM_time_dep_simulation.R`) discounted costs with `d_e` and
  QALYs with `d_c`. The analysis scripts were fixed in 2024, but the PSA
  function was not. This had no numerical effect in the tutorial, because both
  rates are 3%.

* `calculate_ce_out()` of the state-residence model used the WCC and discount
  vectors (`v_wcc`, `v_dwc`, `v_dwe`) from the global environment instead of
  computing them. The model functions of both models also read the strategy
  and state names from the global environment. They now define these
  themselves, so the PSA functions can be used on their own.

## Changes in behaviour

* For an odd number of cycles (such as the tutorial's 75), `gen_wcc()` now
  applies the composite Simpson's 1/3 rule over the first n - 3 cycles and
  Simpson's 3/8 rule over the last three. Both rules integrate polynomials up
  to degree three exactly, and the weights always sum to `n_cycles`.
* `gen_wcc()` now errors when `n_cycles` is not a single positive integer, and
  when `n_cycles = 1` is used with Simpson's rule.

## Infrastructure

* Added unit and regression tests with testthat (`tests/testthat/`) for
  `gen_wcc()`, the base-case results of both models, the PSA functions and the
  PSA input distributions.
* Added a GitHub Actions workflow (`.github/workflows/check.yaml`) that runs
  the tests and both analysis scripts end to end on Ubuntu and macOS.
* Added `renv` (`renv.lock`) to pin package versions, including
  `MethylCapSig` 1.0.1, which has been archived on CRAN.
* Added a `DESCRIPTION` file with the project version and dependencies, and
  `CITATION.cff` with citation metadata for the code and the article.
* Replaced deprecated ggplot2 syntax (`size` → `linewidth` in `geom_line()`,
  `..density..` → `after_stat(density)`).
* Removed a stray anonymous function in `R/Functions.R`.

## Documentation

* Documented two errata in the published article in `manuscript/README.md`
  and the README: the WCC weights, and the cost of Treatment B (`c_trtB`) in
  the parameter table, which lists $12,000 instead of the $13,000 used in the
  text, the code and the published results.
* Rewrote the README: quick start with `renv`, repository structure, the
  corrections above, how to report issues, and acknowledgements.

## Acknowledgements

* Lauren N. Steimle (Georgia Institute of Technology) and GitHub user
  [@ambtnhlr](https://github.com/ambtnhlr) independently found and reported
  the WCC indexing error in the introductory tutorial
  ([DARTH-git/cohort-modeling-tutorial-intro#2](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/issues/2)),
  which led to this fix.
* GitHub user [@conallaque](https://github.com/conallaque) found the
  `c_trtB` parameter-table erratum while replicating both tutorials
  ([DARTH-git/cohort-modeling-tutorial-intro#1](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/issues/1)).

# cohort-modeling-tutorial-timedep 0.2.0 (2022-06-07)

* Code accompanying the resubmitted manuscript and arXiv v2
  (arXiv:2108.13552v2), published as Alarid-Escudero et al. (2023),
  *Medical Decision Making* 43(1):21-41,
  <https://doi.org/10.1177/0272989X221121747>.
  Zenodo: <https://doi.org/10.5281/zenodo.6620902>.
* Later changes on `main` before v0.3.0: README updates, and a fix for
  swapped cost and QALY discount rates in the analysis scripts (2024).

# cohort-modeling-tutorial-timedep 0.1.0 (2021-08-30)

* First release, accompanying the arXiv preprint (arXiv:2108.13552v1).
  Zenodo: <https://doi.org/10.5281/zenodo.5338819>.
