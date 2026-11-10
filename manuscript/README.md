# Manuscript sources (archival)

This folder holds the R Markdown source, LaTeX output, figures and appendices
of the published article:

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P,
> Jalal H. A Tutorial on Time-Dependent Cohort State-Transition Models in R
> using a Cost-Effectiveness Analysis Example. *Medical Decision Making*.
> 2023;43(1):21-41. <https://doi.org/10.1177/0272989X221121747>

The files are kept **unchanged** so that they match the published article.

## Known erratum: within-cycle correction weights

The article uses Simpson's 1/3rd rule for the within-cycle correction (WCC)
and refers to the introductory tutorial for its implementation. The
`within-cycle-vector` code chunk of the code walkthrough appendix
(`Appendix_code_walkthrough_Markov_Tutorial_Part2.Rmd`) and the code used for
the published results assign Simpson's 1/3rd weights by position in the
1-based R vector instead of by time t = 0, 1, ..., n:

* The chunk applies the even/odd rule to `v_cycles <- seq(1, n_cycles + 1)`,
  which are positions, not times. It produces
  **wcc = [1/3, 2/3, 4/3, ..., 1/3]**. The correct vector is
  **wcc = [1/3, 4/3, 2/3, 4/3, ..., 4/3, 1/3]**: 4/3 at odd times
  t = 1, 3, ..., and 2/3 at even times t = 2, 4, .... The fix is to use
  `v_cycles <- 0:n_cycles`.
* The tutorial uses 75 cycles. That is an odd number of intervals, so
  Simpson's 1/3rd rule alone does not apply. The corrected code uses Simpson's
  3/8 rule for the last three cycles.

As a result, the QALYs and costs reported in the article are slightly
underestimated: by about 0.32 QALYs per strategy, with ICERs off by less than
0.6% (for example, the ICER of Strategy B in the simulation-time-dependent
model is $65,611/QALY rather than the published $65,288/QALY). The
conclusions are unchanged.

The corrected implementation is `gen_wcc()` in [`../R/Functions.R`](../R/Functions.R),
from release v0.3.0 on. It is used by both analysis scripts in
[`../analysis/`](../analysis/). See [`../NEWS.md`](../NEWS.md) for details and a
before/after comparison of results.

The error was identified and reported independently by Lauren N. Steimle
(Georgia Institute of Technology) and GitHub user
[ambtnhlr](https://github.com/ambtnhlr)
([DARTH-git/cohort-modeling-tutorial-intro#2](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/issues/2))
in the introductory tutorial, which shares the code. We thank them both.

## Known erratum: cost of Treatment B in the parameter table

The parameter table (Table 2) lists the cost of Treatment B, `c_trtB`, as
**$12,000**. The correct value is **$13,000**. That is the value stated in the
article's text ("costs $13,000 per year"), used in the code, implied by the
PSA distribution in the same row (gamma(86.2, 150.8), whose mean is
86.2 × 150.8 ≈ $13,000), and used to produce all the published results.

The error is only in the table, in `cSTM_Tutorial_TimeDep.Rmd` and the derived
`.tex` and `.pdf` files. The code in this repository is correct, and a unit
test (`tests/testthat/test-psa-params.R`) checks that the PSA distributions
are centred on the base-case values.

This was identified by GitHub user [conallaque](https://github.com/conallaque),
who replicated both tutorials independently
([DARTH-git/cohort-modeling-tutorial-intro#1](https://github.com/DARTH-git/cohort-modeling-tutorial-intro/issues/1)).
We thank them.
