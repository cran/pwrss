## Changes in pwrss v1.2.0
  - Added unit tests for all procedures
  - Added minimum detectable effects
  - Added exact correlation (one sample)
  - Added type S and M errors for two-sided tests, available through `plot()` function
  - Added `clean.parms`, unifying how the input parameters in the results list are returned
  - Added further input checks so that all input parameters are checked
  - Simplified code in `power.z.logistic` and `power.z.poisson`
  - Corrected formula for Hsieh et al. (1998) in `power.z.logistic`
  - Added `target.effect` to determine effect of interest in `power.f.ancova`

## Changes in pwrss v1.0.0
  - Addition 1: ANCOVA procedures with means and contrasts (Shieh)
  - Addition 2: Exact tests for proportions (Fisher & McNemar)
  - Addition 3: Generic binomial test (with plotting feature)
  - Addition 4: Dependent correlations (Steiger)
  - Addition 5: Effect size conversion functions 
  - Extension 1: Generic tests handles non-zero null (useful for equivalence and minimum effect testing). Some specific tests also allow non-zero null. 
  - Extension 2: ANCOVA procedures allow R-squared adjustment (i.e., Keppel and Shieh)
  - Consistent function names, population-consistent argument names (Greek notation) and base R consistency
  - Rigorous error checks and warnings 
  - Comprehensive print utility with three-tier verbosity control
  - Minor bug fixes to `plot()` function

## Changes in pwrss v0.3.2 (not released to CRAN)
  - alternative = "less" now produce correct power rates in `pwrss.z.prop()` function. Thanks to Leszek Gawarecki for reporting the issue. 
  - The default for p0 argument in `pwrss.chisq.gofit()` can now be overwritten. Thanks to Kate Crespi for reporting the issue. 
  - Type 1 and type 2 error plots show and print correct power estimates for two-tailed tests when effect size is near zero. Thanks to dpnichols811 (GitHub profile handle) and Adrian Olszewski for reporting the issue. 
  - `pwrss.z.corr()` print correct errors with alternative = "greater" when in fact it is less (or vice versa). Thanks to Jarrod Hadfield for reporting the issue. 

## Changes in pwrss v0.3.1
  - More detailed vignette (typo fixes)
  - Added web application links (shiny dashboards)
  - Minor bug fixes with `plot()` function

## Changes in pwrss v0.3.0
  - Added `pwrss.np.2means()`, `pwrss.z.logreg()`, `pwrss.z.poisreg()`, `pwrss.chisq.gofit()`, and generic `power.chisq.test()` functions 
  - Improvements to `plot()` function
  - Improvements to documentation
  - Improvements to `pwrss.z.med()` and `pwrss.f.ancova()` console print
  - More detailed and comprehensive examples in vignette
  - Bug fixes (Welch’s t test)
 
## Changes in pwrss v0.2.0
 - Added `pwrss.t.reg()` and `pwrss.z.med()` functions 
 
## pwrss v0.1.0
 - Initial release 

