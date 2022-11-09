<h2> Power and Sample Size Calculation Tools </h2>

To install and load `pwrss`:

`install.packages("pwrss")`

`library(pwrss)`

`pwrss` R package allows statistical power and minimum required sample size calculations for 

  - `(1)` testing a proportion (one-sample) against a constant, 
  - `(2)` testing a mean (one-sample) against a constant,  
  - `(3)` testing difference between two proportions (independent samples),
  - `(4)` testing difference between two means (independent and paired samples),
  - `(5)` testing a correlation (one-sample) against a constant,  
  - `(6)` testing difference between two correlations (independent samples),
  - `(7)` testing an R-squared against zero in linear regression 
  - `(8)` testing an R-squared difference against zero in hierarchical regression,
  - `(9)` testing an eta-squared or f-squared (for main and interaction effects) against zero in analysis of variance (ANOVA) (could be one-way, two-way, and three-way), 
  - `(10)` testing an eta-squared or f-squared (for main and interaction effects) against zero in analysis of covariance (ANCOVA) (could be one-way, two-way, and three-way), 
  - `(11)` testing an eta-squared or f-squared (for between, within, and interaction effects) against zero in one-way repeated measures analysis of variance (RM-ANOVA) (with non-sphericity correction and repeated measures correlation). 
  
  Alternative hypothesis can be formulated as "not equal", "less", "greater",  "non-inferior", "superior", or "equivalent" in `(1)`, `(2)`, `(3)`, and `(4)`;
  as "not equal", "less", or "greater" in `(5)` and `(6)`; but always as "greater" in `(7)`, `(8)`, `(9)`, `(10)`, and `(11)`. 
  
  
**If you find the package and related material useful, along with other references, please cite as:**

- Bulus, M. (2022). pwrss: Power and Sample Size Calculation Tools. R package version 0.1.0. https://CRAN.R-project.org/package=pwrss

- Bulus, M., & Polat, C. (2022, November 7). Istatistiksel Guc Analizine Giris [Introduction to Statistical Power Analysis]. https://doi.org/10.35542/osf.io/tfyxq 




