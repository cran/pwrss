#' Inflate Sample Size for Attrition
#'
#'
#' @param n       sample size.
#' @param rate    attrition rate.
#' @param ceil.n  rounds-up the inflated sample size.
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#'                if \code{2} a more detailed output is given (plus key
#'                parameters and definitions), if \code{0} no output is printed
#'                on the console.
#'
#' @return inflated sample size.
#'
#' @examples
#' inflate.sample(n = 100, rate = 0.05)
#'
#' @export inflate.sample
inflate.sample <- function(n, rate = 0.05, ceil.n = TRUE, verbose = 1) {

  check.sample.size(n)
  verbose <- ensure.verbose(verbose)

  n.adj <- ifelse(ceil.n, ceiling(n / (1 - rate)), n / (1 - rate))

  if (verbose > 0)
    cat(n.adj)

  invisible(n.adj)

} # inflate.sample


#' Conversion from Eta-squared to Cohen's f
#'
#'
#' @param eta.squared (Partial) Eta-squared.
#' @param verbose     \code{1} by default (returns test, hypotheses, and
#'                    results), if \code{2} a more detailed output is given
#'                    (plus key parameters and definitions), if \code{0} no
#'                    output is printed on the console.
#' @return
#'   \item{f}{Cohen's f.}
#'   \item{f.squared}{Cohen's f².}
#'   \item{eta.squared}{(Partial) Eta-squared.}
#'
#' @references
#'   Cohen, J. (1988). *Statistical power analysis for the behavioral sciences*
#'   (2nd ed.). Lawrence Erlbaum Associates.
#'
#' @examples
#' etasq.to.f(eta.squared = 0.01) # small
#' etasq.to.f(eta.squared = 0.06) # medium
#' etasq.to.f(eta.squared = 0.14) # large
#'
#' @export etasq.to.f
etasq.to.f <- function(eta.squared, verbose = 1) {

  check.nonnegative(eta.squared)
  verbose <- ensure.verbose(verbose)

  f.squared <- eta.squared / (1 - eta.squared)

  if (verbose > 0)
    print(c(f.squared = f.squared, f = sqrt(f.squared), eta.squared = eta.squared))

  invisible(list(f.squared = f.squared, f = sqrt(f.squared), eta.squared = eta.squared))

} # etasq.to.f


#' Conversion between Cohen's f and Eta-squared
#'
#'
#' @param f           Cohen's f.
#' @param verbose     \code{1} by default (returns test, hypotheses, and
#'                    results), if \code{2} a more detailed output is given
#'                    (plus key parameters and definitions), if \code{0} no
#'                    output is printed on the console.
#' @return
#'   \item{eta.squared}{(Partial) Eta-squared.}
#'   \item{f.squared}{Cohen's f².}
#'   \item{f}{Cohen's f.}
#'
#' @references
#'   Cohen, J. (1988). *Statistical power analysis for the behavioral sciences*
#'   (2nd ed.). Lawrence Erlbaum Associates.
#'
#' @examples
#' f.to.etasq(f = 0.10) # small
#' f.to.etasq(f = 0.25) # medium
#' f.to.etasq(f = 0.40) # large
#'
#' @export f.to.etasq
f.to.etasq <- function(f, verbose = 1) {

  check.nonnegative(f)
  verbose <- ensure.verbose(verbose)

  f.squared <- f ^ 2
  eta.squared <- f.squared / (1 + f.squared)

  if (verbose > 0)
    print(c(eta.squared = eta.squared, f.squared = f.squared, f = sqrt(f.squared)))

  invisible(list(eta.squared = eta.squared, f.squared = f.squared, f = sqrt(f.squared)))

} # f.to.etasq


#' Conversion from R-squared to Cohen's f
#'
#' @description
#' Helper function to convert between Cohen's f and R-squared.
#'
#' @param r.squared.full    R-squared for the full model.
#' @param r.squared.reduced R-squared for the reduced model.
#' @param verbose           \code{1} by default (returns results), if \code{0}
#'                          no output is printed on the console.
#'
#' @return
#'   \item{f.squared}{Cohen's f-squared.}
#'   \item{f}{Cohen's f.}
#'   \item{r.squared.full}{R-squared for the full model.}
#'   \item{r.squared.reduced}{R-squared for the reduced model.}
#'
#' @references
#'   Cohen, J. (1988). Statistical power analysis for the behavioral sciences
#'   (2nd ed.). Lawrence Erlbaum Associates.
#'
#'   Selya, A. S., Rose, J. S., Dierker, L. C., Hedeker, D., & Mermelstein,
#'   R. J. (2012). A practical guide to calculating Cohen's f2, a measure of
#'   local effect size, from PROC MIXED. *Frontiers in Psychology, 3*, Article
#'   111. https://doi.org/10.3389/fpsyg.2012.00111
#'
#' @examples
#'
#'   rsq.to.f(r.squared.full = 0.02) # small
#'   rsq.to.f(r.squared.full = 0.13) # medium
#'   rsq.to.f(r.squared.full = 0.26) # large
#'
#' @export rsq.to.f
rsq.to.f <- function(r.squared.full, r.squared.reduced = 0, verbose = 0) {

  check.nonnegative(r.squared.full, r.squared.reduced)
  verbose <- ensure.verbose(verbose)

  if (r.squared.full < r.squared.reduced)
    stop("Expecting `r.squared.full` > `r.squared.reduced`.", call. = FALSE)

  f.squared <- (r.squared.full - r.squared.reduced) / (1 - r.squared.full)

  if (verbose > 0)
    print(c(f.squared = f.squared, f = sqrt(f.squared), r.squared.full = r.squared.full, r.squared.reduced = r.squared.reduced))

  invisible(list(f.squared = f.squared, f = sqrt(f.squared),
                 r.squared.full = r.squared.full,
                 r.squared.reduced = r.squared.reduced))
} # rsq.to.f


#' Conversion from Cohen's f to R-squared
#'
#' @description
#' Helper function to convert between Cohen's f and R-squared.
#'
#'
#' @param f              Cohen's f.
#' @param r.squared.full R-squared for the full model.
#' @param verbose        \code{1} by default (returns results), if \code{0} no
#'                       output is printed on the console.
#'
#' @return
#'   \item{f}{Cohen's f.}
#'   \item{f.squared}{Cohen's f-squared.}
#'   \item{r.squared.full}{R-squared for the full model.}
#'   \item{r.squared.reduced}{R-squared for the reduced model.}
#'
#' @references
#'   Cohen, J. (1988). Statistical power analysis for the behavioral sciences
#'   (2nd ed.). Lawrence Erlbaum Associates.
#'
#'   Selya, A. S., Rose, J. S., Dierker, L. C., Hedeker, D., & Mermelstein,
#'   R. J. (2012). A practical guide to calculating Cohen's f2, a measure of
#'   local effect size, from PROC MIXED. *Frontiers in Psychology, 3*, Article
#'   111. https://doi.org/10.3389/fpsyg.2012.00111
#'
#' @examples
#'
#'   f.to.rsq(f = 0.14) # small
#'   f.to.rsq(f = 0.39) # medium
#'   f.to.rsq(f = 0.59) # large
#'
#' @export f.to.rsq
f.to.rsq <- function(f, r.squared.full = NULL, verbose = 0) {

  check.nonnegative(f)
  verbose <- ensure.verbose(verbose)

  f.squared <- f ^ 2
  if (is.null(r.squared.full)) {
    r.squared.full <- f.squared / (1 + f.squared)
    r.squared.reduced <- 0
  } else {
    change.r.squared <- f.squared * (1 - r.squared.full)
    r.squared.reduced <- r.squared.full - change.r.squared
  }

  if (verbose > 0)
    print(c(f.squared = f.squared, f = sqrt(f.squared), r.squared.full = r.squared.full, r.squared.reduced = r.squared.reduced))

  invisible(list(f.squared = f.squared, f = sqrt(f.squared),
                 r.squared.full = r.squared.full,
                 r.squared.reduced = r.squared.reduced))
} # f.to.rsq


#' Conversion from a correlation to a z-value (Fisher's z-transformation)
#'
#'
#' @param rho     correlation
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#'                if \code{2} a more detailed output is given (plus key
#'                parameters and definitions), if \code{0} no output is printed
#'                on the console.
#'
#' @return
#'   \item{z}{z-value}
#'   \item{rho}{correlation}
#'
#' @examples
#' cor.to.z(rho = 0.1)
#' cor.to.z(rho = 0.5)
#' cor.to.z(rho = 0.9)
#' cor.to.z(rho = 0.99)
#'
#' @export cor.to.z
cor.to.z <- function(rho, verbose = 1) {

  check.vector(rho, check.correlation, min.length = 1)
  verbose <- ensure.verbose(verbose)

  z <- vapply(rho, function(v) log((1 + v) / (1 - v)) / 2, numeric(1))

  if (verbose > 0)
    print(c(z = z, rho = rho))

  invisible(list(z = z, rho = rho))

} # cor.to.z()


#' Conversion from a z-value to a correlation (inverse Fisher's z-transformation)
#'
#'
#' @param z       z-value
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#'                if \code{2} a more detailed output is given (plus key
#'                parameters and definitions), if \code{0} no output is printed
#'                on the console.
#'
#' @return
#'   \item{rho}{correlation}
#'   \item{z}{z-value}
#'
#' @examples
#' z.to.cor(z = 0.1003353)
#' z.to.cor(z = 0.5493061)
#' z.to.cor(z = 1.4722193)
#' z.to.cor(z = 2.6466524)
#'
#' @export z.to.cor
z.to.cor <- function(z, verbose = 1) {

  check.vector(z, check.numeric, min.length = 1)
  verbose <- ensure.verbose(verbose)

  rho <- vapply(z, function(v) (exp(2 * v) - 1) / (exp(2 * v) + 1), numeric(1))

  if (verbose > 0)
    print(c(rho = rho, z = z))

  invisible(list(rho = rho, z = z))

} # z.to.cor()


# Cohen (1988, S. 109)
# q < .1: no effect;
# .1 <= q < .3: small effect;
# .3 <= q < .5: intermediate effect;
# q =>.5: large effect

#' Conversion from a correlation Difference to Cohen's q
#'
#'
#' @param rho1    first correlation.
#' @param rho2    second correlation.
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#'                if \code{2} a more detailed output is given (plus key
#'                parameters and definitions), if \code{0} no output is printed
#'                on the console.
#'
#' @return
#'   \item{q}{Cohen's q effect size.}
#'   \item{delta}{correlation difference: rho1 - rho2.}
#'   \item{rho1}{first correlation.}
#'   \item{rho2}{second correlation.}
#'
#' @examples
#' cors.to.q(rho2 = 0.5712027, rho1 = 0.50)
#' cors.to.q(rho2 = 0.6907068, rho1 = 0.50)
#' cors.to.q(rho2 = 0.7815365, rho1 = 0.50)
#'
#' @export cors.to.q
cors.to.q <- function(rho1, rho2, verbose = 1) {

  check.correlation(rho1, rho2)
  verbose <- ensure.verbose(verbose)

  q <- cor.to.z(rho1, FALSE)$z - cor.to.z(rho2, FALSE)$z

  if (verbose > 0)
    print(c(q = q, delta = rho1 - rho2, rho1 = rho1, rho2 = rho2))

  invisible(list(q = q, delta = rho1 - rho2, rho1 = rho1, rho2 = rho2))

} # cors.to.q()


#' Conversion from a Cohen's q to a correlation difference
#'
#'
#' @param q       Cohen's q effect size.
#' @param rho1    first correlation (either rho1 or rho2 needs to be given)
#' @param rho2    second correlation (either rho1 or rho2 needs to be given)
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#'                if \code{2} a more detailed output is given (plus key
#'                parameters and definitions), if \code{0} no output is printed
#'                on the console.
#'
#' @return
#'   \item{q}{Cohen's q effect size.}
#'   \item{delta}{correlation difference: rho1 - rho2.}
#'   \item{rho1}{first correlation.}
#'   \item{rho2}{second correlation.}
#'
#' @examples
#' q.to.cors(q = 0.10, rho1 = 0.50)
#' q.to.cors(q = 0.30, rho1 = 0.50)
#' q.to.cors(q = 0.50, rho1 = 0.50)
#'
#' @export q.to.cors
q.to.cors <- function(q, rho1 = NULL, rho2 = NULL, verbose = 1) {

  check.numeric(q)
  if (!is.null(rho1)) check.correlation(rho1)
  if (!is.null(rho2)) check.correlation(rho2)
  verbose <- ensure.verbose(verbose)

  if (is.null(rho1) && is.null(rho2))
    stop("Both `rho1` and `rho2` can not be NULL.", call. = FALSE)
  if (!is.null(rho1) && !is.null(rho2))
    stop("Exactly one of the `rho1` or `rho2` should be NULL.", call. = FALSE)

  if (is.null(rho1) && !is.null(rho2)) {
    z1 <- cor.to.z(rho2, FALSE)$z - q
    rho1 <- (exp(2 * z1) - 1) / (exp(2 * z1) + 1)
  } else if (!is.null(rho1) && is.null(rho2)) {
    z2 <- cor.to.z(rho1, FALSE)$z + q
    rho2 <- (exp(2 * z2) - 1) / (exp(2 * z2) + 1)
  }

  if (verbose > 0)
    print(c(q = q, delta = rho1 - rho2, rho1 = rho1, rho2 = rho2))

  invisible(list(q = q, delta = rho1 - rho2, rho1 = rho1, rho2 = rho2))

} # q.to.cors()


#' Conversion from Cohen's d to Common Language Effect Size
#'
#'
#' @description
#' Helper function to convert Cohen's d to common language effect size (or vice
#' versa). The result is the probability of superiority for independent
#' samples. It can be interpreted as the probability that a randomly selected
#' observation from Group 1 exceeds a randomly selected observation from Group
#' 2. The rationale is the same for paired-samples and one-sample designs, but
#' the interpretation differs: For paired samples, it can be interpreted as the
#' probability that the difference score (i.e., the score under Condition 1
#' minus the score under Condition 2) is greater than zero for a randomly
#' selected individual. For a one-sample design, it can be interpreted as the
#' probability that a randomly selected observation is greater than the
#' reference value (e.g., 0).
#'
#' @aliases d.to.cles cles.to.d
#'
#'
#' @param d       Cohen's d
#' @param design  character; one of the "independent", "paired", or
#'                "one.sample". The default is "independent".
#' @param cles    common language effect size.
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#'                if \code{2} a more detailed output is given (plus key
#'                parameters and definitions), if \code{0} no output is printed
#'                on the console.
#'
#' @return
#'   \item{d}{Cohen's d}
#'   \item{cles}{common language effect size.}
#'
#' @examples
#' d.to.cles(0.20) # small
#' d.to.cles(0.50) # medium
#' d.to.cles(0.80) # large
#'
#' cles.to.d(0.5562315)
#' cles.to.d(0.6381632)
#' cles.to.d(0.7141962)
#'
#' @export d.to.cles
d.to.cles <- function(d, design = c("independent", "paired", "one.sample"), verbose = 1) {

  check.numeric(d)
  design <- match.arg(design)
  verbose <- ensure.verbose(verbose)

  prob <- stats::pnorm(d / sqrt(ifelse(design == "independent", 2, 1)))

  if (verbose > 0)
    print(c(cles = prob, d = d))

  invisible(list(cles = prob, d = d))

} # d.to.cles

#' @rdname d.to.cles
#' @export cles.to.d
cles.to.d <- function(cles, design = c("independent", "paired", "one.sample"), verbose = 1) {

  check.proportion(cles)
  design <- match.arg(design)
  verbose <- ensure.verbose(verbose)

  d <- sqrt(ifelse(design == "independent", 2, 1)) * stats::qnorm(cles)

  if (verbose > 0)
    print(c(d = d, cles = cles))

  invisible(list(d = d, cles = cles))

} # cles.to.d


#' Conversion from Means and Standard Deviations to Cohen's d
#'
#' @description
#' Helper function to convert means and standard deviations to Cohen's d.
#'
#'
#' @param mu1 mean of the first group.
#' @param mu2 mean of the second group.
#' @param sd1 standard deviation of the first group.
#' @param sd2 standard deviation of the second group.
#' @param n.ratio \code{n1 / n2} ratio (applies to independent samples only).
#' @param paired if \code{TRUE} paired samples
#' @param rho.paired correlation between repeated measures for paired samples
#' (e.g., pretest and post-test).
#' @param n2 integer; sample size in the second group (or for the single group
#' in paired samples).
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#' if \code{2} a more detailed output is given (plus key parameters and
#' definitions), if \code{0} no output is printed on the console.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{d}{Cohen's d}
#'   \item{pooled.sd}{Pooled standard deviation}
#'   \item{var.ratio}{Ratio of the variance in the two groups (applies to independent samples only)}
#'   \item{n1}{Sample size group 1 (applies to independent samples only)}
#'   \item{n2}{Sample size group 2}
#'
#' @examples
#'
#'
#' # means and standard deviations from independent samples
#' means.to.d(mu1 = 20, mu2 = 17.5,
#'            sd1 = 5, sd2 = 15,
#'            n2 = 30, n.ratio = 1)
#'
#' # means and standard deviations from paired samples
#' means.to.d(mu1 = 20, mu2 = 17.5,
#'            sd1 = 5, sd2 = 15,
#'            n2 = 30, n.ratio = 1,
#'            paired = TRUE,
#'            rho.paired = 0.50)
#'
#' @export means.to.d
means.to.d <- function(mu1, mu2 = 0, sd1 = 1, sd2 = 1, n.ratio = 1, n2 = 1e10,
                       paired = FALSE, rho.paired = 0.50, verbose = 1) {

  func.parms <- as.list(environment())

  check.logical(paired)
  check.numeric(mu1, mu2)
  check.correlation(rho.paired)
  check.positive(sd1, sd2, n.ratio)
  check.sample.size(n2)
  verbose <- ensure.verbose(verbose)

  n1 <- n.ratio * n2

  alpha <- 0.05
  var.ratio <- sd1 ^ 2 / sd2 ^ 2
  f.crit <- stats::qf(alpha, df1 = n1 - 1, df2 = n2 - 1, lower.tail = TRUE)
  if (var.ratio <= f.crit || var.ratio >= (1 / f.crit))
    warning("Interpretation of Cohen's d may no longer be valid when variances differ beyond sampling error.\n", call. = FALSE)

  if (paired) {
    pooled.sd <- sqrt(sd1 ^ 2 + sd2 ^ 2 - 2 * sd1 * sd2 * rho.paired)
  } else {
    pooled.sd <- sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2) / (n1 + n2 - 2))
  }

  d <- (mu1 - mu2) / pooled.sd

  if (verbose > 0)
    print(c(d = d))

  invisible(list(parms = func.parms, d = d, pooled.sd = pooled.sd, var.ratio = var.ratio, n1 = n1, n2 = n2))

} # means.to.d


#' Conversion from Probability Difference to Cohen's h
#'
#' @description
#' Helper function to convert probability difference to Cohen's h (and vice
#' versa).
#'
#'
#' @param prob1   Probability of success in the first group, or under the
#'                alternative hypothesis in the one-sample case).
#' @param prob2   Probability of success in the second group, or under the null
#'                hypothesis in the one-sample case).
#' @param h       Cohen's h effect size.
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#'                if \code{2} a more detailed output is given (plus key
#'                parameters and definitions), if \code{0} no output is printed
#'                on the console.
#'
#' @return
#'   \item{h}{Cohen's h effect size.}
#'   \item{prob1}{probability of success in the first group, or under the
#'                alternative hypothesis in the one-sample case).}
#'   \item{prob2}{probability of success in the second group, or under the
#'                null hypothesis in the one-sample case).}
#'
#' @examples
#' probs.to.h(prob1 = 0.56, prob2 = 0.50)
#'
#' @export probs.to.h
probs.to.h <- function(prob1, prob2 = 0.50, verbose = 1) {

  check.proportion(prob1, prob2)
  verbose <- ensure.verbose(verbose)

  h <- 2 * asin(sqrt(prob1)) - 2 * asin(sqrt(prob2))

  if (verbose > 0)
    print(c(h = h, prob1 = prob1, prob2 = prob2))

  invisible(list(h = h, prob1 = prob1, prob2 = prob2))

} # probs.to.h


#' Conversion from joint probabilities to marginal probabilities for the
#' McNemar test applied to paired binary data.
#'
#'
#' @param prob1   (marginal) probability of success in case group (or after).
#' @param prob2   (marginal) probability of success in matched-control group
#'                (or before).
#' @param rho     the correlation between case and matched-control, or after
#'                and before (phi coefficient).
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#'                if \code{2} a more detailed output is given (plus key
#'                parameters and definitions), if \code{0} no output is printed
#'                on the console.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{prob11}{(joint) probability of success in both groups. 'prob11' and
#'                 'prob00' are known as concordant probs.}
#'   \item{prob10}{(joint) probability of success in case (or after) but
#'                 failure in matched control (or before). 'prob10' and
#'                 'prob01' are known as discordant probs.}
#'   \item{prob01}{(joint) probability of failure in case (or after) but
#'                 success in matched control (or before). prob10' and 'prob01'
#'                 are known as discordant probs.}
#'   \item{prob00}{(joint) probability of failure in both groups. 'prob11' and
#'                 'prob00' are known as concordant probs.}
#'
#' @references
#'   Zhang, S., Cao, J., and Ahn, C. (2017). Inference and sample size
#'   calculation for clinical trials with incomplete observations of paired
#'   binary outcomes. *Statistics in Medicine, 36*(4), 581-591.
#'   https://doi.org/10.1002/sim.7168
#'
#' @examples
#' # example data for a matched case-control design
#' # subject  case     control
#' # <int>    <dbl>    <dbl>
#' #   1        1        1
#' #   2        0        1
#' #   3        1        0
#' #   4        0        1
#' #   5        1        1
#' #   ...     ...      ...
#' #   100      0        0
#'
#' # example summary stats
#' # prob1 = mean(case) which is 0.55
#' # prob2 = mean(control) which is 0.45
#' # rho = cor(case, control) which is 0.4141414
#'
#'
#' # example data for a before-after design
#' # subject  before   after
#' # <int>    <dbl>    <dbl>
#' #   1        1        1
#' #   2        0        1
#' #   3        1        0
#' #   4        0        1
#' #   5        1        1
#' #   ...     ...      ...
#' #   100      0        0
#'
#' # example summary stats
#' # prob1 = mean(after) which is 0.55
#' # prob2 = mean(before) which is 0.45
#' # rho = cor(after, before) which is 0.4141414
#'
#' # convert to a 2 x 2 frequency table
#' freqs <- matrix(c(30, 10, 20, 40), nrow = 2, ncol = 2)
#' colnames(freqs) <- c("control_1", "control_0")
#' rownames(freqs) <- c("case_1", "case_0")
#' freqs
#'
#' # convert to a 2 x 2 proportion table
#' props <- freqs / sum(freqs)
#' props
#'
#' # discordant pairs (0 and 1, or 1 and 0) in 'props' matrix
#' # are the sample estimates of prob01 and prob10
#'
#'
#' # we may not have 2 x 2 joint probs
#' # convert marginal probs to joint probs using summary stats
#' jp <- joint.probs.2x2(prob1 = 0.55, # mean of case (or after)
#'                           prob2 = 0.45, # mean of matched control (or before)
#'                           # correlation b/w matched case-control / before-after
#'                           rho = 0.4141414)
#'
#' # required sample size for exact test
#' # assuming prob01 and prob10 are population parameters
#' power.exact.mcnemar(prob01 = jp$prob01,
#'                     prob10 = jp$prob10,
#'                     power = 0.80, alpha = 0.05,
#'                     method = "exact")
#'
#' # convert joint probs to marginal probs and calc phi coefficient (rho)
#' # these values can be used in other procedures
#' marginal.probs.2x2(prob11 = 0.35, # mean of case (or after)
#'                     prob10 = 0.20, # mean of matched control (or before)
#'                     prob01 = 0.10,
#'                     prob00 = 0.35)
#'
#' @export joint.probs.2x2
joint.probs.2x2 <- function(prob1, prob2, rho = 0.50, verbose = 1) {

  func.parms <- as.list(environment())
  verbose <- ensure.verbose(verbose)

  check.proportion(prob1, prob2)
  check.correlation(rho)

  rho.min <- max(
    -sqrt(prob1 * prob2 / ((1 - prob1) * (1 - prob2))),
    -sqrt((1 - prob1) * (1 - prob2) / (prob2 * prob2))
  )

  rho.max <- min(
    sqrt(prob1 * (1 - prob2) / (prob2 * (1 - prob1))),
    sqrt(prob2 * (1 - prob1) / (prob1 * (1 - prob2)))
  )

  if (rho < rho.min || rho > rho.max) {
    stop("Combination of `prob1`, `prob2` and `rho` is not feasible.\n`rho` should be between ",
         round(rho.min, 3), " and ", round(rho.max, 3), call. = FALSE)
  }

  prob11 <- rho * sqrt(prob1 * (1 - prob1) * prob2 * (1 - prob2)) + prob1 * prob2
  prob10 <- prob1 - prob11
  prob01 <- prob2 - prob11
  prob00 <- 1 - (prob11 + prob10 + prob01)

  if (verbose > 0)
    print(c(rho.min = rho.min, rho.max = rho.max, prob11 = prob11, prob10 = prob10, prob01 = prob01, prob00 = prob00))

  invisible(list(parms = func.parms, rho.min = rho.min, rho.max = rho.max,
                 prob11 = prob11, prob10 = prob10, prob01 = prob01, prob00 = prob00))

} # joint.probs.2x2


# no need to export or document - internal function
# find limits of prob2 given prob1 and rho (or prob1 given prob2 and rho)
prob.limits.paired <- function(prob1 = NULL, prob2 = NULL, rho = 0.50,
                               grid.min = 0.0001, grid.max = 0.9999, step = 0.0001) {

  rho.limits <- function(prob1, prob2) {

    rho.min <- max(
      -sqrt(prob1 * prob2 / ((1 - prob1) * (1 - prob2))),
      -sqrt((1 - prob1) * (1 - prob2) / (prob2 * prob2))
    )

    rho.max <- min(
      sqrt(prob1 * (1 - prob2) / (prob2 * (1 - prob1))),
      sqrt(prob2 * (1 - prob1) / (prob1 * (1 - prob2)))
    )

    c(rho.min, rho.max)

  } # rho.limits

  rho.is.feasible <- function(prob1, prob2, rho) {

    limits <- rho.limits(prob1, prob2)
    rho >= limits[1] & rho <= limits[2]

  } # rho.is.feasible

  if (is.null(prob2) && !is.null(prob1)) {

    feasible <- function(x) rho.is.feasible(prob1 = prob1, prob2 = x, rho = rho)

  } else if (is.null(prob1) && !is.null(prob2)) {

    feasible <- function(x) rho.is.feasible(prob1 = x, prob2 = prob2, rho = rho)

  } else {

    stop("Exactly one of the parameters `prob1` or `prob2` must be given, one has to be NULL.", call. = FALSE)

  }

  # find the feasible region
  grid <- seq(grid.min, grid.max, by = step)
  good <- vapply(grid, feasible, logical(1))

  if (!any(good)) {
    stop("No feasible values found. Check that rho is achievable given the fixed probability.", call. = FALSE)
  }

  feasible.vals <- grid[good]

  c(prob.min = min(feasible.vals), prob.max = max(feasible.vals))

} # prob.limits.paired


#' Conversion from marginal probabilities to joint probabilities for the
#' McNemar test applied to paired binary data.
#'
#'
#' @param prob11  (joint) probability of success in both groups. 'prob11' and
#'                'prob00' are known as concordant probs.
#' @param prob10  (joint) probability of success in case (or after) but failure
#'                in matched control (or before). 'prob10' and 'prob01' are
#'                known as discordant probs.
#' @param prob01  (joint) probability of failure in case (or after) but success
#'                in matched control (or before). prob10' and 'prob01' are
#'                known as discordant probs.
#' @param prob00  (joint) probability of failure in both groups. 'prob11' and
#'                'prob00' are known as concordant probs.
#' @param verbose \code{1} by default (returns test, hypotheses, and results),
#'                if \code{2} a more detailed output is given (plus key
#'                parameters and definitions), if \code{0} no output is printed
#'                on the console.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{prob1}{(marginal) probability of success in case group (or after).}
#'   \item{prob2}{(marginal) probability of success in matched-control group
#'                (or before).}
#'   \item{rho}{the correlation between case and matched-control, or after and
#'              before (phi coefficient).}
#'
#' @references
#'   Zhang, S., Cao, J., and Ahn, C. (2017). Inference and sample size
#'   calculation for clinical trials with incomplete observations of paired
#'   binary outcomes. *Statistics in Medicine, 36*(4), 581-591.
#'   https://doi.org/10.1002/sim.7168
#'
#' @examples
#' # example data for a matched case-control design
#' # subject  case     control
#' # <int>    <dbl>    <dbl>
#' #   1        1        1
#' #   2        0        1
#' #   3        1        0
#' #   4        0        1
#' #   5        1        1
#' #   ...     ...      ...
#' #   100      0        0
#'
#' # example summary stats
#' # prob1 = mean(case) which is 0.55
#' # prob2 = mean(control) which is 0.45
#' # rho = cor(case, control) which is 0.4141414
#'
#'
#' # example data for a before-after design
#' # subject  before   after
#' # <int>    <dbl>    <dbl>
#' #   1        1        1
#' #   2        0        1
#' #   3        1        0
#' #   4        0        1
#' #   5        1        1
#' #   ...     ...      ...
#' #   100      0        0
#'
#' # example summary stats
#' # prob1 = mean(after) which is 0.55
#' # prob2 = mean(before) which is 0.45
#' # rho = cor(after, before) which is 0.4141414
#'
#' # convert to a 2 x 2 frequency table
#' freqs <- matrix(c(30, 10, 20, 40), nrow = 2, ncol = 2)
#' colnames(freqs) <- c("control_1", "control_0")
#' rownames(freqs) <- c("case_1", "case_0")
#' freqs
#'
#' # convert to a 2 x 2 proportion table
#' props <- freqs / sum(freqs)
#' props
#'
#' # discordant pairs (0 and 1, or 1 and 0) in 'props' matrix
#' # are the sample estimates of prob01 and prob10
#'
#'
#' # we may not have 2 x 2 joint probs
#' # convert marginal probs to joint probs using summary stats
#' jp <- joint.probs.2x2(prob1 = 0.55, # mean of case (or after)
#'                           prob2 = 0.45, # mean of matched control (or before)
#'                           # correlation b/w matched case-control / before-after
#'                           rho = 0.4141414)
#'
#' # required sample size for exact test
#' # assuming prob01 and prob10 are population parameters
#' power.exact.mcnemar(prob01 = jp$prob01,
#'                     prob10 = jp$prob10,
#'                     power = 0.80, alpha = 0.05,
#'                     method = "exact")
#'
#' # convert joint probs to marginal probs and calc phi coefficient (rho)
#' # these values can be used in other procedures
#' marginal.probs.2x2(prob11 = 0.35, # mean of case (or after)
#'                     prob10 = 0.20, # mean of matched control (or before)
#'                     prob01 = 0.10,
#'                     prob00 = 0.35)
#'
#' @export marginal.probs.2x2
marginal.probs.2x2 <- function(prob11, prob10, prob01, prob00, verbose = 1) {

  func.parms <- as.list(environment())

  check.proportion(prob11, prob10, prob01, prob00)
  verbose <- ensure.verbose(verbose)

  total <- prob11 + prob10 + prob01 + prob00

  if (abs(total - 1) > 1e-6)
    stop("Joint probabilities must sum to 1.", call. = FALSE)

  prob1 <- prob11 + prob10
  prob2 <- prob11 + prob01

  if (prob1 == 1 || prob1 == 0 || prob2 == 1 || prob2 == 0) {

    rho <- NA
    warning("Undefined correlation: division by zero in denominator.", call. = FALSE)

  } else {

    rho <- (prob11 - prob1 * prob2) / sqrt(prob1 * (1 - prob1) * prob2 * (1 - prob2))

  }

  if (verbose > 0)
    print(c(prob1 = prob1, prob2 = prob2, rho = rho))

  invisible(list(parms = func.parms, prob1 = prob1, prob2 = prob2, rho = rho))

} # marginal.probs.2x2


# internal function to get some chisq stat

#' Conversion from Probabilities to Cohen's w
#'
#' @description
#' Helper function to convert (multinomial or product-multinomial)
#' probabilities to Cohen's w.
#'
#'
#' @param prob.matrix      a vector or matrix of cell probabilities under
#'                         the alternative hypothesis
#' @param null.prob.matrix a vector or matrix of cell probabilities under the
#'                         null hypothesis. Calculated automatically when
#'                         \code{prob.matrix} is specified. The default can be
#'                         overwritten by the user via providing a vector of
#'                         the same size or matrix of the same dimensions as
#'                         \code{prob.matrix}.
#' @param verbose          \code{1} by default (returns test, hypotheses, and
#'                         results), if \code{2} a more detailed output is
#'                         given (plus key parameters and definitions), if
#'                         \code{0} no output is printed on the console.
#'
#' @return
#'    \item{w}{Cohen's w effect size. It can be any of Cohen's W, Phi
#'             coefficient, Cramer's V. Phi coefficient is defined as
#'             \code{sqrt(X2 / n)} and Cramer's V is defined as
#'             \code{sqrt(X2 / (n * v))} where \code{v} is
#'             \code{min(nrow - 1, ncol - 1)} and X2 is the chi-square
#'             statistic.}
#' \item{df}{degrees of freedom.}
#'
#' @references
#'   Cohen, J. (1988). *Statistical power analysis for the behavioral sciences*
#'   (2nd ed.). Lawrence Erlbaum Associates.
#'
#' @examples
#'   # ---------------------------------------------------------#
#'   # Example 1: Cohen's W                                     #
#'   # goodness-of-fit test for 1 x k or k x 1 table            #
#'   # How many subjects are needed to claim that               #
#'   # girls choose STEM related majors less than males?       #
#'   # ---------------------------------------------------------#
#'
#'   ## from https://www.aauw.org/resources/research/the-stem-gap/
#'   ## 28 percent of the  workforce in STEM field is women
#'   prob.vector <- c(0.28, 0.72)
#'   null.prob.vector <- c(0.50, 0.50)
#'   probs.to.w(prob.vector, null.prob.vector)
#'
#'   power.chisq.gof(w = 0.44, df = 1,
#'                   alpha = 0.05, power = 0.80)
#'
#'
#'   # ---------------------------------------------------------#
#'   # Example 2: Phi Coefficient (or Cramer's V or Cohen's W)  #
#'   # test of independence for 2 x 2 contingency tables        #
#'   # How many subjects are needed to claim that               #
#'   # girls are underdiagnosed with ADHD?                      #
#'   # ---------------------------------------------------------#
#'
#'   ## from https://web.archive.org/web/20250601054306/https://time.com/growing-up-with-adhd/
#'   ## 5.6 percent of girls and 13.2 percent of boys are diagnosed with ADHD
#'   prob.matrix <- rbind(c(0.056, 0.132),
#'                        c(0.944, 0.868))
#'   colnames(prob.matrix) <- c("Girl", "Boy")
#'   rownames(prob.matrix) <- c("ADHD", "No ADHD")
#'   prob.matrix
#'
#'   probs.to.w(prob.matrix)
#'
#'   power.chisq.gof(w = 0.1302134, df = 1,
#'                   alpha = 0.05, power = 0.80)
#'
#'
#'   # --------------------------------------------------------#
#'   # Example 3: Cramer's V (or Cohen's W)                    #
#'   # test of independence for j x k contingency tables       #
#'   # How many subjects are needed to detect the relationship #
#'   # between depression severity and gender?                 #
#'   # --------------------------------------------------------#
#'
#'   ## from https://doi.org/10.1016/j.jad.2019.11.121
#'   prob.matrix <- cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078),
#'                        c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101))
#'   rownames(prob.matrix) <- c("Normal", "Mild", "Moderate",
#'                              "Severe", "Extremely Severe")
#'   colnames(prob.matrix) <- c("Female", "Male")
#'   prob.matrix
#'
#'   probs.to.w(prob.matrix)
#'
#'   power.chisq.gof(w = 0.03022008, df = 4,
#'                   alpha = 0.05, power = 0.80)
#'
#' @export probs.to.w
probs.to.w <- function(prob.matrix, null.prob.matrix = NULL, verbose = 1) {

  verbose <- ensure.verbose(verbose)

  if (any(prob.matrix < 0) || any(prob.matrix > 1))
    stop("Matrix elements outside of [0, 1] range.", call. = FALSE)

  if (is.null(null.prob.matrix)) {
    if (is.vector(prob.matrix))
      null.prob.matrix <- rep(1 / length(prob.matrix), length(prob.matrix))
    else
      null.prob.matrix <- outer(rowSums(prob.matrix), colSums(prob.matrix)) / sum(prob.matrix)
  } else {
    if (any(null.prob.matrix < 0) || any(null.prob.matrix > 1))
      stop("Matrix elements outside of [0, 1] range.", call. = FALSE)
  }

  if (is.vector(prob.matrix)) {
    if (length(prob.matrix) != length(null.prob.matrix))
      stop("Length of `prob.matrix` and `null.prob.matrix` should match.", call. = FALSE)
    if (sum(prob.matrix) != 1 || sum(null.prob.matrix) != 1)
      stop("Cell probabilities should sum to 1.", call. = FALSE)
  } else if (is.matrix(prob.matrix)) {
    if (any(dim(prob.matrix) != dim(null.prob.matrix)))
      stop("Dimensions for `prob.matrix` and `null.prob.matrix` do not match.", call. = FALSE)
  } else {
    stop("`prob.matrix` must be either a vector or a matrix.", call. = FALSE)
  }

  df <- ifelse(is.vector(prob.matrix), length(prob.matrix) - 1, (nrow(prob.matrix) - 1) * (ncol(prob.matrix) - 1))
  chisq <- sum((prob.matrix - null.prob.matrix) ^ 2 / null.prob.matrix)
  mdf <- ifelse(is.vector(prob.matrix), 1, min(nrow(prob.matrix) - 1, ncol(prob.matrix) - 1))
  w <- sqrt(chisq / (sum(prob.matrix) * mdf))

  if (w > 1)
    warning("w > 1 is unrealistic, please check your input.", call. = FALSE)

  if (verbose > 0)
    print(c(w = w, df = df))

  invisible(list(w = w, df = df, prob.matrix = prob.matrix, null.prob.matrix = null.prob.matrix))

} # probs.to.w


#' Conversion from Means and Standard Deviations to Cohen's f and Eta-squared
#'
#' @description
#' Calculates Cohen's f or Eta-squared for one-way ANOVA/ANCOVA. Set \code{k.cov =
#' 0} for one-way ANOVA (without any pretest or covariate adjustment). Set
#' \code{k.cov > 0} in combination with \code{r.squared > 0} for one-way ANCOVA (with
#' pretest or covariate adjustment).
#'
#'
#' @details
#' Note that R has a partial matching feature which allows you to specify
#' shortened versions of arguments, such as \code{mu} or \code{mu.vec} instead
#' of \code{mu.vector}, or such as \code{k} or \code{k.cov} instead of
#' \code{k.covariates}.
#'
#' @param mu.vector     vector of adjusted means (or estimated marginal means)
#'                      for each level of a factor.
#' @param sd.vector     vector of unadjusted standard deviations for each level
#'                      of a factor.
#' @param n.vector      vector of sample sizes for each level of a factor.
#' @param factor.levels integer; number of levels or groups in each factor. For
#'                      example, for two factors each having two levels or
#'                      groups use e.g. c(2, 2), for three factors each having
#'                      two levels or groups use e.g. c(2, 2, 2)
#' @param r.squared     explanatory power of covariates (R-squared) in the
#'                      ANCOVA model. The default is \code{r.squared = 0},
#'                      which means an ANOVA model would be of interest.
#' @param k.covariates  integer; number of covariates in the ANCOVA model. The
#'                      default is \code{k.covariates = 0}, which means an
#'                      ANOVA model would be of interest.
#' @param verbose       \code{1} by default (returns results), if \code{0} no
#'                      output is printed on the console.
#'
#' @return
#'   \item{f}{Cohen's f}
#'   \item{eta.squared}{(partial) eta-squared.}
#'   \item{df1}{numerator degrees of freedom.}
#'   \item{df2}{denominator degrees of freedom.}
#'   \item{ncp}{non-centrality parameter under alternative.}
#'
#' @references
#'   Keppel, G., & Wickens, T. D. (2004). Design and analysis: A researcher's
#'   handbook (4th ed.). Pearson.
#'
#' @examples
#' means.to.etasq(mu.vector = c(0.50, 0), # marginal means
#'                sd.vector = c(1, 1), # unadjusted standard deviation
#'                n.vector = c(33, 33), # sample size (will be calculated)
#'                k.cov = 1, # number of covariates
#'                r.squared = 0.50)
#'
#' @export means.to.etasq
means.to.etasq <- function(mu.vector, sd.vector, n.vector, k.covariates = 0, r.squared = 0, factor.levels = NULL, verbose = 1) {

  if (!is.vector(mu.vector) || !is.numeric(mu.vector))
    stop("Provide a vector of means (`mu.vector`) with its length equal to number of groups.", call. = FALSE)
  check.vector(mu.vector, check.numeric)
  check.vector(sd.vector, check.positive)
  if (!is.null(n.vector)) check.vector(n.vector, check.sample.size)
  check.same.lengths(mu.vector, sd.vector, n.vector)
  if (is.null(factor.levels)) factor.levels <- length(mu.vector)
  if (length(factor.levels) > 1)
    stop("Factorial designs are not allowed.", call. = FALSE)
  check.factor.level(factor.levels)
  if (length(mu.vector) != factor.levels)
    stop("Length of the vector of means (`mu.vector`) does not match number of levels.", call. = FALSE)
  if (r.squared > 1 || r.squared < 0 || !is.numeric(r.squared) || length(r.squared) != 1)
    stop("R-squared (explanatory power of covariates) takes a value between 0 and 1.", call. = FALSE)
  check.nonnegative(k.covariates)
  if (r.squared > 0 && k.covariates < 1)
    stop("Explanatory power of covariates is expected to be non-zero when number of covariates is non-zero.", call. = FALSE)

  n.total <- sum(n.vector)
  mu_bar <- sum(n.vector * mu.vector) / n.total

  sigma2_pooled <- sum((n.vector - 1) * sd.vector ^ 2) / (n.total - length(mu.vector))

  sigma2_between <- sum(n.vector * (mu.vector - mu_bar) ^ 2) / n.total
  sigma2_error <- sigma2_pooled * (1 - r.squared)

  f.squared <- sigma2_between / sigma2_error
  eta.squared <- sigma2_between / (sigma2_between + sigma2_error)

  u <- prod(factor.levels - 1)
  v <- n.total - length(mu.vector) - k.covariates
  lambda <- f.squared * n.total

  if (verbose > 0)
    print(c(f = sqrt(f.squared), eta.squared = eta.squared,
            df1 = u, df2 = v, ncp = lambda))

  invisible(list(f = sqrt(f.squared), eta.squared = eta.squared,
                 df1 = u, df2 = v, ncp = lambda))

} # means.to.etasq
