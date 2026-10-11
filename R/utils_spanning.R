#' Compute Lexicographic Rank of Last Element in a Vector
#'
#' Computes the rank of the last element in a numeric vector `x` relative to the other elements,
#' using a second vector `uu` for tie-breaking. The function counts how many elements the last entry
#' of `x` is greater than (or tied with but ranked higher based on `uu`), and returns its 1-based rank.
#'
#' @param x A numeric or integer vector. The last element of `x` is the one to rank against the rest.
#' @param uu A numeric or integer vector of the same length as `x`, used to break ties in `x`.
#'
#' @return An integer representing the 1-based lexicographic rank of the last element in `x`.
#' If `x` is of length 0 or 1, returns `NA_integer_`.
#'
#' @keywords internal
#'
#' @noRd
#'
f_ranklex <- function(x, uu) {
  n <- length(x)  # Get the length of vector x

  if (n == 0 || n == 1) {
        return(NA_integer_)
  }

  x_last <- x[n]      # The last element of x (the one to rank)
  uu_last <- uu[n]    # The last element of uu (used for tie-breaking)
  x_others <- x[-n]   # All elements of x except the last
  uu_others <- uu[-n] # All elements of uu except the last

  out <- sum((x_last > x_others) | ((x_last == x_others) & (uu_last > uu_others))) + 1
  return(out)
}


# Cauchy combination of p-values, valid under arbitrary dependence.
#
# A p-value equal to one, a null event for continuous data, is replaced by 1 - 1e-12, as in
# the paper (Section 3.1): otherwise tan(-pi/2) = -1.63e16 outweighs any small p-value.
# p = 0 gives tan(pi/2) = 1.63e16, finite (pi/2 is not exactly representable), so the
# combination saturates rather than overflowing.
#
# Empty input is reachable: callers pass na.omit(<per-asset p-values>), which is
# length zero when every asset is missing. mean(numeric(0)) is NaN, so without
# the guard the function would return a silent NaN where NA is meant.
f_cauchypv <- function(p) {
  if (!length(p)) return(NA_real_)
  p[p >= 1] <- 1 - 1e-12
  out <- 0.5 - atan(mean(tan((0.5 - p) * pi))) / pi
  return(out)
}

#' Maximum-Likelihood Covariance Matrix
#'
#' The covariance with divisor \eqn{T} rather than \eqn{T-1}, as in the
#' likelihood-ratio forms of the spanning tests.
#'
#' @param X Numeric \eqn{T \times p} matrix.
#'
#' @return The \eqn{p \times p} matrix \eqn{T^{-1} \sum_t (x_t - \bar x)(x_t - \bar x)'}.
#'
#' @keywords internal
#'
#' @noRd
#'
f_cov_ml <- function(X) {
  Xc <- scale(X, center = TRUE, scale = FALSE)
  out <- crossprod(Xc) / nrow(X)
  return(out)
}

#' Multiplier Weights of the SCT, One Draw per Statistic
#'
#' Draws the weights \eqn{\kappa_{i,t} = \prod_{l=1}^{k} \kappa_{l,i,t}}, with
#' \eqn{\kappa_{l,i,t} \sim N(1, 1)} independent across factors \eqn{l}, dates
#' \eqn{t} and statistics \eqn{i}: one \eqn{T \times N} matrix for the alpha scores
#' and an independent one for the delta scores.
#'
#' @param Tn Number of observations \eqn{T}.
#' @param N Number of test assets.
#' @param k Number of multiplier factors \eqn{L}. If \code{k <= 0}, the weights are 1.
#' @param cseed Seed of the draw (default 123). The alpha factors are drawn first,
#'   then the delta factors.
#'
#' @return A list with the \eqn{T \times N} matrices \code{A} (alpha scores) and
#'   \code{D} (delta scores); both are the scalar 1 if \code{k <= 0}.
#'
#' @keywords internal
#'
#' @noRd
#'
f_mult <- function(Tn, N, k, cseed = 123) {

  if (k <= 0) return(list(A = 1, D = 1))
  # Seed locally for reproducibility, but restore the caller's RNG state on
  # exit so this helper never leaves a side effect on the global stream.
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }
  set.seed(cseed)
  f_draw <- function() {
    w <- matrix(1, Tn, N)
    for (l in seq_len(k)) w <- w * matrix(rnorm(Tn * N, mean = 1, sd = 1), Tn, N)
    return(w)
  }
  wa <- f_draw()
  wd <- f_draw()
  out <- list(A = wa, D = wd)
  return(out)
}

#' Perform Folded T-tests and Normal Approximations on Matrix Data
#'
#' Splits rows of a matrix into \code{k} folds, applies a summary function (default: column means)
#' within each fold, and computes two sets of p-values for each column: one using Student's t-test
#' across folds, and another using a normal approximation based on fold variability.
#'
#' @param eps A numeric matrix or data frame with observations in rows and variables in columns.
#' @param func A function to summarize each fold (default is \code{colMeans}).
#' @param k A scalar controlling the number of folds as \code{floor(nrow(eps)^k)}. Default is \code{1/1.5}.
#'
#' @return A list with two named numeric vectors:
#' \describe{
#'   \item{student}{P-values from a t-test across the k folds for each column.}
#'   \item{normal}{P-values from a normal approximation using the ratio of column means to standard deviation across folds.}
#' }
#'
#' @keywords internal
#'
#' @noRd
#'
f_ttest <- function(eps, func = colMeans, k = 1/1.5) {

  eps <- as.matrix(eps)
  n <- nrow(eps)
  k <- floor(n^k)
  # The t-statistic below has df = k - 1 and divides by it, so a single fold
  # would give 0/0. This needs an absurdly short sample (T < 8 at the default
  # exponent), but fail loudly rather than return NaN.
  if (k < 2L)
    stop("sample too short for the subseries test: it yields ", k,
         " subseries; at least 2 are required.")
  folds <- cut(seq_len(n), k, labels = FALSE)

  # Fold-wise summaries (k x ncol); func defaults to colMeans
  estvar <- matrix(NA_real_, nrow = k, ncol = ncol(eps))
  for (i in seq_len(k)) {
    estvar[i, ] <- func(eps[folds == i, , drop = FALSE])
  }

  # One-sample two-sided t-test p-value per column, computed directly and
  # vectorised across columns (identical to apply(estvar, 2, t.test)$p.value but
  # without the per-column t.test() overhead): t = mean / (sd / sqrt(k)),
  # df = k - 1, with sd the standard deviation of the k fold summaries.
  m     <- colMeans(estvar)
  sdcol <- sqrt(colSums((estvar - rep(m, each = k))^2) / (k - 1))
  tstat <- m / (sdcol / sqrt(k))
  student <- 2 * stats::pt(-abs(tstat), df = k - 1)

  normal <- 2 * pnorm(-abs(colMeans(eps) / sdcol))

  out <- list(student = student,
              normal = normal)

  return(out)
}

#' Combine Folded T-Test P-values Using the Cauchy Method
#'
#' Applies the \code{\link{f_ttest}} function to a matrix of residuals or estimates, and combines the resulting
#' p-values (from Student's t-test and normal approximation) using the Cauchy combination method.
#'
#' @param eps A numeric matrix or data frame with observations in rows and variables in columns.
#' @param ... Additional arguments passed to \code{\link{f_ttest}} (e.g., \code{func}, \code{k}).
#'
#' @return A named numeric vector of length 2 with combined p-values:
#' \describe{
#'   \item{student}{Combined p-value from the t-test-based column p-values.}
#'   \item{normal}{Combined p-value from the normal-approximation column p-values.}
#' }
#'
#' @keywords internal
#'
#' @noRd
#'
f_testbm <- function(eps, ...) {

  val <- f_ttest(eps, ...)
  out <- c(student = f_cauchypv(val$student),
           normal = f_cauchypv(val$normal))
  return(out)
}
