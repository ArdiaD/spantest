#' Britten-Jones Tangency Portfolio Spanning Test
#'
#' Implements the test proposed by Britten-Jones (1999), which evaluates whether the tangency portfolio
#' formed from the combined set of benchmark and test assets is spanned by the benchmark assets alone.
#' The test uses an F-statistic based on a regression of a constant vector on differences in returns,
#' reflecting whether the benchmark assets span the Sharpe-optimal portfolio of the full asset universe.
#'
#' This test is particularly relevant for assessing if the benchmark portfolio captures the maximum
#' Sharpe ratio achievable when adding new assets.
#'
#' @param R1 Numeric matrix of benchmark asset returns with dimensions (T x N).
#' @param R2 Numeric matrix of test asset returns with dimensions (T x K).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{pval}}{The p-value under the null hypothesis that the tangency portfolio
#'   is spanned by the benchmark assets.}
#'   \item{\code{stat}}{The F-statistic of the test.}
#'   \item{\code{H0}}{A character string describing the null hypothesis tested: \code{"Tangency portfolio is spanned by benchmark"}.}
#' }
#'
#' @references
#' Britten-Jones, M. (1999). "The Sampling Error in Estimates of Mean-Variance Efficient Portfolio Weights."
#' \emph{The Journal of Finance}, 54(2), 655–671.
#'
#' @examples
#' set.seed(321)
#' R1 <- matrix(rnorm(300), 100, 3)  # Benchmark asset returns
#' R2 <- matrix(rnorm(200), 100, 2)  # Test asset returns
#' res <- span_bj(R1, R2)
#' res$stat   # F-statistic
#' res$pval   # p-value
#' res$H0     # Null hypothesis tested
#'
#' @export
span_bj <- function(R1, R2) {
  Rbig <- cbind(R1, R2)  # Combine benchmark (R1) and test (R2)
  t <- nrow(Rbig)
  K <- ncol(R1)          # K = benchmark assets
  N <- ncol(R2)          # N = test assets

  if ((t - K - N) < 1) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  y <- rep(1, t)

  # Take differences between first column and others
  Diff <- sweep(Rbig[, -1, drop = FALSE], 1, Rbig[, 1], FUN = function(x, y) y - x)

  X <- Diff
  XtXi <- crossprod(X)
  XtXi_inv <- solve(XtXi)
  coef <- XtXi_inv %*% crossprod(X, y)
  resid <- y - X %*% coef
  sigma2 <- drop(crossprod(resid) / (t - ncol(X)))

  offset <- (K - 1)
  C <- cbind(matrix(0, N, offset), diag(N))
  num <- as.numeric(
    crossprod(C %*% coef, solve(C %*% XtXi_inv %*% t(C)) %*% (C %*% coef)) / N
  )

  stat <- num / sigma2
  pval <- pf(stat, N, t - ncol(X), lower.tail = FALSE)

  list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "delta = 0")
}

