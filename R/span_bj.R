#' Britten–Jones Tangency-Portfolio Spanning Test (1999)
#'
#' Tests whether the tangency (maximum Sharpe) portfolio of the augmented
#' universe (benchmarks + test assets) is spanned by the benchmark assets
#' alone. Following Britten–Jones (1999), the statistic arises from a
#' regression of a constant on return differences and yields an \eqn{F}
#' test of the tangency-spanning restriction.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value for the \eqn{F}-statistic under the null.}
#'   \item{\code{stat}}{Britten–Jones \eqn{F}-statistic.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"tangency portfolio spanned by benchmark"}.}
#' }
#'
#' @details
#' With \code{X} formed from pairwise differences, the reference distribution is
#' \eqn{F_{N,\ T - ncol(X)}}; here \code{ncol(X)} = \eqn{K + N - 1}.
#' Finite-sample feasibility requires \code{T - (K + N - 1) >= 1}.
#'
#' @references
#' \insertRef{BrittenJones1999}{spantest} \cr
#'
#' @examples
#' set.seed(321)
#' R1 <- matrix(rnorm(300), 100, 3)  # benchmarks: T=100, K=3
#' R2 <- matrix(rnorm(200), 100, 2)  # tests:      T=100, N=2
#' out <- span_bj(R1, R2)
#' out$stat; out$pval; out$H0
#'
#' @family Alpha Spanning Tests
#'
#' @importFrom stats pf
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

