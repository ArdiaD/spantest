#' Kempf–Memmel GMVP Spanning Test
#'
#' Tests whether the Global Minimum Variance Portfolio (GMVP) of the combined
#' (benchmark + test) universe equals the GMVP of the benchmark assets alone.
#' Following Kempf & Memmel (2006), the null assesses whether adding new assets
#' improves the minimum-variance frontier.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times N}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times K}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value for the F-statistic under the null.}
#'   \item{\code{stat}}{F-statistic value.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"GMVP(bmk) = GMVP(full)"}.}
#' }
#'
#' @details
#' The null hypothesis \eqn{H_0} is that augmenting the benchmark set with the
#' test assets does not change the GMVP weights (\eqn{\Delta = 0}), i.e.,
#' the GMVP of the full universe coincides with that of the benchmark subset.
#' The test is implemented via a linear restriction on coefficients in an
#' equivalent regression representation, yielding an \eqn{F}-statistic.
#'
#' @references
#' \insertRef{KempfMemmel2006}{spantest} \cr
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)  # benchmarks: T=100, N=3
#' R2 <- matrix(rnorm(200), 100, 2)  # tests:      T=100, K=2
#' ans <- span_km(R1, R2)
#' ans$pval; ans$stat; ans$H0
#'
#' @family Variance Spanning Tests
#'
#' @importFrom stats pf
#' @export
span_km <- function(R1, R2) {

  R <- cbind(R1, R2)
  TT <- nrow(R)
  p <- ncol(R1)
  p2 <- ncol(R2)

  y <- rep(1, TT)

  Diff <- sweep(R[, -1, drop = FALSE], 1, R[, 1], FUN = function(x, y) y - x)

  X <- Diff
  XtX <- crossprod(X)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) return(NULL))
  if (is.null(XtX_inv)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  beta_hat <- XtX_inv %*% crossprod(X, y)
  residuals <- y - X %*% beta_hat
  sigma2 <- drop(crossprod(residuals) / (TT - ncol(X)))

  offset <- p - 1
  C <- cbind(matrix(0, p2, offset), diag(p2))

  middle <- tryCatch(solve(C %*% XtX_inv %*% t(C)), error = function(e) return(NULL))
  if (is.null(middle)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  theta_part <- C %*% beta_hat
  F_stat <- as.numeric(t(theta_part) %*% middle %*% theta_part / (p2 * sigma2))
  p_val <- pf(F_stat, p2, TT - ncol(X), lower.tail = FALSE)

  return(list(pval = p_val, stat = F_stat, H0 = "delta = 0"))
}
