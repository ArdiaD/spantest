#' Huberman-Kandel Spanning Test (1987)
#'
#' Tests the joint null hypothesis \eqn{H_0: \alpha = 0, \delta = 0} using the
#' mean-variance spanning criterion of Huberman and Kandel (1987). This test evaluates
#' whether adding new assets (R2) improves the efficient frontier spanned by a benchmark (R1).
#'
#' @param R1 Numeric matrix of benchmark returns (T x K)
#' @param R2 Numeric matrix of test asset returns (T x N)
#'
#' @return A list with:
#' \describe{
#'   \item{\code{pval}}{The p-value associated with the test statistic.}
#'   \item{\code{stat}}{The test statistic, F-distributed under the null.}
#'   \item{\code{H0}}{The null hypothesis tested: \code{"alpha = 0 and delta = 0"}.}
#' }
#'
#' @references
#' Huberman, G., & Kandel, S. (1987). "Mean-Variance Spanning." \emph{Journal of Finance}, 42(4), 873–888.
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)  # Benchmark returns
#' R2 <- matrix(rnorm(200), 100, 2)  # Test returns
#' result <- span_hk(R1, R2)
#' print(result$stat)  # F-statistic
#' print(result$pval)  # P-value
#' print(result$H0)    # Null hypothesis
#'
#' # Interpretation:
#' # A p-value below 0.05 suggests rejecting the null that the benchmark spans the extended portfolio.
#'
#' @export
span_hk <- function(R1, R2) {
  R <- cbind(R1, R2)
  T <- nrow(R)
  K <- ncol(R1)
  N <- ncol(R2)

  df1 <- 2 * N
  df2 <- 2 * (T - K - N)
  if (df2 < 1) return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0 and delta = 0"))

  mu <- matrix(colMeans(R), ncol = 1)
  mu1 <- matrix(colMeans(R1), ncol = 1)

  Sigma <- cov(R)
  Sigma1 <- cov(R1)
  iSigma <- tryCatch(solve(Sigma), error = function(e) return(NULL))
  iSigma1 <- tryCatch(solve(Sigma1), error = function(e) return(NULL))
  if (is.null(iSigma) || is.null(iSigma1)) return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0 and delta = 0"))

  a <- t(mu) %*% iSigma %*% mu
  b <- t(mu) %*% iSigma %*% rep(1, K + N)
  c <- t(rep(1, K + N)) %*% iSigma %*% rep(1, K + N)
  d <- a * c - b^2

  a1 <- t(mu1) %*% iSigma1 %*% mu1
  b1 <- t(mu1) %*% iSigma1 %*% rep(1, K)
  c1 <- t(rep(1, K)) %*% iSigma1 %*% rep(1, K)
  d1 <- a1 * c1 - b1^2

  U <- (c1 + d1) / (c + d)
  stat <- ((T - K - N) / N) * (1 / sqrt(U) - 1)
  pval <- pf(stat, df1, df2, lower.tail = FALSE)

  list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "alpha = 0 and delta = 0")
}

