#' F2 Spanning Test for Slope Coefficients (Variance Spanning)
#'
#' Implements the F2 test which examines the null hypothesis
#' \eqn{H_0: \delta = 0}, indicating that the test assets do not improve the
#' minimum variance frontier beyond the benchmark assets.
#'
#' The test compares variance-related parameters of the extended portfolio versus
#' the benchmark portfolio. It assumes i.i.d. normally distributed returns and full-rank covariance matrices.
#'
#' @param R1 A numeric matrix of benchmark asset returns (T observations by K assets).
#' @param R2 A numeric matrix of test asset returns (T observations by N assets).
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value associated with the F2 test statistic under the null hypothesis.}
#'   \item{\code{stat}}{The F2 test statistic measuring variance spanning.}
#'   \item{\code{H0}}{A character string describing the null hypothesis tested: \code{"delta = 0"}.}
#' }
#'
#' @references
#' Huberman, G., & Kandel, S. (1987). "Mean-Variance Spanning." \emph{The Journal of Finance}, 42(4), 873–888. \cr
#' Kan, R., & Zhou, G. (2012). "Tests of Mean-Variance Spanning." \emph{Annals of Economics and Finance}, 13(1), 145–193.
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)  # Benchmark asset returns
#' R2 <- matrix(rnorm(200), 100, 2)  # Test asset returns
#' result <- span_f2(R1, R2)
#' result$stat  # F2 statistic
#' result$pval  # P-value of the test
#' result$H0    # Null hypothesis tested
#'
#' @export
span_f2 <- function(R1, R2) {
  # Combine and define dimensions
  R <- cbind(R1, R2)
  K <- ncol(R1)
  N <- ncol(R2)
  T <- nrow(R)

  df1 <- N
  df2 <- T - K - N + 1
  if (df2 < 1) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  # Compute a, b, c, d for full portfolio (R)
  mu  <- matrix(colMeans(R), ncol = 1)
  V   <- cov(R)
  iV  <- solve(V)
  one <- matrix(1, K + N, 1)

  a  <- t(mu) %*% iV %*% mu
  b  <- t(mu) %*% iV %*% one
  c_ <- t(one) %*% iV %*% one
  d  <- a * c_ - b^2

  # Compute a1, b1, c1, d1 for benchmark only (R1)
  mu1  <- matrix(colMeans(R1), ncol = 1)
  V1   <- cov(R1)
  iV1  <- solve(V1)
  one1 <- matrix(1, K, 1)

  a1 <- t(mu1) %*% iV1 %*% mu1
  b1 <- t(mu1) %*% iV1 %*% one1
  c1 <- t(one1) %*% iV1 %*% one1
  d1 <- a1 * c1 - b1^2

  # Compute F2 statistic
  stat <- (df2 / df1) * (((c_ + d) / (c1 + d1)) * ((1 + a1) / (1 + a)) - 1)
  pval <- pf(stat, df1, df2, lower.tail = FALSE)

  return(list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "delta = 0"))
}
