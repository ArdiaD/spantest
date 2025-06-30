#' Pesaran-Yamagata Test for Cross-Sectional Dependence of Intercepts
#'
#' Implements the Pesaran-Yamagata (2023) test to assess whether intercepts in a multi-factor
#' spanning regression are jointly zero, while accounting for cross-sectional dependence across assets.
#' This test is used to check if the benchmark assets' returns span the test assets in terms of expected returns.
#'
#' @param R1 Numeric matrix of benchmark returns with dimensions (T x K), where T is the number of time observations and K the number of benchmark assets.
#' @param R2 Numeric matrix of test asset returns with dimensions (T x N), where N is the number of test assets.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value of the test under the null hypothesis.}
#'   \item{\code{stat}}{Test statistic value (standardized).}
#'   \item{\code{H0}}{Description of the null hypothesis: \code{"alpha = 0"}.}
#' }
#'
#' @details
#' The Pesaran-Yamagata test extends classical intercept tests by allowing for cross-sectional
#' dependence among test assets. It uses residual covariance structure to adjust the test statistic,
#' improving inference in panels with correlated assets.
#'
#' The null hypothesis is that all intercepts are zero (\eqn{\alpha = 0}), implying that
#' the benchmark assets fully explain the expected returns of the test assets.
#'
#' @references
#' Pesaran, M. H., & Yamagata, T. (2023). "Testing for alpha in linear factor pricing models with a large number of securities."
#' \emph{Journal of Financial Econometrics}, 22(2), 407-460.
#'
#' @examples
#' set.seed(123)
#' # Simulate benchmark and test asset returns
#' R1 <- matrix(rnorm(300), nrow = 100, ncol = 3)  # Benchmark returns
#' R2 <- matrix(rnorm(200), nrow = 100, ncol = 2)  # Test asset returns
#' res <- span_py(R1, R2)
#' res$pval   # P-value of the test
#' res$stat   # Test statistic
#' res$H0     # Null hypothesis being tested
#'
#' # Interpret result
#' if (res$pval < 0.05) {
#'   cat("Reject null: intercepts differ significantly\n")
#' } else {
#'   cat("Fail to reject null: no evidence intercepts differ\n")
#' }
#'
#' @export
span_py <- function(R1, R2) {
  X <- R1
  Y <- R2
  T <- nrow(X)
  N <- ncol(Y)
  K <- ncol(X)

  if ((T - K - N) < 1 || (T - K - 1) <= 4) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))
  }

  XX <- cbind(1, X)
  Bhat1 <- solve(crossprod(XX)) %*% crossprod(XX, Y)
  Ehat1 <- Y - XX %*% Bhat1
  SigmaU <- crossprod(Ehat1) / T

  v <- T - K - 1
  ones <- matrix(1, T, 1)
  X_crossprod <- crossprod(X)
  X_crossprod_inv <- solve(X_crossprod)
  X_ones <- X %*% X_crossprod_inv %*% crossprod(X, ones)
  MX_ones <- ones - X_ones

  num_scalar <- crossprod(ones, MX_ones)[1, 1] * v
  num <- matrix(num_scalar, N, 1)
  diag_SigmaU <- diag(SigmaU)
  t2 <- (Bhat1[1, ]^2) * num / (T * diag_SigmaU)

  pN <- 0.05 / (N - 1)
  thetaN <- qnorm(1 - pN / 2)^2
  rhobar <- 0
  for (i in 2:N) {
    for (j in 1:(i - 1)) {
      temp <- SigmaU[i, j] / sqrt(SigmaU[i, i] * SigmaU[j, j])
      temp2 <- temp^2
      if (v * temp2 >= thetaN) {
        rhobar <- rhobar + temp2
      }
    }
  }
  rhobar <- rhobar * 2 / (N * (N - 1))

  Jalpha2 <- sum(t2 - v / (v - 2)) / sqrt(N)
  den <- (v / (v - 2)) * sqrt(2 * (v - 1) * (1 + (N - 1) * rhobar) / (v - 4))
  stat <- Jalpha2 / den
  pval <- 1 - pnorm(stat)

  list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "alpha = 0")
}

