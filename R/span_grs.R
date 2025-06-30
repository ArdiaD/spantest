#' Gibbons-Ross-Shanken (GRS) Spanning Test
#'
#' Implements the Gibbons, Ross, and Shanken (1989) test for the joint null hypothesis that all intercepts
#' (\eqn{\alpha = 0}) are zero in the multivariate regression of test asset excess returns on benchmark portfolios.
#' The test assumes homoskedasticity and normally distributed errors. It is most reliable when the number of time periods
#' is sufficiently greater than the number of test assets and benchmark factors.
#'
#' @param R1 A numeric matrix of benchmark returns (T x K).
#' @param R2 A numeric matrix of test asset returns (T x N).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{pval}}{The associated p-value under the null.}
#'   \item{\code{stat}}{The GRS F-statistic.}
#'   \item{\code{H0}}{A description of the null hypothesis: \code{"alpha = 0"}.}
#' }
#'
#' @references
#' Gibbons, M. R., Ross, S. A., & Shanken, J. (1989).
#' "A Test of the Efficiency of a Given Portfolio." \emph{Econometrica}, 57(5), 1121–1152.
#'
#' @examples
#' set.seed(42)
#' R1 <- matrix(rnorm(300), nrow = 100, ncol = 3)  # Benchmark returns
#' R2 <- matrix(rnorm(200), nrow = 100, ncol = 2)  # Test asset returns
#' res <- span_grs(R1, R2)
#' res$stat   # GRS test statistic
#' res$pval   # P-value
#' res$H0     # Null hypothesis tested
#'
#' @export
span_grs <- function(R1, R2) {
  # Reassign roles for clarity
  X <- R1  # Benchmark assets
  Y <- R2  # Test assets

  T <- nrow(X)
  K <- ncol(X)
  N <- ncol(Y)

  if ((T - N - K) < 1) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))
  }

  ones <- matrix(1, T, 1)
  XX <- cbind(ones, X)
  Bhat1 <- solve(crossprod(XX)) %*% crossprod(XX, Y)
  SSRu <- crossprod(Y - XX %*% Bhat1)
  SigmaU <- SSRu / T

  ahat <- Bhat1[1, ]
  mufactor <- matrix(colMeans(X), K, 1)
  ssqm <- crossprod(scale(X, center = TRUE, scale = FALSE)) / T
  mufactor_ssqm <- solve(ssqm, mufactor)

  denom <- 1 + crossprod(mufactor, mufactor_ssqm)
  num <- ((T - N - K) / N) * crossprod(ahat, solve(SigmaU, ahat))

  stat <- as.numeric(num / denom)
  pval <- 1 - pf(stat, N, T - N - K)

  list(pval = pval, stat = stat, H0 = "alpha = 0")
}
