#' F1 Spanning Test for Intercepts (Alpha = 0)
#'
#' Implements the F1 test proposed by Kan and Zhou (2012), which tests the null hypothesis
#' \eqn{H_0: \alpha = 0}, meaning the intercepts in the return-generating process of the test assets
#' are jointly zero. This tests whether the mean of the test asset returns lies within the
#' mean-variance frontier of the benchmark assets.
#'
#' @param R1 A numeric matrix of benchmark asset returns (T x K).
#' @param R2 A numeric matrix of test asset returns (T x N).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{pval}}{The associated p-value.}
#'   \item{\code{stat}}{The F1 test statistic.}
#'   \item{\code{H0}}{The null hypothesis tested: \code{"alpha = 0"}.}
#' }
#'
#' @references
#' Kan, R., & Zhou, G. (2012). "Tests of Mean-Variance Spanning." \emph{Annals of Economics and Finance}, 13(1), 145–193.
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)  # Benchmark returns
#' R2 <- matrix(rnorm(200), 100, 2)  # Test asset returns
#' result <- span_f1(R1, R2)
#' result$stat
#' result$pval
#' result$H0
#'
#' @export
span_f1 <- function(R1, R2) {
  R <- cbind(R1, R2)
  T <- nrow(R)
  K <- ncol(R1)
  N <- ncol(R2)

  df1 <- N
  df2 <- T - K - N
  if (df2 < 1) return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))

  mu_full <- matrix(colMeans(R), ncol = 1)
  mu_bench <- matrix(colMeans(R1), ncol = 1)

  Sigma_full <- cov(R)
  Sigma_bench <- cov(R1)

  iSigma_full <- tryCatch(solve(Sigma_full), error = function(e) return(NULL))
  iSigma_bench <- tryCatch(solve(Sigma_bench), error = function(e) return(NULL))
  if (is.null(iSigma_full) || is.null(iSigma_bench)) return(list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0"))

  a <- t(mu_full) %*% iSigma_full %*% mu_full
  a1 <- t(mu_bench) %*% iSigma_bench %*% mu_bench

  stat <- (df2 / df1) * ((a - a1) / (1 + a1))
  pval <- pf(stat, df1, df2, lower.tail = FALSE)

  list(pval = as.numeric(pval), stat = as.numeric(stat), H0 = "alpha = 0")
}
