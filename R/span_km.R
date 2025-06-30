#' Kempf-Memmel GMVP Spanning Test
#'
#' Tests whether the Global Minimum Variance Portfolio (GMVP) constructed from the combined
#' benchmark and test assets is equivalent to the GMVP of the benchmark assets alone.
#' This test, proposed by Kempf and Memmel (2006), evaluates if adding new assets
#' significantly improves the minimum variance frontier.
#'
#' @param R1 Numeric matrix of benchmark returns (T x N), with T time observations and N benchmark assets.
#' @param R2 Numeric matrix of test asset returns (T x K), with K test assets.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{pval}}{P-value corresponding to the F-statistic under the null.}
#'   \item{\code{stat}}{F-statistic value of the test.}
#'   \item{\code{H0}}{Description of the null hypothesis tested: \code{"GMVP of benchmark = GMVP of full asset universe"}.}
#' }
#'
#' @references
#' Kempf, A., & Memmel, C. (2006). "Estimating the Global Minimum Variance Portfolio."
#' \emph{Schmalenbach Business Review}, 58(4), 332–348.
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), nrow = 100, ncol = 3)  # Benchmark assets
#' R2 <- matrix(rnorm(200), nrow = 100, ncol = 2)  # Test assets
#' result <- span_km(R1, R2)
#' print(result$pval)  # p-value of the test
#' print(result$stat)  # F-statistic
#' if (result$pval < 0.05) {
#'   cat("Reject null: GMVP of full set differs from benchmark\n")
#' } else {
#'   cat("Fail to reject null: no evidence of improvement in GMVP\n")
#' }
#'
#' @export
span_km <- function(R1, R2) {
  # Kempf & Memmel (2006) test: H0 = delta = 0
  # Step 1: combine benchmark and test asset returns
  R <- cbind(R1, R2)
  TT <- nrow(R)
  p <- ncol(R1)
  p2 <- ncol(R2)

  # Step 2: y = constant vector of 1s
  y <- rep(1, TT)

  # Step 3: create difference matrix: R[,1] - R[,-1]
  Diff <- sweep(R[, -1, drop = FALSE], 1, R[, 1], FUN = function(x, y) y - x)

  # Step 4: regression of y on Diff (no intercept)
  X <- Diff
  XtX <- crossprod(X)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) return(NULL))
  if (is.null(XtX_inv)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  beta_hat <- XtX_inv %*% crossprod(X, y)
  residuals <- y - X %*% beta_hat
  sigma2 <- drop(crossprod(residuals) / (TT - ncol(X)))

  # Step 5: F-test on the last p2 coefficients
  offset <- p - 1
  C <- cbind(matrix(0, p2, offset), diag(p2))  # p2 × (p-1 + p2)

  middle <- tryCatch(solve(C %*% XtX_inv %*% t(C)), error = function(e) return(NULL))
  if (is.null(middle)) {
    return(list(pval = NA_real_, stat = NA_real_, H0 = "delta = 0"))
  }

  theta_part <- C %*% beta_hat
  F_stat <- as.numeric(t(theta_part) %*% middle %*% theta_part / (p2 * sigma2))
  p_val <- pf(F_stat, p2, TT - ncol(X), lower.tail = FALSE)

  return(list(pval = p_val, stat = F_stat, H0 = "delta = 0"))
}
