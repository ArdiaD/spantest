#' Gungor-Luger Joint Alpha-Delta Spanning Test
#'
#' Implements the joint mean-variance spanning test \eqn{H_0: \alpha = 0, \delta = 0} proposed by
#' Gungor and Luger (2016), using a residual-based subseries Monte Carlo method.
#' This test evaluates whether both intercepts and slope coefficients of test assets are
#' spanned by the benchmark assets, allowing for heteroskedasticity, serial dependence, and
#' time-varying covariance structures.
#'
#' The test computes Least Favorable (LMC) and Balanced (BMC) p-values and provides decision rules
#' (Accept, Reject, Inconclusive) based on these values.
#'
#' @param R1 Numeric matrix (T x N) of benchmark returns.
#' @param R2 Numeric matrix (T x K) of test asset returns.
#' @param control A list of control options:
#' \describe{
#'   \item{\code{seed}}{Random seed for reproducibility (default: 1).}
#'   \item{\code{totsim}}{Number of Monte Carlo simulations (default: 500).}
#'   \item{\code{pval_thresh}}{Significance level threshold for decisions (default: 0.05).}
#'   \item{\code{do_trace}}{Logical to print progress and results (default: TRUE).}
#' }
#'
#' @return A list with:
#' \describe{
#'   \item{\code{pval_LMC}}{Least Favorable Monte Carlo p-value.}
#'   \item{\code{pval_BMC}}{Balanced Monte Carlo p-value.}
#'   \item{\code{stat}}{Observed F-max test statistic.}
#'   \item{\code{Decisions}}{Decision code: 1 = Accept, 0 = Reject, NA = Inconclusive.}
#'   \item{\code{Decisions_string}}{Textual interpretation of the decision.}
#'   \item{\code{H0}}{Null hypothesis tested: \code{"alpha = 0 and delta = 0"}.}
#' }
#'
#' @details
#' This test is recommended for high-dimensional and time-varying volatility environments where classical
#' joint spanning tests are unreliable. The residual-based subseries bootstrap handles
#' conditional heteroskedasticity and autocorrelation without restrictive distributional assumptions.
#'
#' @references
#' Gungor, S., & Luger, R. (2016). "Multivariate Tests of Mean-Variance Efficiency and Spanning With a Large Number of Assets and Time-Varying Covariances."
#' \emph{Journal of Business & Economic Statistics}, 34(2), 161-175.
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)
#' R2 <- matrix(rnorm(200), 100, 2)
#' result <- span_gl_ad(R1, R2, control = list(totsim = 100))
#' result$Decisions_string
#' result$pval_LMC
#'
#' @export
span_gl_ad <- function(R1, R2, control = list()) {

  con <- list(seed = 1, totsim = 500, pval_thresh = 0.05, do_trace = TRUE)
  con[names(control)] <- control
  thresh <- con$pval_thresh

  # Renaming to match former implementation
  X <- R1
  Y <- R2

  totsim <- con$totsim
  seed <- con$seed

  TT <- nrow(X)
  K <- ncol(X)
  N <- ncol(Y)
  ones <- matrix(1,TT,1)

  XX <- matrix(cbind(ones, X), TT, (K + 1))

  # OLS parameter estimates
  XX_crossprod <- crossprod(XX)
  Xtemp <- solve(XX_crossprod)
  Bhat1 <- Xtemp %*% crossprod(XX, Y)

  Ehat1 <- Y - XX %*% Bhat1         # Unrestricted
  SSRu <- crossprod(Ehat1)
  SigmaU <- SSRu / TT


  H <- matrix(0, 2, K + 1)
  H[1, 1] <- 1
  H[2, 2:(K + 1)] <- 1
  C <- matrix(0, 2, N)
  C[2, ] <- 1

  # Precompute constraint components
  HXt <- H %*% Xtemp
  HXtHt <- HXt %*% t(H)  # Critical fix here
  HXtHt_inv <- solve(HXtHt)
  XtH <- crossprod(Xtemp, t(H))  # Xtemp %*% t(H)

  Bhat0 <- Bhat1 - XtH %*% HXtHt_inv %*% (H %*% Bhat1 - C) # Restricted
  Ehat0 <- Y - XX %*% Bhat0
  SSRr <- crossprod(Ehat0)
  SigmaR <- SSRr / TT

  ## MC tests optimization
  diag_SSRr <- diag(SSRr)
  diag_SSRu <- diag(SSRu)
  Fmax <- max((diag_SSRr - diag_SSRu) / diag_SSRu)
  Fmax_actual <- Fmax

  LMCstats <- numeric(totsim)
  BMCstats <- numeric(totsim)
  LMCstats[totsim] <- BMCstats[totsim] <- Fmax

  Ehat0data <- Ehat0
  Bhat0data <- Bhat0

  # Precompute invariant matrices
  Xtemp_XX <- Xtemp %*% t(XX)
  premult <- Xtemp %*% t(H) %*% HXtHt_inv
  set.seed(123456 * seed)

  # Vectorized simulation block
  sign_matrix <- matrix(sign(rnorm(TT * (totsim - 1))), TT, totsim - 1)
  esim_array <- array(Ehat0data, dim = c(TT, N, totsim - 1)) *
    aperm(array(sign_matrix, dim = c(TT, totsim - 1, N)), c(1,3,2))

  Ysim_array <- array(XX %*% Bhat0data, dim = c(TT, N, totsim - 1)) + esim_array

  # Reshape for matrix operations
  Ysim_mat <- matrix(Ysim_array, TT, N * (totsim - 1))
  Bhat1_mat <- Xtemp_XX %*% Ysim_mat
  Ehat1_mat <- Ysim_mat - XX %*% Bhat1_mat
  SSRu_vec <- colSums(Ehat1_mat^2)

  # Constrained estimates
  HB_C <- array(H %*% Bhat1_mat - c(C), dim = c(nrow(H), N, totsim - 1))
  Bhat0_array <- array(Bhat1_mat, dim = c(nrow(Bhat1_mat), N, totsim - 1)) -
    array(premult %*% matrix(HB_C, nrow(H)), dim = c(nrow(Bhat1_mat), N, totsim - 1))

  Ehat0_mat <- Ysim_mat - XX %*% matrix(Bhat0_array, nrow(Bhat0_array))
  SSRr_LMC_vec <- colSums(Ehat0_mat^2)

  # Test statistics
  temp_LMC <- (SSRr_LMC_vec - SSRu_vec) / SSRu_vec
  LMCstats[1:(totsim - 1)] <- matrix(temp_LMC, ncol = N, byrow = TRUE) |> apply(1, max)

  # BMC statistics
  SSRr_BMC_vec <- colSums(matrix(esim_array^2, TT, N * (totsim - 1)))
  temp_BMC <- (SSRr_BMC_vec - SSRu_vec) / SSRu_vec
  BMCstats[1:(totsim - 1)] <- matrix(temp_BMC, ncol = N, byrow = TRUE) |> apply(1, max)

  # Optimized p-value calculation
  uu <- runif(totsim)
  GL_pval_LMC <- (totsim - f_ranklex(LMCstats, uu) + 1) / totsim
  GL_pval_BMC <- (totsim - f_ranklex(BMCstats, uu) + 1) / totsim

  # Decision logic
  Decisions <- if (GL_pval_LMC > thresh) {
    1  # ACCEPT
  } else if (GL_pval_BMC <= thresh) {
    0  # REJECT
  } else {
    NA  # INCONCLUSIVE
  }

  Decisions_string <- if (GL_pval_LMC > thresh) {
    "Accept"
  } else if (GL_pval_BMC <= thresh) {
    "Reject"
  } else {
    "Inconclusive"
  }

  if (con$do_trace) {
    cat('============================================\n')
    cat('F-Max:', Fmax_actual,
        '\nLMC p-value:', GL_pval_LMC,
        '\nBMC p-value:', GL_pval_BMC,
        '\nDecision:', Decisions, '\n')
  }

  out <- list(pval_LMC = GL_pval_LMC, pval_BMC = GL_pval_BMC,
              stat = Fmax_actual, Decisions = Decisions,
              Decisions_string = Decisions_string,
              H0 = "alpha = 0 and delta = 0") # H0 alpha et Delta = 0

  return(out)
}
