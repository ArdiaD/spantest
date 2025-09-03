#' Gungor–Luger Joint Mean–Variance Spanning Test (2016)
#'
#' Tests the joint null \eqn{H_0:\ \alpha = 0,\ \delta = 0} that benchmark assets
#' span both intercepts and slopes of the test assets, allowing for
#' heteroskedasticity, serial dependence, and time-varying covariances.
#' Following Gungor & Luger (2016), the procedure uses a Monte Carlo (MC)
#' test based on an \eqn{F_{\max}} statistic with residual sign-flip simulations,
#' yielding Least-Favorable (LMC) and Balanced (BMC) MC p-values and a
#' three-way decision rule.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times N}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times K}.
#' @param control List of options:
#' \describe{
#'   \item{\code{totsim}}{Number of MC simulations (default \code{500}).}
#'   \item{\code{pval_thresh}}{Significance level for decisions (default \code{0.05}).}
#'   \item{\code{do_trace}}{Logical; print progress (default \code{TRUE}).}
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{pval_LMC}}{Least-Favorable MC p-value.}
#'   \item{\code{pval_BMC}}{Balanced MC p-value.}
#'   \item{\code{stat}}{Observed \eqn{F_{\max}} statistic.}
#'   \item{\code{Decisions}}{Decision code: \code{1} = Accept, \code{0} = Reject, \code{NA} = Inconclusive.}
#'   \item{\code{Decisions_string}}{Text label: \code{"Accept"}, \code{"Reject"}, or \code{"Inconclusive"}.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"alpha = 0 and delta = 0"}.}
#' }
#'
#' @details
#' LMC/BMC follow Gungor & Luger’s MC framework with residual sign-flip draws under
#' the null. The rule is: Accept if \code{pval_LMC > alpha}; Reject if
#' \code{pval_BMC <= alpha}; otherwise Inconclusive. This approach is robust in
#' high-dimensional and time-varying volatility settings where classical joint
#' spanning tests can be unreliable.
#'
#' @references
#' \insertRef{GungorLuger2016}{spantest} \cr
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)
#' R2 <- matrix(rnorm(200), 100, 2)
#' out <- span_gl_ad(R1, R2, control = list(totsim = 100, do_trace = FALSE))
#' out$Decisions_string; out$pval_LMC; out$pval_BMC
#'
#' @family Joint Mean-Variance Spanning Tests
#'
#' @importFrom stats rnorm runif
#' @export
span_gl_ad <- function(R1, R2, control = list()) {

  con <- list(totsim = 500, pval_thresh = 0.05, do_trace = TRUE)
  con[names(control)] <- control
  thresh <- con$pval_thresh

  X <- R1
  Y <- R2

  totsim <- con$totsim

  TT <- nrow(X)
  K <- ncol(X)
  N <- ncol(Y)
  ones <- matrix(1,TT,1)

  XX <- matrix(cbind(ones, X), TT, (K + 1))

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
  HXtHt <- HXt %*% t(H)
  HXtHt_inv <- solve(HXtHt)
  XtH <- crossprod(Xtemp, t(H))

  Bhat0 <- Bhat1 - XtH %*% HXtHt_inv %*% (H %*% Bhat1 - C) # Restricted
  Ehat0 <- Y - XX %*% Bhat0
  SSRr <- crossprod(Ehat0)
  SigmaR <- SSRr / TT

  # MC tests optimization
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
              H0 = "alpha = 0 and delta = 0")

  return(out)
}
