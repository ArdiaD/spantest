#' Gungor–Luger Alpha-Only Spanning Test (2016)
#'
#' Tests the null \eqn{H_0:\ \alpha = 0} that benchmark assets span the mean
#' (intercepts) of the test assets. Following Gungor & Luger (2016), the
#' procedure uses a Monte Carlo (MC) test based on an \eqn{F_{\max}} statistic
#' with residual sign-flip simulations, yielding Least-Favorable (LMC) and
#' Balanced (BMC) MC p-values and a three-way decision rule.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times N}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times K}.
#' @param control List of options:
#' \describe{
#'   \item{\code{seed}}{Integer seed for reproducibility (default \code{1}).}
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
#'   \item{\code{H0}}{Null hypothesis description, \code{"alpha = 0"}.}
#' }
#'
#' @details
#' Accept if \code{pval_LMC > alpha}; Reject if \code{pval_BMC <= alpha};
#' otherwise Inconclusive. The subseries sign-flip MC approach is robust to
#' heteroskedasticity, serial dependence, and heavy tails, making it suitable
#' for high-dimensional settings where classical alpha tests (e.g., GRS) may
#' suffer from size distortions.
#'
#' @references
#' \insertRef{GungorLuger2016}{spantest} \cr
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(300), 100, 3)
#' R2 <- matrix(rnorm(200), 100, 2)
#' out <- span_gl_a(R1, R2, control = list(totsim = 100, do_trace = FALSE))
#' out$Decisions_string; out$pval_LMC; out$pval_BMC
#'
#' @family Alpha Spanning Tests
#'
#' @importFrom stats rnorm runif
#' @export
span_gl_a <- function(R1, R2, control = list()) {

  con <- list(seed = 1, totsim = 500, do_trace = TRUE, pval_thresh = 0.05)
  con[names(control)] <- control
  thresh <- con$pval_thresh

  X <- R1
  Y <- R2

  K <- ncol(X)
  N <- ncol(Y)
  TT <- nrow(X)
  totsim <- con$totsim
  ones <- matrix(1, TT, 1)

  XX <- cbind(ones, X)

  XX_crossprod <- crossprod(XX)
  Xtemp <- solve(XX_crossprod)
  Bhat1 <- Xtemp %*% crossprod(XX, Y)

  Ehat1 <- Y - XX %*% Bhat1
  SSRu <- crossprod(Ehat1)
  SigmaU <- SSRu / TT

  H <- matrix(0, 1, K + 1)
  H[1, 1] <- 1
  C <- matrix(0, 1, N)

  HXt <- H %*% Xtemp
  HXtHt <- HXt %*% t(H)
  HXtHt_inv <- solve(HXtHt)
  HB_minus_C <- H %*% Bhat1 - C

  Bhat0 <- Bhat1 - Xtemp %*% t(H) %*% HXtHt_inv %*% HB_minus_C
  Ehat0 <- Y - XX %*% Bhat0
  SSRr <- crossprod(Ehat0)

  diag_SSRr <- diag(SSRr)
  diag_SSRu <- diag(SSRu)
  temp <- (diag_SSRr - diag_SSRu) / diag_SSRu
  Fmax_actual <- max(temp)

  LMCstats <- numeric(totsim)
  BMCstats <- numeric(totsim)
  LMCstats[totsim] <- BMCstats[totsim] <- Fmax_actual

  premult <- Xtemp %*% t(H) %*% HXtHt_inv
  Xtemp_XX <- Xtemp %*% t(XX)

  set.seed(123456 * con$seed)

  sim_count <- totsim - 1
  sign_mat <- matrix(sign(rnorm(TT * sim_count)), TT, sim_count)
  esim_array <- array(Ehat0, dim = c(TT, N, sim_count)) *
    aperm(array(sign_mat, dim = c(TT, sim_count, N)), c(1, 3, 2))

  Ysim_array <- array(XX %*% Bhat0, dim = c(TT, N, sim_count)) + esim_array
  Ysim_mat <- matrix(Ysim_array, TT, N * sim_count)
  Bhat1_mat <- Xtemp_XX %*% Ysim_mat
  Ehat1_mat <- Ysim_mat - XX %*% Bhat1_mat
  SSRu_vec <- colSums(Ehat1_mat^2)

  HB_C <- array(H %*% Bhat1_mat - c(C), dim = c(nrow(H), N, sim_count))
  Bhat0_array <- array(Bhat1_mat, dim = c(nrow(Bhat1_mat), N, sim_count)) -
    array(premult %*% matrix(HB_C, nrow(H)), dim = c(nrow(Bhat1_mat), N, sim_count))
  Ehat0_mat <- Ysim_mat - XX %*% matrix(Bhat0_array, nrow(Bhat0_array))
  SSRr_LMC_vec <- colSums(Ehat0_mat^2)

  temp_LMC <- (SSRr_LMC_vec - SSRu_vec) / SSRu_vec
  LMCstats[1:sim_count] <- matrix(temp_LMC, ncol = N, byrow = TRUE) |> apply(1, max)

  SSRr_BMC_vec <- colSums(matrix(esim_array^2, TT, N * sim_count))
  temp_BMC <- (SSRr_BMC_vec - SSRu_vec) / SSRu_vec
  BMCstats[1:sim_count] <- matrix(temp_BMC, ncol = N, byrow = TRUE) |> apply(1, max)

  uu <- runif(totsim)
  GL_pval_LMC <- (totsim - f_ranklex(LMCstats, uu) + 1) / totsim
  GL_pval_BMC <- (totsim - f_ranklex(BMCstats, uu) + 1) / totsim

  Decisions <- if (GL_pval_LMC > thresh) {
    1
  } else if (GL_pval_BMC <= thresh) {
    0
  } else {
    NA
  }

  Decisions_string <- if (GL_pval_LMC > thresh) {
    "Accept"
  } else if (GL_pval_BMC <= thresh) {
    "Reject"
  } else {
    "Inconclusive"
  }

  if (con$do_trace) {
    cat("============================================\n")
    cat("F-max:", Fmax_actual, "\n")
    cat("LMC p-value:", GL_pval_LMC, "\n")
    cat("BMC p-value:", GL_pval_BMC, "\n")
    cat("Decision:", Decisions_string, "\n")
  }

  GL <- list(
    pval_LMC = GL_pval_LMC,
    pval_BMC = GL_pval_BMC,
    stat = Fmax_actual,
    Decisions = Decisions,
    Decisions_string = Decisions_string,
    H0 = "alpha = 0"
  )

  return(GL)
}
