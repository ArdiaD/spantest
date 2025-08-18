#' Subseries-Based Cauchy Combination Test (SCT) for Spanning Hypotheses
#'
#' Computes robust p-values for testing spanning-related linear restrictions using the residual-based
#' subseries method. This function supports high-dimensional inference on hypotheses such as
#' \eqn{H_0^\delta}, \eqn{H_0^\alpha}, and their joint null \eqn{H_0^{\alpha, \delta}} using
#' a simulation-based approximation and aggregation via the Cauchy Combination Test (CCT).
#'
#' @param u A numeric vector of outcomes (e.g., returns or residuals) of length \eqn{T}.
#' @param x A numeric matrix of regressors (T x K), where the first column is used as baseline and others are compared.
#' @param ks A numeric vector of subseries exponents (e.g., \code{1/3}); each value defines subseries size as \code{floor(T^k)}. Default is \code{c(1/3)}.
#' @param L A numeric vector controlling the strength of randomization applied to residual scores. Default is \code{c(0, 2)}.
#'
#' @return A named numeric vector of CCT p-values. Each name encodes the test type, L value, and subseries size:
#' \describe{
#'   \item{CCTd}{Test of \eqn{H_0^\delta}: no directional (slope) deviation.}
#'   \item{CCTad}{Joint test of \eqn{H_0^{\alpha,\delta}}: no intercept or slope deviation.}
#'   \item{CCTa}{Test of \eqn{H_0^\alpha}: no intercept deviation.}
#' }
#' Names are suffixed with the values of \code{L} and the subseries exponent index, e.g., \code{CCTa_L2_k1}.
#'
#' @details
#' This function builds score vectors from OLS residuals and applies randomized weightings
#' (via \code{f_prods}) to simulate perturbations. These perturbed scores are passed through a subseries
#' t-test pipeline (via \code{f_testbm}), and resulting p-values are aggregated using the
#' Cauchy Combination Test (Liu & Xie, 2020), which is valid under arbitrary dependence.
#'
#' The methodology is motivated by Carlstein's (1986) subseries inference and extends it using
#' residual-based oracle approximations that remain valid when K increases with T.
#'
#' @references
#' \insertRef{ArdiaSessinou2025}{spantest} \cr
#'
#' @importFrom stats cov na.omit pf pnorm pt qnorm rnorm runif sd t.test
#' @importFrom utils head
#'
#' @keywords internal
#'
#' @noRd
#'
f_getpv <- function(u, x, ks = c(1/3), L = c(0, 2)) {

  stop("span_as NOT YET AVAILABLE")
  result <- NA
  return(result)
}

#' Ardia and Sessinou (2024) Subseries-Based Cauchy Combination Test (SCT) for Spanning
#'
#' Computes robust p-values for linear spanning restrictions using a residual-based
#' subseries procedure with Cauchy Combination Test (CCT) aggregation. Supports high-dimensional
#' inference for \eqn{H_0^\delta} (variance/spread slopes), \eqn{H_0^\alpha} (intercepts),
#' and the joint null \eqn{H_0^{\alpha,\delta}}.
#'
#' @param bench Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param test  Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#' @param control Optional list passed to internal computation:
#' \describe{
#'   \item{\code{ks}}{Numeric vector of subseries exponents (block size approximately \code{floor(T^k)}); default \code{c(1/3)}.}
#'   \item{\code{L}}{Numeric vector of perturbation scales for randomized projections; default \code{c(0, 2)}.}
#' }
#'
#' @return A named list of global (combined) p-values. Names encode hypothesis and settings:
#' \itemize{
#'   \item \code{CCTd_L{L}_k{i}} — variance (slope) spanning, \eqn{\delta = 0};
#'   \item \code{CCTa_L{L}_k{i}} — alpha spanning, \eqn{\alpha = 0};
#'   \item \code{CCTad_L{L}_k{i}} — joint mean–variance spanning, \eqn{\alpha = 0,\ \delta = 0}.
#' }
#'
#' @details
#' For each \code{k} in \code{ks}, data are partitioned into overlapping subseries of length \code{floor(T^k)}.
#' Residual perturbations controlled by \code{L} generate test statistics robust to serial and
#' cross-sectional dependence and conditional heteroskedasticity. Resulting sub-p-values are aggregated
#' by the Cauchy Combination Test (CCT), which remains valid under dependence and retains power in high dimensions.
#'
#' @references
#' \insertRef{ArdiaSessinou2025}{spantest} \cr
#'
#' @examples
#' set.seed(123)
#' bench <- matrix(rnorm(300), 100, 3)
#' test  <- matrix(rnorm(200), 100, 2)
#' out <- span_as(bench, test)
#' out$CCTa_L0_k1; out$CCTd_L0_k1; out$CCTad_L0_k1
#'
#' @family Alpha Spanning Tests
#' @family Variance Spanning Tests
#' @family Joint Mean-Variance Spanning Tests
#'
#' @importFrom stats na.omit
#' @importFrom utils head
#' @keywords internal
#'
#' @noRd
#'
span_as <- function(bench, test, control = list()) {

  stop("span_as NOT YET AVAILABLE")
  combined_results <- NA

  return(as.list(combined_results))
}
