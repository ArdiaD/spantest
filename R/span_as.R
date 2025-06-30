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
#' Carlstein, E. (1986). "The Use of Subseries Values for Estimating the Variance of a General Statistic from a Stationary Sequence." \emph{Annals of Statistics}, 14(3), 1171–1179.
#' Liu, Y., & Xie, J. (2020). "Cauchy Combination Test: A Powerful Test With Analytic p-Value Calculation Under Arbitrary Dependency Structures." \emph{JASA}, 115(529), 393–402.
#' Chernozhukov, V., et al. (2021). "High-Dimensional CLTs, Causal Inference, and Robust Hypothesis Testing." \emph{Annals of Statistics}.
#' Gungor, S., & Luger, R. (2016). "A New Test of Mean–Variance Spanning." \emph{Econometric Reviews}, 35(6), 1146–1171.
#'
#' @keywords internal
#'
#' @noRd
#'
f_getpv <- function(u, x, ks = c(1/3), L = c(0, 2)) {
  one <- rep(1, nrow(x))
  y <- u - x[, 1]

  if (ncol(x) == 1) {
    x_centered <- matrix(nrow = nrow(x), ncol = 0)  # empty
    x_augmented <- x  # only one column
  } else {
    x_centered <- sweep(x[, -1, drop = FALSE], 1, x[, 1], "-")
    x_augmented <- cbind(x[, 1], x_centered)
  }

  x <- x_augmented

  X_main <- cbind(1, x)
  qr_main <- qr(X_main)
  res <- qr.resid(qr_main, y)

  xy_centered <- cbind(y, x_centered)

  f_compute_resid <- function(X, y) {
    qr <- qr(X)
    qr.resid(qr, y)
  }

  ew <- f_compute_resid(cbind(1, xy_centered), x[, 1])
  ew1 <- f_compute_resid(cbind(y, x), one)
  ew3 <- f_compute_resid(xy_centered, x[, 1])

  res_ew <- res * ew
  score <- score2 <- matrix(res_ew, ncol = 1)

  X_reg3 <- x
  qr_reg3 <- qr(X_reg3)
  res3 <- qr.resid(qr_reg3, y)
  score3 <- matrix(res3 * ew3, ncol = 1)

  # Score processing function
  process_scores <- function(score_mat, sn) {
    scoreb <- score_mat * f_prods(score_mat, sn)
    vapply(ks, function(k) f_testbm(scoreb, k = k)[1], numeric(1))
  }

  # Delta=0
  test1 <- lapply(L, function(sn_value) {
    process_scores(score, sn = sn_value)
  })

  # Delta=alpha=0
  score_combo <- cbind(score2, res * ew1)
  test2 <- lapply(L, function(sn_value) {
    process_scores(score_combo, sn = sn_value)
  })

  # Delta=0/alpha=0
  test3 <- lapply(L, function(sn_value) {
    process_scores(score3, sn = sn_value)
  })

  # Alpha=0
  score_alpha <- matrix(res * ew1, ncol = 1)
  test4 <- lapply(L, function(sn_value) {
    process_scores(score_alpha, sn = sn_value)
  })

  #test5 <- apply(cbind(test3, test4), 1, f_cauchypv)

  #result <- c(test1, test2, test3, test4, test5)
  result <- unlist(c(test1, test2, test4))

  # Dynamic naming
  f_ksuffix <- function(ks) {
    paste0("_k", seq_along(ks))
  }

  suffixes <- f_ksuffix(ks)
  #test_names <- paste0(rep(c("CCTd", "CCTad", "CCTd.a", "CCTa", "CCTd.a.a"), each = length(ks)), suffixes)
  test_names <- paste0(
    # Test types (CCTd, CCTad, CCTa)
    rep(c("CCTd", "CCTad", "CCTa"),
        each = length(L) * length(ks)),

    # L values
    "_L", rep(L,
              each = length(ks),
              times = length(c("CCTd", "CCTad", "CCTa"))),

    # k indices
    "_k", rep(seq_along(ks),
              times = length(L) * length(c("CCTd", "CCTad", "CCTa")))
  )

  names(result) <- test_names
  return(result)
}

#' Subseries-Based Cauchy Combination Test (SCT) for Spanning Hypotheses
#'
#' Computes robust p-values for testing spanning-related linear restrictions using the residual-based
#' subseries method. This function supports high-dimensional inference on hypotheses such as
#' \eqn{H_0^\delta}, \eqn{H_0^\alpha}, and their joint null \eqn{H_0^{\alpha, \delta}} using
#' a simulation-based approximation and aggregation via the Cauchy Combination Test (CCT).
#'
#' @param bench A numeric matrix of benchmark asset returns (T x K), used as the spanning base.
#' @param test A numeric matrix of test asset returns (T x N), to be evaluated for spanning.
#' @param control Optional named list of control parameters passed to f_getpv (internal function):
#' \describe{
#'   \item{\code{ks}}{A numeric vector of subseries exponents controlling block size scaling, \code{c(1/3)} (default).}
#'   \item{\code{L}}{A numeric vector of perturbation scaling constants controlling randomized projection behavior \code{c(0, 2)} (default).}
#' }
#'
#' @return A named list of global p-values for each hypothesis type and configuration:
#' \describe{
#'   \item{\code{CCTd}}{Global p-value for the null \eqn{\delta = 0} (no slope deviation).}
#'   \item{\code{CCTa}}{Global p-value for the null \eqn{\alpha = 0} (no intercept deviation).}
#'   \item{\code{CCTad}}{Global p-value for the joint null \eqn{\alpha = 0, \delta = 0}.}
#' }
#' Each entry is suffixed by the subseries settings used, e.g., \code{CCTa_L0_k1}.
#'
#' @details
#' The function performs residual-based inference by dividing the data into overlapping subseries
#' of length approximately \code{floor(T^k)} for each exponent \code{k} in \code{ks}. Randomized
#' perturbations of residuals controlled by \code{L} generate test statistics robust to serial
#' and cross-sectional dependence. P-values from these subseries tests are combined using the
#' Cauchy Combination Test (CCT), which controls for dependence between tests and retains good
#' power in high dimensions.
#'
#' The naming convention of the output encodes the test type:
#' \itemize{
#'   \item \code{CCTd}: minimum-variance (delta) spanning test,
#'   \item \code{CCTa}: max Sharpe (alpha) spanning test,
#'   \item \code{CCTad}: Full mean variance spanning (MVS), joint test of alpha and delta.
#' }
#' Suffixes specify the perturbation scale \code{L} and block size index \code{k}.
#'
#' This method is particularly useful in large asset universes where classical F-tests or GRS
#' tests become invalid due to dimensionality and complex dependence structures.
#' However, computation time increases with sample size and number of simulations.
#'
#' @references
#' Ardia, D., & Sessinou, M. (2025). "High-Dimensional Spanning Tests Using Residual-Based Subseries Methods." \emph{Working paper}.\cr
#' Liu, Y., & Xie, J. (2020). "Cauchy Combination Test: A Powerful Test With Analytic p-Value Calculation Under Arbitrary Dependency Structures." \emph{JASA}, 115(529), 393–402.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' bench <- matrix(rnorm(300), nrow = 100, ncol = 3)  # Benchmark assets
#' test  <- matrix(rnorm(200), nrow = 100, ncol = 2)  # Test assets
#' span_as(bench, test)
span_as <- function(bench, test, control = list()) {
  # Set control parameters
  con <- list(ks = c(1/3), L = c(0, 2))
  con[names(control)] <- control
  k_values <- con$ks
  l_values <- con$L

  # Generate explicit template names
  test_types <- c("CCTd", "CCTad", "CCTa")
  template_names <- paste0(
    rep(test_types, each = length(l_values) * length(k_values)),
    "_L", rep(l_values, each = length(k_values), times = length(test_types)),
    "_k", rep(seq_along(k_values), times = length(l_values) * length(test_types))
  )

  # Preallocate matrix with proper dimensions
  pval_matrix <- matrix(NA,
                        nrow = length(template_names),
                        ncol = ncol(test),
                        dimnames = list(template_names, NULL))

  # Populate matrix with p-values
  for (i in seq_len(ncol(test))) {
    pvals <- f_getpv(test[, i], bench, ks = k_values, L = l_values)

    # Validate dimensions
    if (!identical(names(pvals), template_names)) {
      stop(paste("Name mismatch in column", i,
                 "\nExpected:", paste(head(template_names), collapse = " "),
                 "\nGot:", paste(head(names(pvals)), collapse = " ")))
    }

    pval_matrix[, i] <- unlist(pvals)
  }

  # Combine p-values
  f <- function(x) f_cauchypv(na.omit(x))
  combined_results <- apply(pval_matrix, 1, f)

  return(as.list(combined_results))
}
