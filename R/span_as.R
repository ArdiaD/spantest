#' Subseries-Based Cauchy Combination Test (SCT) for Spanning Hypotheses
#'
#' Computes robust p-values for testing spanning-related linear restrictions using the residual-based
#' subseries method. This function supports high-dimensional inference on hypotheses such as
#' \eqn{H_0^\delta}, \eqn{H_0^\alpha}, and their joint null \eqn{H_0^{\alpha, \delta}} using
#' a simulation-based approximation and aggregation via the Cauchy Combination Test (CCT).
#'
#' @param u A numeric vector of outcomes (e.g., returns or residuals) of length \eqn{T}.
#' @param x A numeric matrix of regressors (T x K), where the first column is used as baseline and others are compared.
#' @param ks A numeric vector of subseries exponents (e.g., \code{1/3}). Each value
#'   sets the NUMBER of subseries to \code{floor(T^k)}, so the sample is split into
#'   that many contiguous blocks of length about \code{T / floor(T^k)}. At
#'   \code{T = 250}, \code{k = 1/3} gives 6 blocks of about 42 observations and
#'   \code{k = 2/3} gives 39 blocks of about 6 -- a larger exponent means MORE and
#'   SHORTER blocks, not longer ones. This matches \eqn{l_T = \lfloor T^\psi \rfloor}
#'   in Ardia and Sessinou (2025), where \eqn{l_T} indexes the blocks and
#'   \eqn{l_T b_T = T}. Default is \code{c(1/3)}.
#' @param L A numeric vector controlling the strength of randomization applied to residual scores. Default is \code{c(0, 2)}.
#' @param wN,wcol The weights of this asset are column \code{wcol} of a draw for
#'   \code{wN} assets, so that a loop over the columns reproduces \code{f_getpv_batch()};
#'   a single asset by default.
#'
#' @return A named numeric vector of CCT p-values. Each name encodes the test type, L value, and subseries-exponent index:
#' \describe{
#'   \item{CCTd}{Test of \eqn{H_0^\delta}: no directional (slope) deviation.}
#'   \item{CCTad}{Joint test of \eqn{H_0^{\alpha,\delta}}: no intercept or slope deviation.}
#'   \item{CCTa}{Test of \eqn{H_0^\alpha}: no intercept deviation.}
#' }
#' Names are suffixed with the values of \code{L} and the subseries exponent index, e.g., \code{CCTa_L2_k1}.
#'
#' @details
#' This function builds score vectors from OLS residuals and applies randomized weightings
#' (via \code{f_mult}, one draw per statistic) to simulate perturbations. These perturbed scores are passed through a subseries
#' t-test pipeline (via \code{f_testbm}), and resulting p-values are aggregated using the
#' Cauchy Combination Test (Liu & Xie, 2020), which is valid under arbitrary dependence.
#'
#' The methodology is motivated by Carlstein's (1986) subseries inference and extends it using
#' residual-based oracle approximations that remain valid when K increases with T.
#'
#' @references
#' \insertRef{ArdiaSessinou2025}{spantest} \cr
#'
#' @importFrom stats cov na.omit pf pnorm pt qnorm rnorm runif
#'
#' @keywords internal
#'
#' @noRd
#'
f_getpv <- function(u, x, ks = c(1/3), L = c(0, 2), seed = 123L, wN = 1L, wcol = 1L) {
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

  res_ew <- res * ew
  score <- score2 <- matrix(res_ew, ncol = 1)

  # Score processing: randomized perturbation then subseries Cauchy p-value. Each
  # statistic has its own weights: hyp names the weights ("D" delta, "A" alpha) of
  # each column of score_mat; one draw from seed.
  f_process_scores <- function(score_mat, hyp, sn) {
    w <- 1
    if (sn > 0L) {
      W <- f_mult(nrow(score_mat), wN, sn, cseed = seed)
      w <- vapply(hyp, function(h) W[[h]][, wcol], numeric(nrow(score_mat)))
    }
    scoreb <- score_mat * w
    out <- vapply(ks, function(k) f_testbm(scoreb, k = k)[1], numeric(1))
    return(out)
  }

  # Delta = 0
  test1 <- lapply(L, function(sn_value) {
    f_process_scores(score, hyp = "D", sn = sn_value)
  })

  # Delta = alpha = 0 (joint)
  score_combo <- cbind(score2, res * ew1)
  test2 <- lapply(L, function(sn_value) {
    f_process_scores(score_combo, hyp = c("D", "A"), sn = sn_value)
  })

  # Alpha = 0
  score_alpha <- matrix(res * ew1, ncol = 1)
  test4 <- lapply(L, function(sn_value) {
    f_process_scores(score_alpha, hyp = "A", sn = sn_value)
  })

  result <- unlist(c(test1, test2, test4))

  # Dynamic naming: CCT{d,ad,a} x L x k
  test_names <- paste0(
    rep(c("CCTd", "CCTad", "CCTa"),
        each = length(L) * length(ks)),
    "_L", rep(L,
              each = length(ks),
              times = length(c("CCTd", "CCTad", "CCTa"))),
    "_k", rep(seq_along(ks),
              times = length(L) * length(c("CCTd", "CCTad", "CCTa")))
  )

  names(result) <- test_names
  return(result)
}

#' Vectorized per-asset CCT p-values (batch over the test cross-section)
#'
#' Computes the same per-asset subseries CCT p-values as \code{f_getpv()} but for
#' all \eqn{N} test assets at once. Because every test asset shares the same
#' benchmark design, the benchmark-only QR decompositions are formed once and the
#' three swap regressions are obtained by Frisch--Waugh partialling, turning the
#' original \eqn{O(N)} loop of full QR factorizations into a handful of shared
#' factorizations plus vectorized arithmetic. It is numerically equivalent to
#' looping \code{f_getpv(test[, j], bench, wN = N, wcol = j)} over the columns of
#' \code{test}.
#'
#' @param bench Numeric \eqn{T \times K} matrix of benchmark returns.
#' @param test  Numeric \eqn{T \times N} matrix of test-asset returns.
#' @param ks,L  As in \code{f_getpv()}.
#' @param seed  Seed of the draw of the weights when \code{L > 0}.
#'
#' @return A named list; each element is a length-\eqn{N} vector of per-asset
#' p-values, named as \code{CCT{d,ad,a}_L{L}_k{i}}.
#'
#' @keywords internal
#'
#' @noRd
#'
f_getpv_batch <- function(bench, test, ks = c(1/3), L = c(0, 2), seed = 123L) {
  x  <- bench
  Tn <- nrow(x)
  K  <- ncol(x)
  N  <- ncol(test)
  x1 <- x[, 1]
  one <- rep(1, Tn)
  Y  <- test - x1                                   # T x N: y_j = test_j - x1

  if (K == 1) {
    xc   <- matrix(nrow = Tn, ncol = 0)
    xaug <- x
  } else {
    xc   <- sweep(x[, -1, drop = FALSE], 1, x1, "-")
    xaug <- cbind(x1, xc)
  }

  # Benchmark-only QR factorizations (shared across all test assets)
  qr_main <- qr(cbind(1, xaug))   # residualize on [1, x]
  qr_bd   <- qr(cbind(1, xc))     # FWL base for ew  (add y as the extra column)
  qr_ba   <- qr(xaug)             # FWL base for ew1 (add y as the extra column)

  # res_j = residual of y_j on [1, x]  (batched over all assets)
  res <- qr.resid(qr_main, Y)

  # ew_j  = residual of x1 on [1, xc, y_j]  via Frisch-Waugh on base [1, xc]
  x1_p <- qr.resid(qr_bd, x1)
  yp1  <- qr.resid(qr_bd, Y)
  b1   <- colSums(yp1 * x1_p) / colSums(yp1 * yp1)
  ew   <- x1_p - sweep(yp1, 2, b1, "*")

  # ew1_j = residual of 1 on [x, y_j]  via Frisch-Waugh on base x
  one_p <- qr.resid(qr_ba, one)
  yp2   <- qr.resid(qr_ba, Y)
  b2    <- colSums(yp2 * one_p) / colSums(yp2 * yp2)
  ew1   <- one_p - sweep(yp2, 2, b2, "*")

  D <- res * ew    # delta scores, T x N
  A <- res * ew1   # alpha scores, T x N

  # For each (L, k): student p-values for every asset via one f_ttest per score
  # matrix (columns are independent), then Cauchy-combine the two for CCTad.
  # For L > 0 every statistic has its own weights, independent across assets and
  # between the alpha and the delta score of an asset (f_mult), one draw from seed.
  out <- list()
  for (sn in L) {
    W  <- f_mult(Tn, N, sn, cseed = seed)
    Dw <- D * W$D
    Aw <- A * W$A
    for (ki in seq_along(ks)) {
      sD <- unname(f_ttest(Dw, k = ks[ki])$student)
      sA <- unname(f_ttest(Aw, k = ks[ki])$student)
      out[[paste0("CCTd_L",  sn, "_k", ki)]] <- sD
      out[[paste0("CCTa_L",  sn, "_k", ki)]] <- sA
      out[[paste0("CCTad_L", sn, "_k", ki)]] <- vapply(seq_len(N),
                                                       function(j) f_cauchypv(c(sD[j], sA[j])),
                                                       numeric(1))
    }
  }
  return(out)
}

#' Ardia and Sessinou (2025) Subseries-Based Cauchy Combination Test (SCT) for Spanning
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
#'   \item{\code{ks}}{Numeric vector of subseries exponents; each sets the NUMBER of blocks to \code{floor(T^k)}, each of length about \code{T / floor(T^k)}; default \code{c(1/3)}.}
#'   \item{\code{L}}{Numeric vector of the number of multiplier factors in the weights of each statistic (0: no weights); default \code{c(0, 2)}. See \sQuote{Weights}.}
#'   \item{\code{seed}}{Seed of the draw of the weights when \code{L > 0}; default \code{123}. A call is one draw; a simulation should pass a different seed in each replication; see \sQuote{One draw per call} and \sQuote{Simulations}.}
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
#' For each \code{k} in \code{ks}, the sample is partitioned into \code{floor(T^k)}
#' contiguous, NON-overlapping subseries -- an exact partition of \code{1:T}, not a
#' moving block.
#' Residual perturbations controlled by \code{L} generate test statistics robust to serial and
#' cross-sectional dependence and conditional heteroskedasticity. Resulting sub-p-values are aggregated
#' by the Cauchy Combination Test (CCT), which remains valid under dependence and retains power in high dimensions.
#'
#' @section Weights:
#' When \code{L > 0} the score of each statistic is multiplied by its own weights
#' \eqn{\kappa_{i,t} = \prod_{l=1}^{L} \kappa_{l,i,t}}, \eqn{\kappa_{l,i,t} \sim N(1, 1)},
#' drawn independently of the data and of each other across factors, dates, assets
#' and between the alpha and the delta score of an asset. Since version 1.4-3 this
#' is the definition of the theory; up to version 1.4-2 one weight vector was
#' shared by every asset and by both scores. Sharing it leaves the cross-sectional
#' dependence of the statistics intact, while independent weights shrink it, which
#' matters when the residuals carry a pervasive common factor: at \code{L = 2},
#' \eqn{T = 250}, \eqn{K} = 2 and 10 and \eqn{N} = 100 and 1000, the three tests
#' reject a true null in 6.0--8.5\% of 1000 samples with the shared vector and in
#' 2.4--4.8\% with independent weights. The results at \code{L = 0} are unchanged.
#'
#' @section One draw per call:
#' At \code{L > 0} a call returns the test for one draw of the weights, from
#' \code{seed}: the alpha weights first, then the delta weights. The test is valid
#' for any draw --- its size is correct --- but the p-value it returns on one data
#' set varies from draw to draw, since a Cauchy average of independent terms does
#' not concentrate as \eqn{N} grows; report the seed with the result. R scrambles
#' the seed when initialising the Mersenne-Twister, so draws from consecutive
#' seeds behave as independent draws. Up to version 1.4-2 a control \code{B}
#' merged \code{B} draws per asset with the Cauchy rule; it was removed in 1.4-3.
#'
#' @section Simulations:
#' A Monte Carlo study should draw new weights in each replication, for instance
#' \code{seed = base + r} in replication \code{r}. With the same seed in every
#' replication, every cell of the study is computed on one set of weights: the
#' results are one realisation of the weights, common to all cells, and the Monte
#' Carlo standard errors understate their uncertainty. The weights enter through
#' sums over each subseries, so the effect on size is small but not nil (about one
#' point at \eqn{T = 250}, in either direction), and the effect on power can be
#' larger: in the simulation study of Ardia and Sessinou, the weights of seed 123
#' overstated the power at \code{L = 2} by 10 to 25 points at moderate
#' alternatives.
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
#' @export
span_as <- function(bench, test, control = list()) {

  # Set control parameters
  if ("B" %in% names(control))
    stop("control$B was removed in spantest 1.4-3: a call is one draw of the weights; ",
         "pass another seed for another draw.", call. = FALSE)
  con <- list(ks = c(1/3), L = c(0, 2), seed = 123L)
  con[names(control)] <- control
  k_values <- con$ks
  l_values <- con$L
  # seed must be a whole number: as.integer() would silently truncate a fraction
  stopifnot(length(con$seed) == 1L, is.finite(con$seed),
            isTRUE(all.equal(con$seed, round(con$seed))))

  # Generate explicit template names
  test_types <- c("CCTd", "CCTad", "CCTa")
  template_names <- paste0(
    rep(test_types, each = length(l_values) * length(k_values)),
    "_L", rep(l_values, each = length(k_values), times = length(test_types)),
    "_k", rep(seq_along(k_values), times = length(l_values) * length(test_types))
  )

  # Per-asset p-values for the whole cross-section, computed in batch.
  pv <- f_getpv_batch(bench, test, ks = k_values, L = l_values, seed = as.integer(con$seed))

  # Combine per-asset p-values across the cross-section via the Cauchy method.
  combined_results <- vapply(template_names,
                             function(nm) f_cauchypv(na.omit(pv[[nm]])),
                             numeric(1))
  n_na <- max(vapply(template_names, function(nm) sum(is.na(pv[[nm]])), numeric(1)))
  if (n_na > 0)
    warning(sprintf("span_as(): %d of %d test assets have no p-value and are left out of the combination.",
                    as.integer(n_na), ncol(test)), call. = FALSE)

  return(as.list(combined_results))
}
