#' Randomized Alpha Test of Massacci, Sarno, Trapani and Vallarino
#'
#' Implements the randomized test of \insertCite{MassacciEtAl2026;textual}{spantest}
#' for the joint null that every intercept is zero in a linear factor pricing
#' model. The test is built from equation-by-equation estimation, needs no
#' covariance matrix, and allows \eqn{N} to grow faster than \eqn{T}. Both the
#' one-shot test of their Theorem 3.1 and the derandomized decision rule of their
#' Section 3.2 are returned.
#'
#' @param R1 Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param R2 Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#' @param control Optional list:
#' \describe{
#'   \item{\code{nu}}{Number of finite moments the data are assumed to admit,
#'     \eqn{\nu \ge 4}; default \code{4}. See \sQuote{Choosing nu}.}
#'   \item{\code{tau}}{Nominal level at which the derandomized rule decides;
#'     default \code{0.05}. It does not affect \code{pval}.}
#'   \item{\code{B}}{Number of replications the derandomized threshold is built
#'     for; default \code{floor(log(N)^2)}, their guideline. The replications
#'     themselves are not simulated (see \sQuote{Details}).}
#'   \item{\code{seed}}{Seed of the single perturbation draw behind \code{stat}
#'     and \code{pval}; default \code{123}. The caller's RNG state is restored.
#'     A simulation must pass a different seed in each replication; see
#'     \sQuote{Simulations}.}
#' }
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{pval}}{P-value of the one-shot test, from the Gumbel limit.}
#'   \item{\code{stat}}{The one-shot statistic \eqn{Z_{N,T}}.}
#'   \item{\code{H0}}{Null hypothesis description, \code{"alpha = 0"}.}
#'   \item{\code{crit}}{Critical value \eqn{c_\tau} of the one-shot test.}
#'   \item{\code{Q}}{Derandomized quantity \eqn{Q_{N,T,\infty}(\tau)}, in closed form.}
#'   \item{\code{logQ}}{Its logarithm, which stays finite when \code{Q} underflows.}
#'   \item{\code{reject}}{Decision of the derandomized rule: \code{TRUE} when
#'     \code{Q < (1 - tau) - B^(-1/4)}.}
#'   \item{\code{nu}, \code{B}}{The settings used.}
#' }
#' All components except \code{H0} are \code{NA} when the test does not apply.
#'
#' @details
#' Write \eqn{\hat\alpha_i} for the OLS intercept of test asset \eqn{i} on the
#' benchmarks and \eqn{\hat s_{NT}} for the pooled root mean squared residual
#' over the whole panel. The test rescales the intercepts into
#' \deqn{\psi_i = \left(T^{1/\nu}\,|\hat\alpha_i| / \hat s_{NT}\right)^{\nu/2}
#'             = \sqrt{T}\,\left(|\hat\alpha_i| / \hat s_{NT}\right)^{\nu/2},}
#' which drifts to zero under the null and diverges under the alternative,
#' perturbs them with independent standard normals, and takes the largest,
#' \eqn{Z_{N,T} = \max_i (\psi_i + \omega_i)}. Under the null \eqn{Z_{N,T}} is
#' the maximum of \eqn{N} standard normals, so
#' \eqn{a_N^{-1}(Z_{N,T} - b_N)} is asymptotically Gumbel with the usual
#' extreme-value constants, giving both \code{crit} and \code{pval}. Only rates
#' of convergence are used, which is what frees the procedure from estimating a
#' covariance matrix and lets \eqn{N} outgrow \eqn{T}.
#'
#' Because the perturbation does not vanish asymptotically, two researchers
#' applying the one-shot test to the same data can reach different conclusions.
#' Their remedy is to repeat the perturbation \eqn{B} times and record the share
#' \eqn{Q_{N,T,B}(\tau)} of replications that do not reject, then decide against
#' the null when that share falls below \eqn{(1-\tau) - f(B)} with
#' \eqn{f(B) = B^{-1/4}}. We evaluate the share exactly rather than by
#' simulation: the perturbations are independent of each other and of the data,
#' so conditionally on the sample
#' \deqn{Q_{N,T,\infty}(\tau) = \prod_{i=1}^{N} \Phi(c_\tau - \psi_i),}
#' computed here on the log scale. This is the \eqn{B \to \infty} limit of their
#' average, so it removes the simulation error and the residual dependence on the
#' draws while leaving the rule itself untouched; \eqn{B} still enters through
#' the threshold \eqn{f(B)}. It also costs one pass over the cross-section
#' instead of \eqn{B} regressions.
#'
#' The Gumbel limit is an approximation in \eqn{N}, and a poor one when \eqn{N}
#' is small: at \eqn{N = 2} the one-shot test rejects about a quarter of the time
#' at a nominal 5%. The function therefore returns \code{NA} below
#' \eqn{N = 10}, and \eqn{N} in the hundreds is where the procedure is meant to
#' operate.
#'
#' @section Choosing nu:
#' \eqn{\nu} is the number of finite moments assumed, entering only through the
#' exponent \eqn{\nu/2}; the theory needs \eqn{\nu \ge 4} and the admissible
#' growth of \eqn{N} relative to \eqn{T} widens with it. It can be estimated with
#' a tail-index estimator, or bounded below by testing \eqn{E|y_{i,t}|^{\nu_0}}
#' for a candidate \eqn{\nu_0}. The default here is the smallest value the theory
#' allows, \eqn{\nu = 4}, which is what the authors recommend when \eqn{T} is too
#' short for reliable tail inference; their simulations use \eqn{\nu = 5}.
#'
#' @section Simulations:
#' The size of the one-shot test is a probability over the data and the
#' perturbation together, so a Monte Carlo study must draw a new perturbation in
#' each replication, for instance \code{seed = base + r} in replication \code{r}.
#' With the same seed in every replication, all replications share one vector
#' \eqn{\omega} and the rejection rate measures the test conditionally on that
#' draw. Under the null at small \eqn{K} the \eqn{\psi_i} are close to zero and
#' the statistic is essentially \eqn{\max_i \omega_i}, so that rate is near 0 or
#' near 1 whatever the data: with the default seed and \eqn{N = 1000},
#' \eqn{\max_i \omega_i = 3.24} against a critical value of 3.98, and the
#' rejection rate is 0 instead of about 4%. The derandomized rule is computed in
#' closed form and does not depend on the seed.
#'
#' @references
#' \insertRef{MassacciEtAl2026}{spantest}
#'
#' @examples
#' set.seed(123)
#' R1 <- matrix(rnorm(3 * 250), 250, 3)     # benchmarks: T=250, K=3
#' R2 <- matrix(rnorm(100 * 250), 250, 100) # tests:      T=250, N=100
#' out <- span_mstv(R1, R2)
#' out$pval; out$Q; out$reject
#'
#' @family Alpha Spanning Tests
#'
#' @importFrom stats pnorm qr qr.solve
#' @export
span_mstv <- function(R1, R2, control = list()) {

  con <- list(nu = 4, tau = 0.05, B = NULL, seed = 123L)
  con[names(control)] <- control
  stopifnot(
    "nu must be a single number >= 4" =
      length(con$nu) == 1L && is.finite(con$nu) && con$nu >= 4,
    "tau must be a single number in (0, 1)" =
      length(con$tau) == 1L && is.finite(con$tau) && con$tau > 0 && con$tau < 1,
    "seed must be a single whole number" =
      length(con$seed) == 1L && is.finite(con$seed) &&
      isTRUE(all.equal(con$seed, round(con$seed)))
  )

  R1 <- as.matrix(R1)
  R2 <- as.matrix(R2)
  Tn <- nrow(R2)
  N  <- ncol(R2)
  K  <- ncol(R1)
  stopifnot("R1 and R2 must have the same number of rows" = nrow(R1) == Tn)

  na_out <- list(pval = NA_real_, stat = NA_real_, H0 = "alpha = 0",
                 crit = NA_real_, Q = NA_real_, logQ = NA_real_,
                 reject = NA, nu = con$nu, B = NA_integer_)

  # The reference distribution is the Gumbel limit of the maximum of N normals.
  # It is an asymptotic statement in N and is worthless for a handful of assets:
  # at N = 2 the one-shot test rejects 26% of the time at a nominal 5%. Returning
  # NA is more useful than a number no one should act on.
  if (N < 10L || Tn - K - 1L < 1L) return(na_out)

  X  <- cbind(1, R1)
  # One QR for the whole cross-section: the N equations share their regressors,
  # so this is a single decomposition rather than N regressions, and it avoids
  # squaring the conditioning of X as the normal equations would at large K.
  qrX <- qr(X)
  if (qrX$rank < ncol(X)) return(na_out)
  bh <- qr.solve(qrX, R2)
  E  <- R2 - X %*% bh

  # s_NT pools the residuals across the whole panel rather than per asset: any
  # scale that removes the unit of measurement serves, and the cross-sectional
  # average smooths away spikes in individual variances.
  s_NT <- sqrt(mean(E^2))
  if (!is.finite(s_NT) || s_NT <= 0) return(na_out)

  # (T^{1/nu} |alpha| / s)^{nu/2} = sqrt(T) (|alpha| / s)^{nu/2}: the same
  # quantity with one power per element instead of two, and no intermediate that
  # overflows when nu is large.
  psi <- sqrt(Tn) * (abs(bh[1L, ]) / s_NT)^(con$nu / 2)

  bN <- sqrt(2 * log(N)) - (log(log(N)) + log(4 * pi)) / (2 * sqrt(2 * log(N)))
  aN <- bN / (1 + bN^2)
  cr <- bN - aN * log(-log(1 - con$tau))

  # The one-shot statistic is one realisation of a randomized test, so the draw
  # is seeded; the caller's stream is restored so that anything run afterwards --
  # in particular the sign-flip tests of the simulation study -- is unaffected.
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = globalenv())
  set.seed(as.integer(con$seed))
  Z <- max(psi + stats::rnorm(N))
  if (has_seed) {
    assign(".Random.seed", old_seed, envir = globalenv())
  } else {
    rm(".Random.seed", envir = globalenv())
  }

  # Derandomization in closed form. The perturbations are independent across
  # assets and of the sample, so conditionally on the data
  #   P(max_i (psi_i + w_i) <= c) = prod_i Phi(c - psi_i),
  # which is the B -> infinity limit of the share of replications that do not
  # reject. Accumulated on the log scale: the product runs over N terms and
  # underflows to exactly zero well before its logarithm ceases to be informative.
  B    <- if (is.null(con$B)) max(1L, as.integer(floor(log(N)^2))) else as.integer(con$B)
  logQ <- sum(stats::pnorm(cr - psi, log.p = TRUE))
  Q    <- exp(logQ)

  list(pval   = as.numeric(1 - exp(-exp(-(Z - bN) / aN))),
       stat   = as.numeric(Z),
       H0     = "alpha = 0",
       crit   = as.numeric(cr),
       Q      = as.numeric(Q),
       logQ   = as.numeric(logQ),
       reject = isTRUE(Q < (1 - con$tau) - B^(-1 / 4)),
       nu     = con$nu,
       B      = B)
}
