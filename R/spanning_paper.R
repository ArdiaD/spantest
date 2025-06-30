### !!! THE FUNCTIONS IN THIS FILE ARE PROPRIETARY AND SHOULD NOT BE USED BY PACKAGE USERS !!! ###

#' Portfolio Spanning Tests: Huberman-Kandel (HK) and Cauchy Combination Tests (CCT)
#'
#' Performs portfolio spanning tests to assess whether the test asset returns in \code{R2}
#' are mean-variance spanned by the benchmark asset returns in \code{R1}.
#' This includes the Huberman and Kandel (1987) test, F-type statistics, and p-value combination
#' using the Cauchy Combination Test (CCT).
#'
#' Null hypotheses tested:
#' \itemize{
#'   \item \code{HK}: \eqn{\alpha = 0} and \eqn{\delta = 0} — full spanning (Huberman & Kandel, 1987)
#'   \item \code{F1}: \eqn{\alpha = 0} — intercept test only
#'   \item \code{F2}: \eqn{\delta = 0} — slope (directional) spanning test
#'   \item \code{CCTa}: t-test for intercept \eqn{\alpha = 0}
#'   \item \code{CCTd}: t-test for slope \eqn{\delta = 0}
#'   \item \code{CCTad}: Cauchy combination of \code{CCTa} and \code{CCTd}
#' }
#'
#' @param R1 A numeric matrix (T x K) of benchmark asset returns.
#' @param R2 A numeric matrix (T x N) of test asset returns to be evaluated for spanning.
#'
#' @return A named list where each element corresponds to a test result and contains:
#' \describe{
#'   \item{pval}{The p-value of the test.}
#'   \item{stat}{The test statistic (for HK, F1, and F2).}
#'   \item{H0}{A description of the null hypothesis being tested.}
#' }
#' !!! This function MUST not be used by package users !!!
#'
#' @references
#' Huberman, G., & Kandel, S. (1987). "Mean-Variance Spanning." \emph{Journal of Finance}, 42(4), 873–888.
#'
#' Liu, Y., Xie, J., & Li, R. (2020). "Cauchy Combination Test: A Powerful Test With Analytic p-Value Calculation Under Arbitrary Dependency Structures." \emph{Journal of the American Statistical Association}, 115(529), 393–402.
#'
#' @keywords internal
#'
#' @noRd
#'
f_hk <- function(R1, R2) {
      # Combine returns and precompute dimensions
      R <- cbind(R1, R2)
      K <- ncol(R1)
      N <- ncol(R2)
      TT <- nrow(R1)

      # CCT: Optimized multivariate regression
      Y_cct <- R2 - R1[, 1]
      X_cct <- cbind(-R1[, 1], R1[, -1] - R1[, 1])
      X_mat <- cbind(1, X_cct)

      # Precompute QR decomposition once
      qr_X <- qr(X_mat)
      Q <- qr.Q(qr_X)
      R_inv <- solve(qr.R(qr_X))
      se_factors <- sqrt(rowSums(R_inv^2))  # Efficient SE calculation

      # Initialize storage
      p_values <- matrix(nrow = 2, ncol = ncol(Y_cct))

      # Vectorized coefficient estimation (all columns at once)
      coef <- R_inv %*% (t(Q) %*% Y_cct)

      # Vectorized residuals and sigma calculation
      resid <- Y_cct - X_mat %*% coef
      sigma <- sqrt(colSums(resid^2) / (TT - ncol(X_mat)))

      # Vectorized t-statistics and p-values for first 2 coefficients
      t_stats <- coef[1:2, ] / (se_factors[1:2] %o% sigma)
      p_values <- 2 * pt(abs(t_stats), TT - ncol(X_mat), lower.tail = FALSE)

      # Cauchy combination
      CCTa_pval = f_cauchypv(p_values[1, ])
      CCTd_pval = f_cauchypv(p_values[2, ])

      # Portfolio efficiency with optimized matrix ops
      mu <- matrix(colMeans(R))
      ones <- matrix(1, K + N)
      V <- cov(R)
      iV <- solve(V)

      a <- crossprod(mu, iV) %*% mu
      b <- crossprod(mu, iV) %*% ones
      c_val <- crossprod(ones, iV) %*% ones
      d_val <- a * c_val - crossprod(b)

      V11 <- cov(R1)
      iV11 <- solve(V11)
      mu1 <- matrix(colMeans(R1))
      ones1 <- matrix(1, K)

      a1 <- crossprod(mu1, iV11) %*% mu1
      b1 <- crossprod(mu1, iV11) %*% ones1
      c1 <- crossprod(ones1, iV11) %*% ones1
      d1 <- a1 * c1 - crossprod(b1)

      # Final calculations
      U <- (c1 + d1) / (c_val + d_val)
      HK_stat <- ((TT - K - N) / N) * (1 / sqrt(U) - 1)
      F1_stat <- ((TT - K - N) / N) * ((a - a1) / (1 + a1))
      F2_stat <- ((TT - K - N + 1) / N) * (((c_val + d_val) / (c1 + d1)) * ((1 + a1) / (1 + a)) - 1)

      # P-values
      HK_pval = pf(HK_stat, 2 * N, 2 * (TT - K - N), lower.tail = FALSE) # H0 alpha et delta =0
      F1_pval = pf(F1_stat, N, (TT - K - N), lower.tail = FALSE) # H0 alpha = 0
      F2_pval = pf(F2_stat, N, (TT - K - N + 1), lower.tail = FALSE) # H0 delta = 0

      ### !!!FIXME: Why do we return outpv and perform cauchy transformation in this function ?

      CCTad_pval <- f_cauchypv(c(CCTa_pval, CCTd_pval))

      # Faire que des outputs comme ça
      out <- list(HK = list(pval = HK_pval, stat = HK_stat, H0 = "alpha = 0 and delta = 0"),
                  F1 = list(pval = F1_pval, stat = F1_stat, H0 = "alpha = 0"),
                  F2 = list(pval = F2_pval, stat = F2_stat, H0 = "delta = 0"),
                  CCTa = list(pval = CCTa_pval, H0 = "alpha = 0"),
                  CCTd = list(pval = CCTd_pval, H0 = "delta = 0"),
                  CCTad = list(pval = CCTad_pval, H0 = "alpha = 0 and delta = 0")
                  )

      return(out)
}

f_alternative_tests <- function(R1, R2) {
  # R1: TT×p   bench returns
  # R2: TT×p2  test returns

  X <- R2
  Y <- R1

  Rbig <- cbind(R1, R2)
  TT   <- nrow(Rbig)
  p    <- ncol(R1)
  p2   <- ncol(R2)

  # Dependent y = first column of Rbig
  y <- Rbig[,1]

  # Build the FULL difference matrix: R[,1] - R[,-1]
  # this will be TT × (p-1 + p2)
  Diff <- sweep(Rbig[,-1, drop=FALSE], 1, y, FUN = function(x,y) y - x)

  # 1) GMVP regression: intercept + *all* differences
  X1    <- cbind(1, Diff)                # TT × (1 + p-1 + p2)
  k1    <- ncol(X1)
  XtX1  <- crossprod(X1)
  XtX1i <- solve(XtX1)
  theta1 <- XtX1i %*% crossprod(X1, y)   # (1 + p-1 + p2) × 1
  resid1 <- y - X1 %*% theta1
  sigma2_1 <- drop(crossprod(resid1)/(TT - k1))

  # Contrast to pick *only* the last p2* coefficients of theta1:
  # offset = 1 + (p-1)  (i.e. after intercept + bench‐bench diffs)
  offset <- 1 + (p-1)
  C1 <- cbind(
    matrix(0, p2, offset),    # zeros for the first offset columns
    diag(p2)                  # identity on the last p2 columns
  )                           # p2 × (1 + p-1 + p2)

  num1 <- as.numeric(
    crossprod(C1 %*% theta1,
              solve(C1 %*% XtX1i %*% t(C1)) %*% (C1 %*% theta1)
    ) / p2
  )
  F_gmvp <- num1 / sigma2_1
  p_gmvp <- pf(F_gmvp, p2, TT - k1, lower.tail = FALSE)


  # 2) Tangency:  y2 = constant 1 regressed on *the same* differences, no intercept
  y2    <- rep(1, TT)
  X2    <- Diff                        # TT × (p-1 + p2)
  XtX2  <- crossprod(X2)
  XtX2i <- solve(XtX2)
  gamma2 <- XtX2i %*% crossprod(X2, y2) # (p-1+p2) × 1
  resid2 <- y2 - X2 %*% gamma2
  sigma2_2 <- drop(crossprod(resid2)/(TT - ncol(X2)))

  # Contrast to pick the *last* p2 columns of gamma2
  # gamma2 has length (p-1 + p2)
  offset2 <- (p-1)
  C2 <- cbind(
    matrix(0, p2, offset2),
    diag(p2)
  )                                     # p2 × (p-1 + p2)

  num2 <- as.numeric(
    crossprod(C2 %*% gamma2,
              solve(C2 %*% XtX2i %*% t(C2)) %*% (C2 %*% gamma2)
    ) / p2
  )
  F_tang <- num2 / sigma2_2
  p_tang <- pf(F_tang, p2, TT - ncol(X2), lower.tail = FALSE)

  out <- list(
    BJ = list(pval = p_gmvp, stat = F_gmvp, H0 = "alpha = 0"),
    KM = list(pval = p_tang, stat = F_tang, H0 = "delta = 0")
  )

  return(out)

}

#' Portfolio Spanning Tests: Gungor-Luger, GRS, HK, and PY
#'
#' Performs a comprehensive suite of portfolio spanning tests to evaluate whether the asset returns in \code{R2}
#' are mean-variance spanned by the benchmark returns in \code{R1}. This includes:
#'
#' \itemize{
#'   \item \strong{GRS test} — Gibbons, Ross, and Shanken (1989): tests \eqn{\alpha = 0} using the F-distribution.
#'   \item \strong{HK test} — Huberman and Kandel (1987): joint test of \eqn{\alpha = 0} and \eqn{\delta = 0}.
#'   \item \strong{PY test} — Pesaran and Yamagata (2023): tests cross-sectional dependence in intercepts.
#'   \item \strong{GL test} — Gungor and Luger (2016): subseries-based max-F test using Monte Carlo (LMC) and bootstrap (BMC) approximations.
#' }
#'
#' The function returns p-values, test statistics, and decisions based on a threshold (default 0.05),
#' and supports both exact and fast (MC-only) modes via the \code{control} parameter.
#'
#' @param R1 A numeric matrix of benchmark returns (T x N).
#' @param R2 A numeric matrix of test asset returns (T x K).
#' @param control A named list of control parameters. Possible fields:
#' \describe{
#'   \item{\code{seed}}{Random seed (default: 1).}
#'   \item{\code{totsim}}{Number of simulations (default: 500).}
#'   \item{\code{do_trace}}{Whether to print results (default: \code{TRUE}).}
#'   \item{\code{pval_thresh}}{Significance threshold for decision (default: 0.05).}
#'   \item{\code{do_fast}}{If \code{TRUE}, skips GRS, PY, and HK tests to save time (default: \code{FALSE}).}
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{GL}}{Results from the Gungor & Luger Fmax test: LMC/BMC p-values, statistic, decision, and hypothesis.}
#'   \item{\code{GRS}}{List with statistic, p-value, and null hypothesis for GRS test (if applicable).}
#'   \item{\code{HK}}{Output from the Huberman-Kandel test (via \code{f_hk()}).}
#'   \item{\code{PY}}{List with statistic, p-value, and null hypothesis for the Pesaran-Yamagata test.}
#' }
#'
#' !!! This function MUST not be used by package users !!!
#' @keywords internal
#' @references
#' Gibbons, M. R., Ross, S. A., & Shanken, J. (1989). "A test of the efficiency of a given portfolio." \emph{Econometrica}, 57(5), 1121–1152.
#' Huberman, G., & Kandel, S. (1987). "Mean-Variance Spanning." \emph{Journal of Finance}, 42(4), 873–888.
#' Kan, R., & Zhou, G. (2012). "Tests of mean–variance spanning." \emph{Annals of Economics and Finance}, 13(1), 145–193.
#' Gungor, S., & Luger, R. (2016). "A New Test of Mean–Variance Spanning." \emph{Econometric Reviews}, 35(6), 1146–1171.
#' Pesaran, M. H., & Yamagata, T. (2023). "On Testing Cross-Sectional Dependence of the Intercepts." \emph{Working Paper}.
#'
#' @keywords internal
#'
#' @noRd
#'
f_span_gl_a <- function(R1, R2, control = list()) {

      con <- list(seed = 1, totsim = 500, do_trace = TRUE, pval_thresh = 0.05, do_fast = FALSE)
      con[names(control)] <- control
      thresh <- con$pval_thresh

      X <- R1
      Y <- R2

      K <- ncol(X)
      N <- ncol(Y)
      TT <- nrow(X)
      totsim <- con$totsim
      ones <- matrix(1, TT, 1)

      XX <- matrix(cbind(ones, X), TT, (K + 1))

      # OLS parameter estimates
      XX_crossprod <- crossprod(XX)
      Xtemp <- solve(XX_crossprod)
      Bhat1 <- Xtemp %*% crossprod(XX, Y)

      Ehat1 <- Y - XX %*% Bhat1         # Unrestricted
      SSRu <- crossprod(Ehat1)
      SigmaU <- SSRu / TT

      H <- matrix(0, 1, K + 1)
      H[1, 1] <- 1
      C <- matrix(0, 1, N)

      # Precompute reusable expressions
      HXt <- H %*% Xtemp
      HXtHt <- HXt %*% t(H)
      HXtHt_inv <- solve(HXtHt)
      HB_minus_C <- H %*% Bhat1 - C

      Bhat0 <- Bhat1 - Xtemp %*% t(H) %*% HXtHt_inv %*% HB_minus_C  # Restricted
      Ehat0 <- Y - XX %*% Bhat0
      SSRr <- crossprod(Ehat0)
      SigmaR <- SSRr / TT

      # Initialize output variables
      GRS_pval <- HK <- PY_pval <- NA

      # Skip GRS and Pesaran-Yamagata tests if do_fast = TRUE
      if (!con$do_fast) {
            # GRS test
            if ((TT - K - N) >= 1) {
                  ahat <- Bhat1[1, ]
                  shat <- SigmaU

                  mufactor <- matrix(colMeans(X), K, 1)
                  X_centered <- scale(X, center = TRUE, scale = FALSE)
                  ssqm <- crossprod(X_centered) / TT

                  mufactor_ssqm <- solve(ssqm, mufactor)
                  GRS_denom <- 1 + crossprod(mufactor, mufactor_ssqm)
                  GRS_num <- ((TT - N - K) / N) * crossprod(ahat, solve(shat, ahat))
                  GRS_stat <- GRS_num / GRS_denom
                  GRS_pval <- 1 - pf(GRS_stat, N, (TT - N - K))
                  HK <- f_hk(X, Y)
            }

            # Pesaran-Yamagata tests
            v <- TT - K - 1
            X_crossprod <- crossprod(X)
            X_crossprod_inv <- solve(X_crossprod)
            X_ones <- X %*% X_crossprod_inv %*% crossprod(X, ones)
            MX_ones <- ones - X_ones

            num_scalar <- crossprod(ones, MX_ones)[1, 1] * v
            num <- matrix(num_scalar, N, 1)

            diag_SigmaU <- diag(SigmaU)
            t2 <- (Bhat1[1, ]^2) * num / (TT * diag_SigmaU)

            pN <- 0.05 / (N - 1)
            thetaN <- qnorm(1 - pN / 2)^2
            rhobar <- 0

            for (i in 2:N) {
                  for (j in 1:(i - 1)) {
                        temp <- SigmaU[i, j] / sqrt(SigmaU[i, i] * SigmaU[j, j])
                        temp2 <- temp^2
                        if (v * temp2 >= thetaN) {
                              rhobar <- rhobar + temp2
                        }
                  }
            }
            rhobar <- rhobar * 2 / (N * (N - 1))

            Jalpha2 <- sum(t2 - v / (v - 2)) / sqrt(N)
            den <- (v / (v - 2)) * sqrt(2 * (v - 1) * (1 + (N - 1) * rhobar) / (v - 4))
            PY_stat <- Jalpha2 / den
            PY_pval <- 1 - pnorm(PY_stat)
      }

      ## MC Fmax tests (always computed)
      diag_SSRr <- diag(SSRr)
      diag_SSRu <- diag(SSRu)
      temp <- (diag_SSRr - diag_SSRu) / diag_SSRu
      Fmax <- max(temp)
      Fmax_actual <- Fmax

      LMCstats <- matrix(0, totsim, 1)
      LMCstats[totsim] <- Fmax
      BMCstats <- matrix(0, totsim, 1)
      BMCstats[totsim, 1] <- Fmax

      Ehat0data <- Ehat0
      Bhat0data <- Bhat0

      premult <- Xtemp %*% t(H) %*% HXtHt_inv
      Xtemp_XX <- Xtemp %*% t(XX)
      set.seed(123456 * con$seed)

      sim_count <- totsim - 1
      sign_mat <- matrix(sign(rnorm(TT * sim_count)), TT, sim_count)
      esim_array <- array(Ehat0data, dim = c(TT, N, sim_count)) *
            aperm(array(sign_mat, dim = c(TT, sim_count, N)), c(1, 3, 2))

      Ysim_array <- array(XX %*% Bhat0data, dim = c(TT, N, sim_count)) + esim_array

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
      LMCstats[1:sim_count] <- matrix(temp_LMC, ncol = N, byrow = TRUE) %>% apply(1, max)

      SSRr_BMC_vec <- colSums(matrix(esim_array^2, TT, N * sim_count))
      temp_BMC <- (SSRr_BMC_vec - SSRu_vec) / SSRu_vec
      BMCstats[1:sim_count] <- matrix(temp_BMC, ncol = N, byrow = TRUE) %>% apply(1, max)

      uu <- runif(totsim)

      temp <- f_ranklex(LMCstats, uu)
      GL_pval_LMC <- (totsim - temp + 1) / totsim

      temp <- f_ranklex(BMCstats, uu)
      GL_pval_BMC <- (totsim - temp + 1) / totsim

      if (con$do_trace) {
            print('============================================')
            Decisions <- Decisions_string <- ""
            print('F-Max')
            print(c('------', 'F-max           =', Fmax_actual))

            if (GL_pval_LMC > thresh) {
                  Decisions <- 1
                  Decisions_string <- "Accept"
                  print(c('ACCEPT ----------', 'LMCp-value       =', GL_pval_LMC))
            }

            if (GL_pval_BMC <= thresh) {
                  Decisions <- 0
                  Decisions_string <- "Reject"
                  print(c('REJECT ----------', 'BMC p-value      =', GL_pval_BMC))
            }

            if ((GL_pval_LMC <= thresh) & (GL_pval_BMC > thresh)) {
                  Decisions <- NA
                  Decisions_string <- "Inconclusive"
                  print(c('INCONCLUSIVE ----', 'LMC p-value      =', GL_pval_LMC))
                  print(c('INCONCLUSIVE ----', 'BMC p-value      =', GL_pval_BMC))
            }

      } else {
            Decisions <- Decisions_string <- ""

            if (GL_pval_LMC > thresh) {
                  Decisions <- 1
                  Decisions_string <- "Accept"
            }

            if (GL_pval_BMC <= thresh) {
                  Decisions <- 0
                  Decisions_string <- "Reject"
            }

            if ((GL_pval_LMC <= thresh) & (GL_pval_BMC > thresh)) {
                  Decisions <- NA
                  Decisions_string <- "Inconclusive"
            }

      }

      if (con$do_fast) {
        PY <- GRS <- HK <- NA

      } else {
        # Pesaran Yamagata
        PY <- list(pval = PY_pval, stat = PY_stat, H0 = "alpha = 0")

        if (TT - K - N >= 1) {
          # Gibbons et al.
          GRS <- list(pval = GRS_pval, stat = GRS_stat, H0 = "alpha = 0")

          # Huberman & Kandel (already formatted in the f_hk)
          HK <- HK
        } else {
          GRS <- HK <- NA
        }
      }

      # Gungor & Luger
      GL <- list(
        pval_LMC = GL_pval_LMC,
        pval_BMC = GL_pval_BMC,
        stat      = Fmax_actual,
        Decisions = Decisions,
        Decisions_string = Decisions_string,
        H0        = "alpha = 0"
      )

      if (con$do_fast || (TT - K - N < 1)) {

        out <- list(PY = PY, GRS = NA, HK = NA, GL = GL, BJ = NA, KM = NA)

        return(out)
      }

      alt <- f_alternative_tests(R1, R2)
      BJ <- alt$BJ
      KM <- alt$KM

      out <- list(PY = PY, GRS = GRS, HK = HK, GL = GL, BJ = BJ, KM = KM)

      return(out)
}

#' Gungor–Luger Spanning Test (Joint Null: Alpha = 0 and Delta = 0)
#'
#' Performs the spanning test proposed by Gungor and Luger (2016), targeting the joint null hypothesis
#' that both the intercepts and slope shifts of a set of test assets with respect to benchmark assets are zero:
#' \eqn{H_0: \alpha = 0 \text{ and } \delta = 0}.
#'
#' This implementation uses the F-max statistic across individual assets and computes p-values using both
#' Monte Carlo (LMC) and bootstrap (BMC) resampling strategies. The test accounts for joint deviations
#' from the benchmark space and supports reproducible simulation-based inference.
#'
#' @param R1 A numeric matrix of benchmark returns (T x N).
#' @param R2 A numeric matrix of test asset returns (T x K).
#' @param control An optional named list of control parameters:
#' \describe{
#'   \item{\code{seed}}{Integer seed for reproducibility (default: \code{1}).}
#'   \item{\code{totsim}}{Number of simulations for p-value estimation (default: \code{500}).}
#'   \item{\code{pval_thresh}}{Threshold for determining statistical decision (default: \code{0.05}).}
#'   \item{\code{do_trace}}{Logical; if \code{TRUE}, prints intermediate outputs and decision (default: \code{TRUE}).}
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{pval_LMC}}{P-value computed via Monte Carlo simulation.}
#'   \item{\code{pval_BMC}}{P-value computed via bootstrap-based simulation.}
#'   \item{\code{stat}}{The observed F-max test statistic.}
#'   \item{\code{Decisions}}{Numeric decision: 1 = accept null, 0 = reject, \code{NA} = inconclusive.}
#'   \item{\code{Decisions_string}}{Character representation of the decision ("Accept", "Reject", "Inconclusive").}
#'   \item{\code{H0}}{A description of the tested null hypothesis.}
#' }
#'
#' !!! This function MUST not be used by package users !!!
#' @keywords internal
#'
#' @references
#' Gungor, S., & Luger, R. (2016). "A New Test of Mean–Variance Spanning." \emph{Econometric Reviews}, 35(6), 1146–1171.
#'
#' @keywords internal
#'
#' @noRd
#'
f_span_gl_ad <- function(R1, R2, control = list()) {

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
      LMCstats[1:(totsim - 1)] <- matrix(temp_LMC, ncol = N, byrow = TRUE) %>% apply(1, max)

      # BMC statistics
      SSRr_BMC_vec <- colSums(matrix(esim_array^2, TT, N * (totsim - 1)))
      temp_BMC <- (SSRr_BMC_vec - SSRu_vec) / SSRu_vec
      BMCstats[1:(totsim - 1)] <- matrix(temp_BMC, ncol = N, byrow = TRUE) %>% apply(1, max)

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
