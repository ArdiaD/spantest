#include <RcppArmadillo.h>
using namespace arma;

// Simulated F-max statistics for the Gungor-Luger (2016) sign-flip Monte Carlo
// spanning test. The random sign-flips and the tie-breaking uniforms are drawn
// in R (so the RNG stream is unchanged); this routine only performs the per-
// simulation restricted/unrestricted SSR and F-max computation. Each simulation
// is processed on its own T x N buffers, so the large T x (N * nsim)
// intermediates of the fully-vectorised R implementation are never formed.
//
// Inputs (all pre-computed in R):
//   XX        T x (K+1)   design matrix [1, benchmarks]
//   Xtemp_XX  (K+1) x T   (X'X)^{-1} X'
//   XXB0      T x N       XX %*% Bhat0 (restricted fitted values)
//   Ehat0     T x N       restricted residuals
//   H         nH x (K+1)  restriction matrix (nH = 1 for alpha, 2 for joint)
//   C         nH x N      restriction RHS
//   premult   (K+1) x nH  (X'X)^{-1} H' (H (X'X)^{-1} H')^{-1}
//   sign_mat  T x nsim    +/-1 sign-flips
//   ssrr_bmc  1 x N       colSums(Ehat0^2) (balanced-MC restricted SSR, constant)
//
// [[Rcpp::export]]
Rcpp::List gl_sim_stats(const arma::mat& XX,
                        const arma::mat& Xtemp_XX,
                        const arma::mat& XXB0,
                        const arma::mat& Ehat0,
                        const arma::mat& H,
                        const arma::mat& C,
                        const arma::mat& premult,
                        const arma::mat& sign_mat,
                        const arma::rowvec& ssrr_bmc) {
  const uword nsim = sign_mat.n_cols;
  vec LMC(nsim), BMC(nsim);
  for (uword s = 0; s < nsim; ++s) {
    mat esim = Ehat0;
    esim.each_col() %= sign_mat.col(s);             // esim[t,n] = Ehat0[t,n]*sign[t,s]
    mat Ysim  = XXB0 + esim;                         // T x N
    mat Bhat1 = Xtemp_XX * Ysim;                     // (K+1) x N  (unrestricted)
    mat Ehat1 = Ysim - XX * Bhat1;                   // T x N
    rowvec SSRu = sum(Ehat1 % Ehat1, 0);             // 1 x N
    mat Bhat0s = Bhat1 - premult * (H * Bhat1 - C);  // (K+1) x N  (restricted)
    mat Ehat0s = Ysim - XX * Bhat0s;                 // T x N
    rowvec SSRr = sum(Ehat0s % Ehat0s, 0);           // 1 x N
    LMC[s] = ((SSRr     - SSRu) / SSRu).max();
    BMC[s] = ((ssrr_bmc - SSRu) / SSRu).max();
  }
  return Rcpp::List::create(Rcpp::Named("LMC") = LMC,
                            Rcpp::Named("BMC") = BMC);
}
