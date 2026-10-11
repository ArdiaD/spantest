# Changes in Version 1.4-4 (DA)
- The Cauchy combination replaces a p-value equal to one by 1 - 1e-12, as the paper
  states (Section 3.1). An exact one, a null event for continuous data, otherwise
  contributes tan(-pi/2) = -1.63e16 and outweighs any small p-value: p = (1e-14, 1)
  combined to 1, and now to 2e-14. Values below one are unchanged, so the result of
  any data set without an exact one is unchanged bit for bit.
- span_as() warns when some test assets have no p-value (NA) and are left out of the
  combination, with their number. They were dropped silently.
- span_gl_a(), span_gl_ad(): pval_LMC is the local Monte Carlo p-value and pval_BMC
  the bounds Monte Carlo p-value of Gungor and Luger (2016); the documentation called
  them "Least-Favorable" and "Balanced". stat is the largest asset-level ratio
  (SSR_r - SSR_u)/SSR_u, without the constant of the F form.
- span_simulate(): the documentation says that the preset numbers are not the paper's
  DGP numbers (the AR and AR-GARCH blocks are swapped), and that names reproduce it.

# Changes in Version 1.4-3 (DA)
- span_as(): at L > 0 the multiplier weights are now drawn per statistic: an
  independent T x N matrix for the alpha scores and another for the delta
  scores, each a product of L independent N(1, 1) factors, so the weights are
  independent across factors, dates, assets and the two scores. This is how the
  theory defines them (drawn afresh for each asset-level statistic). Up to 1.4-2
  one T-vector was shared by every asset and by both scores.
  A shared vector leaves the cross-sectional dependence of the asset-level
  statistics intact; independent weights shrink it (by about 2^-L for scores
  without autocorrelation). Under a pervasive common factor in the residuals
  (span_simulate(gamma = ), gamma_i ~ U(0.7, 0.9)), at L = 2, T = 250, K = 2 and
  10, N = 100 and 1000, 1000 replications, the three tests reject a true null in
  6.0-8.5% of the samples with the shared vector and in 2.4-4.8% with
  independent weights (5.4-9.6% at L = 0). Under approximately sparse dependence
  the two differ by about one point of size; alpha and joint power gain 2-7
  points and delta power at N = 1000 loses 3-5.
  Results at L = 0 are unchanged bit for bit (checked against 1.4-2 on 60 data
  sets); every result at L > 0 changes, including with the same seed. The draw
  of a seed puts the alpha factors first, then the delta factors.
  Independent weights do not make the global p-value less dependent on the
  draw: a Cauchy average of independent terms does not concentrate as N grows.
  Cost: T x N x 2L normal draws per call. span_as() at T = 250, N = 1000,
  three exponents and L = 0, 1, 2 takes 0.12 s instead of 0.08 s.
- span_as(): the control B is removed: a call is one draw of the weights, from
  seed (default 123), and passing B is an error. B > 1 merged B draws per asset
  with the Cauchy rule, which treats B views of one asset's returns as B
  independent pieces of evidence; in the one design where it was measured, it
  was mildly liberal per asset.
  Results with B = 1, the default, are unchanged: that draw already came from
  seed. f_getpv_batch() and f_getpv() lose their B argument too.
- span_hk(), span_f1(), span_f2(): the covariance matrices are now the
  maximum-likelihood ones (divisor T, new internal f_cov_ml()) instead of cov()
  (divisor T - 1), as in Kan and Zhou (2012) and in the papers that use these
  tests. span_hk() is then the likelihood-ratio form of the Huberman-Kandel test,
  and span_f1() returns the GRS statistic exactly, since F1, the first step of
  Kan and Zhou's step-down test, is the GRS test (checked to 3e-12; span_grs() and
  span_bj() were already identical). The statistics move by at most 1.4-1.7% at
  T = 60 to 250, and size by about 0.2 point. span_f2()'s documentation now says
  that it tests delta = 0 given alpha = 0, and points to span_km() for delta = 0
  with free intercepts.
- span_py(): the condition N <= T - K - 1 is removed. The statistic uses the
  asset-level t-statistics and the pairwise residual correlations, never the
  inverse of the N x N residual covariance, and Pesaran and Yamagata design the
  test for large N, N > T included. The thresholded correlation sum is
  vectorised (identical results where the test was already defined; 0.14 s at
  N = 1000, T = 250).
- span_mstv(): the derandomized rule now follows the authors' replication code
  (github.com/PierValla95/RandomizedAlphaTest, Empirics/Analysis_MSTV_nu4.R): B =
  floor(log(N)^2) perturbation draws from seed, after the one-shot draw, each a
  block of N normals; Q is the share of draws whose maximum does not exceed the
  critical value, and the rule rejects when Q < (1 - tau) - B^(-1/4). Up to 1.4-2
  Q was its B -> infinity limit in closed form, prod_i Phi(c - psi_i). That limit
  is a different rule in practice: with B between 15 and 50 the simulated share is
  coarse, and a size-adjusted MSTV-D built on the closed form reaches 99-100%
  power where the authors' rule reaches 18-73% (DGP12, K <= 10, a = 0.2). With
  identical draws, Q and the decision equal those of the co-author's
  implementation of the authors' rule (40 data sets). logQ is no longer returned.
- span_mstv(): the default nu is 5, the value of the authors' simulations,
  instead of 4. The statistic, its pooled scale (divisor N T, as in their code)
  and the critical value are unchanged.
- span_simulate(): the benchmarks R1 are now the simulated process itself,
  as in the paper's equation (Toeplitz correlation rho_factor, the AR/GARCH
  dynamics and the innovation law of the preset), and the test assets load
  on them with B_1j = (1 + ncp)(2 - K) and B_kj = 1 + ncp for k >= 2, so that
  alpha_j = ncp and delta_j = -ncp as before. Up to 1.4-2 the benchmarks were
  [z1, z_{-1} + z1] with loadings 1 + ncp on z: the same loadings on the
  benchmarks, but benchmarks with correlations of 0.91-0.95 and variances of
  about 1, 3.6 and 3.3 instead of the Toeplitz structure. The random draws are
  unchanged. On paired panels the alpha tests return the same p-values (GRS,
  PY and the SCT of alpha at L = 0 and 2, to 1e-13), while the delta and joint
  tests change; KM is 1.5-2.5 points more oversized at K <= 10 with the
  process as benchmarks (David's decision, so that the code is the equation).
- f_mult() replaces the internal f_prods(). f_getpv() gains wN and wcol so
  that a loop over the columns reproduces f_getpv_batch() (tested at L = 0, 1
  and 2).
- Tests: test-span-as-weights.R checks the structure and moments of the
  weights, the caller's random-number stream, two identical assets, the L = 0
  results against 1.4-2, and that at L = 2 each test fires on its own
  hypothesis only. The tests of B give way to tests of the seed: the default,
  its effect at L = 0 (none) and L = 2, its validation, and the error on B.

# Changes in Version 1.4-2 (DA)
- span_mstv(): new test, the randomized alpha test of Massacci, Sarno, Trapani
  and Vallarino (2026, forthcoming in JASA). It estimates equation by equation,
  needs no covariance matrix and allows N to grow faster than T, so it is the
  natural competitor to span_as() on the alpha side when N >> T. Both their
  one-shot test (Theorem 3.1, with a p-value from the Gumbel limit) and their
  derandomized decision rule (Section 3.2) are returned.
  The derandomized quantity is computed in CLOSED FORM rather than by simulating
  B perturbations: the draws are independent of each other and of the data, so
  conditionally on the sample the share of replications that do not reject is
  exactly prod_i Phi(c_tau - psi_i). That is the B -> infinity limit of their
  average, so it removes the simulation error and the dependence on the draws
  while leaving the rule untouched, and it costs one pass over the cross-section
  instead of B regressions. Accumulated on the log scale, since the product runs
  over N terms and underflows long before its logarithm stops being informative.
  Returns NA below N = 10: the reference distribution is the Gumbel limit of the
  maximum of N normals, and at N = 2 the one-shot test rejects 26% of the time
  at a nominal 5%.
  Checked against the authors' own replication code: bit-identical statistics
  over 300 replications of their data-generating process, and their published
  rejection frequencies are reproduced (their Table 4.1 power, 93.5% published
  against 93.3% here).
- span_simulate(): new `gamma` argument, a pervasive common component in the
  idiosyncratic terms (eps_it + gamma_i g_t). The existing `rho_error` gives a
  Toeplitz decay, which is approximately sparse and inside the assumptions of
  the tests; `gamma` does not die out across the cross-section and is there to
  leave that regime, which is where the randomized alpha test is designed to
  operate. NULL (the default) draws no random numbers and the draw is placed
  after the K + N existing ones, so every earlier result is reproduced bit for
  bit; a regression test pins this down across four processes and three (K, N)
  shapes including K = 100, N = 1000.
- span_mstv(), span_as(): new section 'Simulations'. Both tests are randomized
  and both take a fixed default seed, so a Monte Carlo study that keeps the
  default measures them conditionally on one draw, shared by every replication.
  For span_mstv() this is fatal: at small K the one-shot statistic is
  essentially the maximum of the perturbation vector, and with seed 123 at
  N = 1000 that maximum is 3.24 against a critical value of 3.98, so the
  measured size is 0 where a new draw per replication gives about 4%. For
  span_as() the weights enter through subseries sums and the effect is about one
  point of size either way. Documentation only; no result of the package changes.

# Changes in Version 1.4-1 (DA)
- span_as(): the documentation of `ks` said `floor(T^k)` was the subseries SIZE
  and that the subseries overlapped. Both were the reverse of the code, which
  matches the paper: `floor(T^k)` is the NUMBER of blocks, each of length about
  T / floor(T^k), forming an exact partition of 1:T. No result changes; the
  inverted reading had however already propagated into the manuscript's own
  explanation of its sensitivity appendix.
- span_simulate(): the twelve presets are now a clean 3 x 4 factorial -- one
  innovation law per family (normal; standardised t with 5 df; standardised
  skew-t with 4 df and xi = 0.9) repeated across the four dynamics -- so moving
  down a column changes the DYNAMICS and nothing else. It did not before: the two
  Student presets without GARCH were RAW t_5, of variance 5/3, while the two with
  GARCH were standardised t_4. Standardisation is required for the GARCH
  recursion, since omega/(1 - alpha - beta) is the unconditional variance only
  for a unit-variance innovation, and it had been applied where it was needed and
  omitted where it merely mattered. The consequence was that rows within the
  Student family of a size table were not comparable with one another.
  RESULTS: iid-ST and AR-ST are unaffected -- their change is a pure rescaling and
  the tests are scale-invariant, verified to 2e-13 across 27 columns and 40 draws.
  GARCH-ST and AR-GARCH-ST do change, because their degrees of freedom move from
  4 to 5; anything computed from those two presets needs regenerating.
- span_simulate(): new argument `scale`, and the twelve presets now carry unit
  PROCESS variance rather than merely unit innovations. Standardising the
  innovations is not enough: an AR(1) with unit innovations has variance
  1/(1 - phi^2) = 1.042, so the AR and AR-GARCH rows carried 4% more noise than
  the iid ones, and the twelve differed in level as well as in the shape of the
  dependence. The GARCH parameters had always been chosen so that
  omega/(1 - alpha - beta) is exactly 1 -- the design already normalised the
  process where someone had thought about it -- and this completes that intent.
  `scale = "innovation"` remains the default outside the presets, so a bare
  span_simulate(dynamics = "ar", ar = 0.2) is still the textbook object; the
  presets set `scale = "process"`. Rescaling leaves the autocorrelation function
  untouched, so `ar` keeps its meaning (measured: 0.202 either way).
  RESULTS: no size cell changes. At ncp = 0 this rescales the whole system and
  the tests are scale-invariant -- verified across the six AR/AR-GARCH presets,
  30 p-value columns and 25 draws, maximum discrepancy 2e-12. Only the
  signal-to-noise ratio under the alternative moves, by 2.1%, which on an ncp
  grid stepped by 0.1 is well inside the Monte Carlo error of a 500-replication
  rejection rate.
- span_dgp_table() reports `scale` as well, so its columns are enough to rebuild
  a preset by hand. A table that describes a process only partly is how the
  previous discrepancy survived: the reconstruction silently took a default the
  preset does not use.
- span_simulate(): the preset definitions are no longer written out twice. They
  had been restated as a list inside span_simulate() and again inside
  span_dgp_table(), and the two had drifted -- which is how the discrepancy above
  survived. There is now one table, .DGP_PRESETS, and span_simulate() reads it.
  tests/testthat/test-span_simulate.R checks that each innovation family shares
  one df and one standardisation across all four dynamics, and that every
  heavy-tailed preset really has unit variance rather than merely being flagged
  as standardised.

# Changes in Version 1.4-0 (DA)
- span_as(): new control parameters `B` and `seed`, de-randomising the L > 0
  test. For L > 0 the score is multiplied by a random weight vector that is
  SHARED by every test asset, so its effect does not average out across the
  cross-section: a single draw is a single realisation of the test. The test is
  valid for any fixed draw (its size is correct), but the p-value returned on one
  data set can move by orders of magnitude from draw to draw -- most visibly when
  T is short, where floor(T^(1/3)) leaves few subseries. `B > 1` draws B
  independent weight vectors and merges the per-asset p-values with the Cauchy
  rule, which is valid under arbitrary dependence and is not carried by a single
  extreme draw. `B = 1` (the default) reproduces the previous behaviour exactly;
  `B = 100` or more is recommended for any reported application. The seed of the
  first draw is now user-visible (`seed`, default 123, draw b uses seed + b - 1)
  instead of hard-coded. Consecutive seeds were kept, rather than an explicit
  substream generator, so that `B = 1` reproduces the historical result exactly;
  the draws they produce are empirically indistinguishable from independent
  (mean absolute pairwise correlation 0.051 over 100 seeds at T = 250, against
  0.050 expected), and this is documented under "Choosing B".
- span_simulate(): `dgp` now accepts a process NAME as well as a preset number
  -- "AR-GARCH-SKST" instead of 9 -- and the two are exactly equivalent. The
  twelve names are exported as `DGP_NAMES` and described by the new
  `span_dgp_table()`, which lists each preset with its innovation law, dynamics
  and degrees of freedom. Numbers remain accepted, so existing code is
  unaffected -- including a preset handed over as a string, since "9" and 9
  select the same preset. Name matching is exact after trimming and
  case-folding, not partial: "AR" is a prefix of six of the twelve names.
  The motivation is that a preset number carries no meaning and can
  be transposed silently: presets 7-9 are the AR-GARCH family and 10-12 the AR
  family, an ordering that does not match every convention in use. A name
  cannot be transposed, and a wrong one is an error rather than a wrong process.
- span_simulate(): `dgp` must be a whole number. Previously `dgp = 9.7` passed
  through as.integer() and silently ran preset 9, a process the argument does
  not name. Whole numbers stored as doubles (`dgp = 9`) are unaffected.
- span_as(): `B` and `seed` must now be whole numbers. Previously a fractional
  `B = 1.9` passed validation and was silently truncated to a single draw.
- span_gl_a() / span_gl_ad(): a singular benchmark design (e.g. collinear
  benchmarks) now returns NA p-values instead of raising, matching the
  convention already used by span_grs(), span_f2() and span_py().
- The Cauchy combination now returns NA, not a silent NaN, when it is handed no
  usable p-values (every test asset missing).
- The subseries t-test now fails with an explicit message when the sample is too
  short to form two subseries, instead of dividing by zero degrees of freedom.

# Changes in Version 1.3-0 (DA)
- span_gl_a() / span_gl_ad(): the C++ sign-flip kernel now forms the restricted
  residuals via Ehat0s = Ehat1 + (XX %*% premult) %*% (H %*% Bhat1 - C), avoiding
  the per-simulation T x (K+1) x N restricted-least-squares matmul. This removes
  one of the three large per-simulation matrix products, cutting the GL kernel
  time by roughly 20% (about 1.2x at K = 100, N = 200, T = 250), with
  bit-identical p-values.
- span_gl_a() and span_gl_ad() now default to do_trace = FALSE (no console
  output unless requested).
- span_grs(), span_f2() and span_py() now return NA (like the other tests)
  instead of raising an error when a covariance/design matrix is singular
  (e.g. collinear benchmarks).
- span_simulate() validates its arguments (dimensions, correlation ranges,
  degrees of freedom, GARCH stationarity).
- span_as() is faster: the internal subseries t-test (f_ttest) no longer calls
  stats::t.test() once per column (whose input-checking dominated the runtime).
  The one-sample two-sided p-values are now computed directly and vectorised
  across columns, giving identical results (to floating-point) at about 2.5x
  lower cost for span_as().
- NEW: span_simulate() generates benchmark and test-asset returns for spanning
  size/power studies from a factor model with a controllable spanning violation
  (ncp). Innovations may be normal, Student-t, or skew-t, with iid, AR(1),
  GARCH(1,1), or AR-GARCH dynamics; the `dgp` argument selects the twelve
  processes of Ardia and Sessinou (2025). It is a fast, validated drop-in for
  ad-hoc fGarch::garchSim() loops (roughly 30-45x faster for the GARCH DGPs) and
  matches the reference processes in distribution, dynamics, and test size/power.
  The skew-t innovations use a base-R re-implementation of the Fernandez-Steel
  standardised skew-t (identical draws to fGarch::rsstd given the same RNG
  state), so the package has no external simulation dependency. Its GARCH(1,1) recursion runs
  in C++ (bit-for-bit identical to the R loop; floating-point contraction is
  disabled so the fused multiply-add does not alter the roundings), which about
  halves the data-generation time for the GARCH DGPs.
- span_gl_a() and span_gl_ad() now run their sign-flip Monte Carlo simulations
  in C++ (via Rcpp / RcppArmadillo). The random signs and the tie-breaking
  uniforms are still drawn in R, so the RNG stream -- and hence every p-value and
  decision -- is unchanged (verified identical to 1.2-1 across dimensions,
  totsim, and seeds); only the per-simulation SSR / F_max computation moved to a
  streaming C++ kernel that never forms the large T x (N*totsim) intermediates.
  About 3x faster than the 1.2-1 R implementation at large K and N, with much
  lower memory use, which also improves the parallel scaling of simulation
  studies. Installing from source now requires a C++ compiler (Rcpp and
  RcppArmadillo were added as build dependencies).

# Changes in Version 1.2-1 (DA)
- span_as() is now vectorized across the test cross-section: benchmark-only QR
  decompositions are formed once and the swap regressions are obtained by
  Frisch-Waugh partialling, replacing the per-asset loop of full QR
  factorizations. Output is numerically identical to 1.2-0 (verified); span_as
  is roughly 12-15x faster for moderate-to-large test sets.
- span_gl_a() and span_gl_ad() streamlined. The sign-flip simulations are now
  assembled directly as T x (N*totsim) matrices (recycling the residuals and
  fitted values across simulations, gathering the sign columns) instead of
  building a T x N x totsim array and transposing it with aperm(), which removed
  the dominant cost at large N. In addition the balanced-MC restricted SSR is
  constant across sign-flips (it equals the raw restricted SSR, since squaring
  removes the flipped sign) so it is computed once; redundant array/matrix
  reshapes in the constrained-estimate step were collapsed to matrix products;
  and unused Sigma allocations were removed. Output is bit-for-bit identical to
  1.2-0 (verified via identical()); GL is roughly 1.5-1.7x faster, most visibly
  for large test cross-sections (e.g. K = 100, N = 200).

# Changes in Version 1.2-0 (DA)
- NEW: span_as(), the Ardia-Sessinou subseries-based Cauchy Combination Test (CCT)
  for high-dimensional spanning, is now available and exported. It returns combined
  p-values for the alpha (CCTa), variance/slope (CCTd), and joint (CCTad) nulls,
  robust to serial/cross-sectional dependence and conditional heteroskedasticity.
  Internal helper f_getpv() implements the per-asset computation.
- span_bj now implements the Britten-Jones (1999) tangency-spanning test correctly
  (regression of ones on raw returns); the previous version had no power against
  alpha alternatives
- span_py now returns NA instead of erroring when there is a single test asset (N = 1)
- f_prods no longer overwrites the caller's global RNG state (.Random.seed is
  saved and restored)
- span_gl_ad trace now prints the decision label instead of the numeric code
- tail p-values in span_grs and span_py use lower.tail = FALSE for better precision
- added @importFrom Rdpack reprompt to silence the R CMD check NOTE on unused Imports
- documentation: benchmark/test dimensions in span_km, span_gl_a, and span_gl_ad
  aligned with the rest of the package (R1 is T x K, R2 is T x N)
- added correctness (size/power) tests for span_bj and an N = 1 test for span_py
- added tests for span_as and f_getpv

# Changes in Version 1.1-3 (DA)
- seed is finally removed from GL functions

# Changes in Version 1.1-2 (BS)
- DESCRIPTION modified according to CRAN guidelines
- seed is now an input for GL functions
- error fixed in the test file of span_bj

# Changes in Version 1.1-1 (BS)
- Function span_km and span_bj were testing the same

# Changes in Version 1.1-0 (DA)
- First release public version
