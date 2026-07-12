# Changes in Version 1.2-1 (DA)
- span_as() is now vectorized across the test cross-section: benchmark-only QR
  decompositions are formed once and the swap regressions are obtained by
  Frisch-Waugh partialling, replacing the per-asset loop of full QR
  factorizations. Output is numerically identical to 1.2-0 (verified); span_as
  is roughly 12-15x faster for moderate-to-large test sets.

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
