test_that("B = 1 reproduces the historical single-draw behaviour", {
  set.seed(7)
  x <- matrix(rnorm(200 * 3), 200, 3)
  y <- matrix(rnorm(200 * 6), 200, 6)

  # the default must stay exactly what it was before B/seed existed
  a <- span_as(x, y)
  b <- span_as(x, y, control = list(B = 1L, seed = 123L))
  expect_identical(a, b)

  pv1 <- spantest:::f_getpv_batch(x, y, ks = 1/3, L = c(0, 2))
  pv2 <- spantest:::f_getpv_batch(x, y, ks = 1/3, L = c(0, 2), B = 1L, seed = 123L)
  expect_identical(pv1, pv2)
})

test_that("L = 0 is unaffected by B (there is no perturbation to merge)", {
  set.seed(8)
  x <- matrix(rnorm(150 * 2), 150, 2)
  y <- matrix(rnorm(150 * 4), 150, 4)

  p1 <- span_as(x, y, control = list(L = 0, B = 1L))
  p9 <- span_as(x, y, control = list(L = 0, B = 9L))
  expect_identical(p1, p9)
})

test_that("B > 1 merges draws and stabilises the L > 0 p-value", {
  set.seed(9)
  x <- matrix(rnorm(150 * 2), 150, 2)
  y <- matrix(rnorm(150 * 8), 150, 8)

  merged <- span_as(x, y, control = list(L = 2, B = 20L))
  expect_true(all(vapply(merged, function(p) is.finite(p) && p > 0 && p <= 1, logical(1))))

  # the merged value must not be pinned to any single draw
  draws <- vapply(1:20, function(s)
    span_as(x, y, control = list(L = 2, B = 1L, seed = s))$CCTad_L2_k1, numeric(1))
  expect_gt(stats::sd(draws), 0)                       # draws genuinely differ
  expect_gte(merged$CCTad_L2_k1, min(draws))           # merged is not the luckiest draw
})

test_that("seed is honoured and B is validated", {
  set.seed(10)
  x <- matrix(rnorm(120 * 2), 120, 2)
  y <- matrix(rnorm(120 * 3), 120, 3)

  expect_false(isTRUE(all.equal(span_as(x, y, control = list(L = 2, seed = 1L))$CCTad_L2_k1,
                                span_as(x, y, control = list(L = 2, seed = 999L))$CCTad_L2_k1)))
  expect_error(span_as(x, y, control = list(B = 0)))
  expect_error(span_as(x, y, control = list(B = NA)))
})

test_that("control validation rejects non-integer B and seed", {
  set.seed(12)
  x <- matrix(rnorm(120 * 2), 120, 2); y <- matrix(rnorm(120 * 3), 120, 3)
  expect_error(span_as(x, y, control = list(B = 1.9)))      # would silently truncate to 1
  expect_error(span_as(x, y, control = list(seed = 1.5)))
  expect_error(span_as(x, y, control = list(B = c(1, 2))))
  expect_silent(span_as(x, y, control = list(B = 3L, seed = -7L)))
})

test_that("degenerate inputs return NA rather than NaN or an error", {
  expect_true(is.na(spantest:::f_cauchypv(numeric(0))))     # every asset missing
  expect_false(is.nan(spantest:::f_cauchypv(numeric(0))))
  # singular benchmark: GL now matches the package's NA convention
  set.seed(13)
  x <- matrix(rnorm(100), 100, 1); x <- cbind(x, x)         # collinear
  y <- matrix(rnorm(300), 100, 3)
  expect_true(is.na(span_gl_a(x, y, control = list(totsim = 20))$pval_LMC))
  expect_true(is.na(span_gl_ad(x, y, control = list(totsim = 20))$pval_LMC))
  # too-short sample fails loudly instead of returning NaN
  expect_error(spantest:::f_ttest(matrix(rnorm(6), 6, 1), k = 1/3), "subseries")
})
