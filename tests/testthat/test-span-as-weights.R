# The multiplier weights of span_as() at L > 0: one draw per statistic (version 1.4-3).

test_that("f_mult draws an independent T x N matrix for the alpha and for the delta scores", {
  W <- spantest:::f_mult(50L, 4L, 2L, cseed = 1L)
  expect_identical(dim(W$A), c(50L, 4L))
  expect_identical(dim(W$D), c(50L, 4L))
  expect_false(isTRUE(all.equal(W$A, W$D)))                  # alpha and delta differ
  expect_false(isTRUE(all.equal(W$A[, 1], W$A[, 2])))        # assets differ
  expect_identical(spantest:::f_mult(50L, 4L, 2L, cseed = 1L), W)   # the seed fixes the draw
  expect_identical(spantest:::f_mult(50L, 4L, 0L), list(A = 1, D = 1))
})

test_that("f_mult leaves the caller's random-number stream untouched", {
  set.seed(42); a <- runif(3)
  set.seed(42); invisible(spantest:::f_mult(30L, 5L, 2L, cseed = 9L)); b <- runif(3)
  expect_identical(a, b)
})

test_that("the weights have the moments of a product of L independent N(1, 1)", {
  W <- spantest:::f_mult(4000L, 50L, 2L, cseed = 3L)
  w <- c(W$A, W$D)
  expect_equal(mean(w), 1, tolerance = 0.02)                 # E[kappa] = 1
  expect_equal(mean(w^2), 4, tolerance = 0.05)               # E[kappa^2] = 2^L
  expect_lt(abs(stats::cor(W$A[, 1], W$A[, 2])), 0.06)       # independent across assets
  expect_lt(abs(stats::cor(W$A[, 1], W$D[, 1])), 0.06)       # and between the two scores
})

test_that("two identical assets share the L = 0 p-value but not the L = 2 one", {
  set.seed(5)
  x <- matrix(rnorm(250 * 2), 250, 2)
  y <- matrix(rnorm(250), 250, 1) + x[, 1]
  pv <- spantest:::f_getpv_batch(x, cbind(y, y), ks = 1/3, L = c(0, 2))
  for (h in c("CCTa", "CCTd", "CCTad")) {
    expect_identical(pv[[paste0(h, "_L0_k1")]][1], pv[[paste0(h, "_L0_k1")]][2], info = h)
    expect_false(isTRUE(all.equal(pv[[paste0(h, "_L2_k1")]][1], pv[[paste0(h, "_L2_k1")]][2])), info = h)
  }
})

test_that("the results at L = 0 are those of version 1.4-2", {
  set.seed(20261005)
  x <- matrix(rnorm(250 * 3), 250, 3)
  y <- matrix(rnorm(250 * 8), 250, 8) + x[, 1]
  r <- unlist(span_as(x, y, control = list(ks = c(1/3, 1/2), L = c(0, 2))))
  ref <- c(CCTd_L0_k1 = 0.44816218271603897, CCTd_L0_k2 = 0.23901130551202276,
           CCTad_L0_k1 = 0.514253950872117, CCTad_L0_k2 = 0.4220984486685691,
           CCTa_L0_k1 = 0.5791549786353043, CCTa_L0_k2 = 0.66540083497502844)
  expect_equal(r[names(ref)], ref, tolerance = 1e-12)
})

test_that("at L = 2 each test responds to the hypothesis it claims and ignores the other", {
  skip_on_cran()
  set.seed(20261005)
  T <- 250L; K <- 3L; N <- 10L; NREP <- 200L
  gen <- function(alpha, colsum) {
    R1 <- matrix(rnorm(T * K), T, K)
    B  <- matrix(runif(K * N), K, N); B <- sweep(B, 2, colSums(B), "/") * colsum
    list(R1 = R1, R2 = sweep(R1 %*% B, 2, alpha, "+") + matrix(rnorm(T * N, sd = 0.05), T, N))
  }
  rate <- function(alpha, colsum) {
    p <- vapply(seq_len(NREP), function(i) {
      d <- gen(alpha, colsum)
      r <- span_as(d$R1, d$R2, control = list(L = 2L, seed = i))
      c(a = r$CCTa_L2_k1, d = r$CCTd_L2_k1)
    }, numeric(2))
    rowMeans(p < 0.05)
  }
  ra <- rate(rep(0.03, N), 1.00)                             # alpha != 0, delta = 0
  rd <- rate(rep(0, N), 1.25)                                # delta != 0, alpha = 0
  expect_gt(ra[["a"]], 0.90); expect_lt(ra[["d"]], 0.15)
  expect_gt(rd[["d"]], 0.90); expect_lt(rd[["a"]], 0.15)
})
