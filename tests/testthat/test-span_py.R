test_that("span_py returns expected structure", {
  set.seed(123)
  R1 <- matrix(rnorm(300), 100, 3)
  R2 <- matrix(rnorm(200), 100, 2)

  result <- span_py(R1, R2)

  expect_type(result, "list")
  expect_named(result, c("pval", "stat", "H0"))
  expect_type(result$pval, "double")
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_type(result$stat, "double")
  expect_type(result$H0, "character")
  expect_equal(result$H0, "alpha = 0")
})

test_that("span_py handles different dimensions correctly", {
  set.seed(42)
  R1 <- matrix(rnorm(400), 100, 4)
  R2 <- matrix(rnorm(300), 100, 3)

  result <- span_py(R1, R2)

  expect_true(is.numeric(result$stat))
  expect_true(is.numeric(result$pval))
  expect_true(length(result$stat) == 1)
  expect_true(length(result$pval) == 1)
})

test_that("span_py returns NA for a single test asset (N = 1)", {
  set.seed(1)
  R1 <- matrix(rnorm(300), 100, 3)
  R2 <- matrix(rnorm(100), 100, 1)

  result <- span_py(R1, R2)

  expect_true(is.na(result$pval))
  expect_true(is.na(result$stat))
})


test_that("span_py is defined when N exceeds T, and its correlation sum is the pairwise one", {
  set.seed(3)
  T <- 120L; K <- 3L; N <- 300L
  R1 <- matrix(rnorm(T * K), T, K)
  R2 <- R1 %*% matrix(runif(K * N), K, N) + matrix(rnorm(T * N), T, N)
  result <- span_py(R1, R2)
  expect_true(is.finite(result$stat))
  expect_true(result$pval >= 0 && result$pval <= 1)

  # the vectorised average of thresholded squared correlations equals the double loop
  E <- R2 - cbind(1, R1) %*% solve(crossprod(cbind(1, R1)), crossprod(cbind(1, R1), R2))
  S <- crossprod(E) / T; v <- T - K - 1; th <- qnorm(1 - 0.05 / (N - 1) / 2)^2
  loop <- 0
  for (i in 2:N) for (j in 1:(i - 1)) {
    r2 <- S[i, j]^2 / (S[i, i] * S[j, j])
    if (v * r2 >= th) loop <- loop + r2
  }
  rho2 <- (S / sqrt(outer(diag(S), diag(S))))[upper.tri(S)]^2
  expect_equal(sum(rho2[v * rho2 >= th]), loop, tolerance = 1e-12)
})
