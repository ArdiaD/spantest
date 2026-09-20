test_that("span_mstv returns the documented structure", {
  set.seed(123)
  R1 <- matrix(rnorm(250 * 3), 250, 3)
  R2 <- matrix(rnorm(250 * 100), 250, 100)

  result <- span_mstv(R1, R2)

  expect_type(result, "list")
  expect_named(result, c("pval", "stat", "H0", "crit", "Q", "logQ",
                         "reject", "nu", "B"))
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_true(result$Q >= 0 && result$Q <= 1)
  expect_equal(result$Q, exp(result$logQ))
  expect_type(result$reject, "logical")
  expect_equal(result$H0, "alpha = 0")
  expect_equal(result$B, floor(log(100)^2))   # their guideline, B = floor(log N)^2
})

test_that("the statistic is the one in the paper", {
  # An independent transcription of equations (3.1)-(3.4) and the extreme-value
  # constants, written out rather than factored, so that a change in span_mstv()
  # has to be matched here deliberately.
  set.seed(4)
  Tn <- 300L; K <- 4L; N <- 60L
  R1 <- matrix(rnorm(Tn * K), Tn, K)
  R2 <- matrix(rnorm(Tn * N), Tn, N)
  nu <- 4; tau <- 0.05; seed <- 123L

  X    <- cbind(1, R1)
  bhat <- solve(crossprod(X), crossprod(X, R2))
  res  <- R2 - X %*% bhat
  sNT  <- sqrt(sum(res^2) / (N * Tn))
  psi  <- (Tn^(1 / nu) * abs(bhat[1, ]) / sNT)^(nu / 2)   # the (T^{1/nu} .)^{nu/2} form
  bN   <- sqrt(2 * log(N)) - (log(log(N)) + log(4 * pi)) / (2 * sqrt(2 * log(N)))
  aN   <- bN / (1 + bN^2)
  cr   <- bN - aN * log(-log(1 - tau))
  set.seed(seed)
  Zref <- max(psi + rnorm(N))

  out <- span_mstv(R1, R2, control = list(nu = nu, tau = tau, seed = seed))

  expect_equal(out$stat, Zref)
  expect_equal(out$crit, cr)
  expect_equal(out$pval, 1 - exp(-exp(-(Zref - bN) / aN)))
})

test_that("the closed-form Q is the limit of their B replications", {
  # Q is the share of perturbation draws that do not reject. span_mstv() computes
  # it exactly, as prod_i Phi(c - psi_i); this is the average they simulate.
  set.seed(11)
  R1 <- matrix(rnorm(250 * 3), 250, 3)
  R2 <- matrix(rnorm(250 * 200), 250, 200)
  R2[, 1:5] <- R2[, 1:5] + 0.15                 # a few non-zero intercepts

  out <- span_mstv(R1, R2)

  X    <- cbind(1, R1)
  bhat <- qr.solve(qr(X), R2)
  psi  <- sqrt(250) * (abs(bhat[1, ]) / sqrt(mean((R2 - X %*% bhat)^2)))^(4 / 2)
  set.seed(2)
  Qmc <- mean(replicate(20000, max(psi + rnorm(200)) <= out$crit))

  expect_equal(out$Q, Qmc, tolerance = 0.01)
})

test_that("span_mstv leaves the caller's RNG stream untouched", {
  # The simulation study runs sign-flip tests after this one; a shifted stream
  # would silently change their results.
  set.seed(5)
  R1 <- matrix(rnorm(250 * 3), 250, 3)
  R2 <- matrix(rnorm(250 * 50), 250, 50)

  set.seed(99); before <- rnorm(5)
  set.seed(99); invisible(span_mstv(R1, R2)); after <- rnorm(5)

  expect_equal(before, after)
})

test_that("span_mstv returns NA where its reference distribution does not apply", {
  set.seed(6)
  R1 <- matrix(rnorm(250 * 3), 250, 3)
  R2 <- matrix(rnorm(250 * 50), 250, 50)

  # Below N = 10 the Gumbel limit for the maximum of N normals is worthless.
  small <- span_mstv(R1, R2[, 1:9, drop = FALSE])
  expect_true(is.na(small$pval))
  expect_true(is.na(small$stat))
  expect_equal(small$H0, "alpha = 0")

  # Fewer observations than parameters.
  short <- span_mstv(R1[1:3, , drop = FALSE], R2[1:3, , drop = FALSE])
  expect_true(is.na(short$pval))
})

test_that("span_mstv rejects when the intercepts are not zero", {
  set.seed(8)
  Tn <- 250L; N <- 200L
  R1 <- matrix(rnorm(Tn * 3), Tn, 3)
  R2 <- matrix(rnorm(Tn * N), Tn, N)
  R2[, 1:10] <- R2[, 1:10] + 1                  # a sparse, large alternative

  out <- span_mstv(R1, R2)

  expect_lt(out$pval, 0.01)
  expect_true(out$reject)
})

test_that("nu and tau are validated", {
  set.seed(9)
  R1 <- matrix(rnorm(250 * 3), 250, 3)
  R2 <- matrix(rnorm(250 * 50), 250, 50)

  expect_error(span_mstv(R1, R2, control = list(nu = 3)), "nu")
  expect_error(span_mstv(R1, R2, control = list(tau = 1)), "tau")
  expect_error(span_mstv(R1[-1, ], R2), "same number of rows")
})
