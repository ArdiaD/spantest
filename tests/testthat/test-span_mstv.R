test_that("span_mstv returns the documented structure", {
  set.seed(123)
  R1 <- matrix(rnorm(250 * 3), 250, 3)
  R2 <- matrix(rnorm(250 * 100), 250, 100)

  result <- span_mstv(R1, R2)

  expect_type(result, "list")
  expect_named(result, c("pval", "stat", "H0", "crit", "Q",
                         "reject", "nu", "B"))
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_true(result$Q >= 0 && result$Q <= 1)
  expect_equal(result$Q * result$B, round(result$Q * result$B))   # a share of B draws
  expect_equal(result$nu, 5)                                        # the default
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

test_that("Q is the share of the B draws that do not reject, as in the authors' code", {
  # The derandomized rule draws B blocks of N normals after the one-shot draw, from
  # the same seed, and rejects when the share Q falls below (1 - tau) - B^(-1/4).
  set.seed(11)
  R1 <- matrix(rnorm(250 * 3), 250, 3)
  R2 <- matrix(rnorm(250 * 200), 250, 200)
  R2[, 1:5] <- R2[, 1:5] + 0.15                 # a few non-zero intercepts

  out <- span_mstv(R1, R2, control = list(seed = 7L))

  X    <- cbind(1, R1)
  bhat <- qr.solve(qr(X), R2)
  psi  <- sqrt(250) * (abs(bhat[1, ]) / sqrt(mean((R2 - X %*% bhat)^2)))^(5 / 2)
  B    <- floor(log(200)^2)
  set.seed(7)
  invisible(rnorm(200))                                   # the one-shot draw
  Zb   <- apply(psi + matrix(rnorm(200 * B), 200, B), 2, max)
  Qref <- mean(Zb <= out$crit)

  expect_equal(out$Q, Qref)
  expect_equal(out$reject, Qref < 0.95 - B^(-1 / 4))
  # another seed, other draws
  out2 <- span_mstv(R1, R2, control = list(seed = 8L, B = 1000L))
  expect_equal(out2$B, 1000L)
  expect_equal(out2$Q * 1000, round(out2$Q * 1000))
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
  expect_error(span_mstv(R1, R2, control = list(B = 0)), "B")
  expect_error(span_mstv(R1, R2, control = list(B = 2.5)), "B")
  expect_error(span_mstv(R1[-1, ], R2), "same number of rows")
})
