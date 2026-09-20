test_that("span_simulate returns correct dimensions", {
  set.seed(1)
  sim <- span_simulate(n = 200, K = 3, N = 8, ncp = 0)
  expect_type(sim, "list")
  expect_equal(dim(sim$R1), c(200L, 3L))
  expect_equal(dim(sim$R2), c(200L, 8L))
  expect_true(all(is.finite(sim$R1)))
  expect_true(all(is.finite(sim$R2)))
})

test_that("span_simulate is reproducible given the RNG seed", {
  set.seed(42); a <- span_simulate(120, 2, 5, ncp = 0.1, dgp = 4)
  set.seed(42); b <- span_simulate(120, 2, 5, ncp = 0.1, dgp = 4)
  expect_identical(a, b)
})

test_that("all DGP presets 1:12 run and return finite returns", {
  for (d in 1:12) {
    set.seed(d)
    sim <- span_simulate(150, 3, 6, ncp = 0, dgp = d)
    expect_equal(dim(sim$R2), c(150L, 6L), info = paste("dgp", d))
    expect_true(all(is.finite(sim$R2)), info = paste("dgp", d))
  }
})

test_that("the twelve presets are one innovation law per family, four dynamics", {
  # The grid is a 3 x 4 factorial: moving down a column must change the DYNAMICS
  # and nothing else. It did not always -- the two t presets without GARCH were
  # raw t_5 (variance 5/3) while the two with GARCH were standardised t_4, so the
  # rows of the manuscript's size tables were not comparable within the Student
  # family. This test is what stops that returning.
  tb <- span_dgp_table()
  for (fam in unique(tb$innovation)) {
    r <- tb[tb$innovation == fam, ]
    expect_equal(length(unique(r$df)), 1L, info = paste(fam, "df"))
    expect_equal(length(unique(r$standardize)), 1L, info = paste(fam, "standardize"))
    expect_setequal(r$dynamics, c("iid", "garch", "ar-garch", "ar"))
  }
  # and every heavy-tailed preset really is unit-variance, not merely flagged so
  for (d in tb$preset[tb$innovation != "normal"]) {
    set.seed(1)
    v <- var(as.numeric(span_simulate(1e5, 1, 1, dgp = d)$R1[, 1]))
    expect_lt(abs(v - 1), 0.15, label = paste0("var(preset ", d, ") = ", round(v, 3)))
  }
})

test_that("every preset carries unit PROCESS variance, not just unit innovations", {
  # Standardising the innovations is not enough: an AR(1) with unit innovations
  # has variance 1/(1 - phi^2) = 1.042, so the AR rows used to carry 4% more
  # noise than the iid ones and the twelve DGPs differed in level as well as in
  # shape. The GARCH parameters had always been chosen so that
  # omega/(1 - alpha - beta) = 1; scale = "process" completes that intent.
  for (d in seq_len(12)) {
    set.seed(7)
    v <- var(as.numeric(span_simulate(1e5, 1, 1, dgp = d)$R1[, 1]))
    expect_lt(abs(v - 1), 0.05, label = paste0("var(preset ", d, ") = ", round(v, 3)))
  }
})

test_that("scale = 'process' rescales without touching the dependence or the null", {
  # phi keeps its meaning: rescaling leaves the autocorrelation function alone.
  set.seed(3); a <- span_simulate(5e4, 1, 1, dynamics = "ar")                    # innovation
  set.seed(3); b <- span_simulate(5e4, 1, 1, dynamics = "ar", scale = "process")
  expect_lt(abs(stats::acf(a$R1[, 1], 1, plot = FALSE)$acf[2] -
                stats::acf(b$R1[, 1], 1, plot = FALSE)$acf[2]), 0.01)
  expect_lt(abs(var(a$R1[, 1]) - 1 / (1 - 0.2^2)), 0.05)   # textbook default
  expect_lt(abs(var(b$R1[, 1]) - 1), 0.05)

  # and under the spanning null it changes no test result, because it rescales
  # the whole system and the tests are scale-invariant
  set.seed(5); x <- span_simulate(250, 2, 8, ncp = 0, dynamics = "ar-garch",
                                  innovation = "skew-t", df = 4, scale = "innovation")
  set.seed(5); y <- span_simulate(250, 2, 8, ncp = 0, dynamics = "ar-garch",
                                  innovation = "skew-t", df = 4, scale = "process")
  expect_equal(span_hk(x$R1, x$R2)$pval, span_hk(y$R1, y$R2)$pval, tolerance = 1e-8)
  expect_equal(unlist(span_as(x$R1, x$R2, control = list(ks = 1/3, L = 2L))),
               unlist(span_as(y$R1, y$R2, control = list(ks = 1/3, L = 2L))),
               tolerance = 1e-8)
})

test_that("dgp must be a single integer in 1:12", {
  expect_error(span_simulate(50, 2, 3, dgp = 13))
  expect_error(span_simulate(50, 2, 3, dgp = 0))
})

test_that("ncp = 0 gives near-nominal GRS size (spanning null holds)", {
  skip_on_cran()
  p <- vapply(seq_len(300), function(s) {
    set.seed(s)
    sim <- span_simulate(250, 3, 8, ncp = 0, dgp = 1)
    span_grs(sim$R1, sim$R2)$pval
  }, numeric(1))
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.03)
})

test_that("skew-t innovation path works", {
  set.seed(1)
  sim <- span_simulate(120, 2, 4, innovation = "skew-t", dynamics = "iid",
                       df = 4, xi = 0.9)
  expect_equal(dim(sim$R2), c(120L, 4L))
  expect_true(all(is.finite(sim$R2)))
})

test_that("gamma = NULL reproduces the simulation as it was before the argument existed", {
  # The pervasive-factor hook draws one extra series, placed AFTER the K + N
  # draws that build Z and E. When gamma is NULL it draws nothing at all, so
  # every cell of the published size/power study must come back unchanged. The
  # reference here is the generative core transcribed without the hook.
  f_ref <- function(n, K, N, ncp, dgp, rho_factor, rho_error, burnin = 500L) {
    p <- spantest:::.DGP_PRESETS[match(toupper(dgp), toupper(spantest:::.DGP_PRESETS$name)), ]
    om <- 0.1; al <- 0.1; be <- 0.8; ar <- 0.2
    v <- if (p$dynamics %in% c("garch", "ar-garch")) om / (1 - al - be) else 1
    if (p$dynamics %in% c("ar", "ar-garch")) v <- v / (1 - ar^2)
    sd_proc <- sqrt(v)
    one_series <- function() {
      m <- n + burnin
      z <- switch(p$innovation,
        "normal" = rnorm(m),
        "t"      = { r <- rt(m, p$df); if (p$standardize) r / sqrt(p$df / (p$df - 2)) else r },
        "skew-t" = spantest:::f_rsstd(m, p$df, p$xi))
      s <- switch(p$dynamics,
        "iid"      = z,
        "garch"    = spantest:::garch_filter(z, om, al, be),
        "ar"       = as.numeric(filter(z, ar, method = "recursive")),
        "ar-garch" = as.numeric(filter(spantest:::garch_filter(z, om, al, be), ar, method = "recursive")))
      s[(burnin + 1L):m] / sd_proc
    }
    Z <- matrix(0, n, K); for (j in seq_len(K)) Z[, j] <- one_series()
    E <- matrix(0, n, N); for (j in seq_len(N)) E[, j] <- one_series()
    z   <- Z %*% chol(toeplitz(rho_factor^(seq_len(K) - 1L)))
    eps <- E %*% chol(toeplitz(rho_error^(seq_len(N) - 1L)))
    y <- sweep(z %*% matrix(1 + ncp, K, N) + eps, 2, rep(ncp, N), "+")
    list(R1 = cbind(z[, 1], z[, -1, drop = FALSE] + z[, 1]), R2 = y)
  }

  for (d in c("iid-N", "GARCH-ST", "AR-SKST", "AR-GARCH-N")) {
    for (cfg in list(c(n = 250, K = 2, N = 2), c(n = 250, K = 5, N = 100),
                     c(n = 120, K = 100, N = 300))) {
      set.seed(11)
      got <- span_simulate(n = cfg[["n"]], K = cfg[["K"]], N = cfg[["N"]], ncp = 0.2,
                           dgp = d, rho_factor = 0.8, rho_error = 0.5)
      set.seed(11)
      ref <- f_ref(cfg[["n"]], cfg[["K"]], cfg[["N"]], 0.2, d, 0.8, 0.5)
      expect_identical(got$R1, ref$R1, info = paste(d, cfg[["K"]], cfg[["N"]]))
      expect_identical(got$R2, ref$R2, info = paste(d, cfg[["K"]], cfg[["N"]]))
    }
  }
})

test_that("gamma adds a common component that does not die out across the cross-section", {
  n <- 2000L; N <- 50L
  set.seed(3)
  a <- span_simulate(n = n, K = 2L, N = N, ncp = 0, dgp = "iid-N",
                     rho_factor = 0.8, rho_error = 0)
  set.seed(3)
  b <- span_simulate(n = n, K = 2L, N = N, ncp = 0, dgp = "iid-N",
                     rho_factor = 0.8, rho_error = 0, gamma = 1)

  # gamma acts on the idiosyncratic terms, so look at them: the test assets all
  # load on the same factors and are strongly correlated whatever gamma does.
  # With rho_error = 0 the residuals are independent across assets; gamma = 1
  # adds the same unit-variance series to each, taking every pairwise residual
  # correlation from (almost) zero to (almost) one half -- including between the
  # first asset and the last, which is what "does not die out" means.
  f_resid_cor <- function(s) {
    X <- cbind(1, s$R1)
    cor(s$R2 - X %*% qr.solve(qr(X), s$R2))
  }
  cor_a <- f_resid_cor(a); cor_b <- f_resid_cor(b)
  off <- upper.tri(cor_a)
  expect_lt(mean(abs(cor_a[off])), 0.05)
  expect_gt(min(cor_b[off]), 0.3)
  expect_equal(mean(cor_b[off]), 0.5, tolerance = 0.05)
  expect_equal(cor_b[1, N], cor_b[1, 2], tolerance = 0.1)

  # A vector of loadings is accepted, a wrong length is not.
  set.seed(3)
  expect_silent(span_simulate(n = 100L, K = 2L, N = N, dgp = "iid-N", gamma = runif(N, 0.7, 0.9)))
  expect_error(span_simulate(n = 100L, K = 2L, N = N, dgp = "iid-N", gamma = c(0.1, 0.2)),
               "length N")
})
