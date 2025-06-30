test_that("f_hk returns full suite of test outputs", {
  set.seed(106)
  R1 <- matrix(rnorm(300), nrow = 100, ncol = 3)
  R2 <- matrix(rnorm(200), nrow = 100, ncol = 2)

  result <- f_hk(R1, R2)
  keys <- c("HK", "F1", "F2", "CCTa", "CCTd", "CCTad")
  expect_type(result, "list")
  expect_true(all(keys %in% names(result)))

  for (key in keys) {
    expect_type(result[[key]], "list")
    expect_true("pval" %in% names(result[[key]]))
    expect_true(result[[key]]$pval >= 0 && result[[key]]$pval <= 1)
  }
})
