# test_that("span_as returns aggregated p-values", {
#   set.seed(1)
#   bench <- matrix(rnorm(300), nrow = 100, ncol = 3)
#   test  <- matrix(rnorm(200), nrow = 100, ncol = 2)
#
#   result <- span_as(bench, test, control = list(ks = c(1/3), L = c(0)))
#
#   expect_type(result, "list")
#   expect_true(all(grepl("^CCT[ad]{1,2}_L\\d+_k\\d+$", names(result))))
#   expect_true(all(sapply(result, function(x) is.numeric(x) && x >= 0 && x <= 1)))
# })
