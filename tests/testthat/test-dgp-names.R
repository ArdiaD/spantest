test_that("every name maps to its preset and produces identical data", {
  tab <- span_dgp_table()
  expect_equal(nrow(tab), 12L)
  expect_equal(tab$name, DGP_NAMES)

  for (i in seq_len(12L)) {
    set.seed(11); by_num  <- span_simulate(80, 2, 3, dgp = i)
    set.seed(11); by_name <- span_simulate(80, 2, 3, dgp = DGP_NAMES[i])
    expect_identical(by_num, by_name,
                     info = paste0("preset ", i, " = ", DGP_NAMES[i]))
  }
})

test_that("names are case-insensitive and tolerate surrounding space", {
  set.seed(3); ref <- span_simulate(60, 2, 2, dgp = 9)
  for (v in c("AR-GARCH-SKST", "ar-garch-skst", "Ar-Garch-SkSt", "  AR-GARCH-SKST  ")) {
    set.seed(3); expect_identical(span_simulate(60, 2, 2, dgp = v), ref, info = v)
  }
})

test_that("the table's stated process matches what the preset actually sets", {
  # A name must not drift from the process it claims to denote. Rebuild each
  # preset from its tabulated innovation/dynamics/df/standardize and require the
  # simulated data to match.
  tab <- span_dgp_table()
  for (i in seq_len(12L)) {
    set.seed(5)
    a <- span_simulate(80, 2, 3, dgp = tab$name[i])
    set.seed(5)
    b <- span_simulate(80, 2, 3, innovation = tab$innovation[i],
                       dynamics = tab$dynamics[i], df = tab$df[i],
                       standardize = tab$standardize[i], xi = 0.9)
    expect_identical(a, b, info = paste0(tab$name[i], " (preset ", i, ")"))
  }
})

test_that("the AR and AR-GARCH blocks are not transposed", {
  # The specific mix-up this naming exists to prevent.
  tab <- span_dgp_table()
  expect_equal(tab$dynamics[tab$name == "AR-GARCH-SKST"], "ar-garch")
  expect_equal(tab$dynamics[tab$name == "AR-SKST"],       "ar")
  expect_equal(tab$preset[tab$name == "AR-GARCH-SKST"],   9L)
  expect_equal(tab$preset[tab$name == "AR-SKST"],        12L)
})

test_that("bad dgp values are rejected with a usable message", {
  expect_error(span_simulate(60, 2, 2, dgp = "nope"), "unknown")
  expect_error(span_simulate(60, 2, 2, dgp = 13))
  expect_error(span_simulate(60, 2, 2, dgp = 0))
  expect_error(span_simulate(60, 2, 2, dgp = c("iid-N", "GARCH-N")))
})
