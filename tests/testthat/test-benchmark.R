test_that("benchmark_gs_boundaries returns one tidy row per analysis", {
  b <- benchmark_gs_boundaries(k = 3, type = "OF")

  expect_s3_class(b, "data.frame")
  expect_identical(nrow(b), 3L)
  expect_identical(
    names(b),
    c(
      "type", "k", "analysis", "information_rate", "rpact_z", "gsdesign_z",
      "abs_diff", "tolerance", "agrees"
    )
  )
  expect_identical(b$analysis, 1:3)
  expect_type(b$agrees, "logical")
})

test_that("rpact and gsDesign agree on classical O'Brien-Fleming bounds", {
  for (k in 2:5) {
    b <- benchmark_gs_boundaries(k = k, type = "OF", tolerance = 1e-4)
    expect_true(all(b$agrees), info = paste("OF k =", k))
    expect_lt(max(b$abs_diff), 1e-5)
  }
})

test_that("rpact and gsDesign agree on classical Pocock bounds", {
  for (k in 2:5) {
    b <- benchmark_gs_boundaries(k = k, type = "P", tolerance = 1e-4)
    expect_true(all(b$agrees), info = paste("Pocock k =", k))
    expect_lt(max(b$abs_diff), 1e-5)
  }
})

test_that("the spending-function variants also agree within tolerance", {
  for (type in c("asOF", "asP")) {
    b <- benchmark_gs_boundaries(k = 4, type = type, tolerance = 1e-4)
    expect_true(all(b$agrees), info = type)
  }
})

test_that("the benchmark reproduces the boundaries the designs actually use", {
  b <- benchmark_gs_boundaries(k = 3, alpha = 0.025, beta = 0.1, type = "asOF")
  g <- design_group_sequential(
    scenario("ercp_acute_cholangitis"),
    k = 3, type_of_design = "asOF"
  )

  expect_equal(b$rpact_z, g$detail$efficacy_z)
})

test_that("classical bounds agree more tightly than spending approximations", {
  classical <- max(benchmark_gs_boundaries(k = 5, type = "OF")$abs_diff)
  spending <- max(benchmark_gs_boundaries(k = 5, type = "asOF")$abs_diff)

  expect_lt(classical, spending)
})

test_that("an impossibly tight tolerance is reported as disagreement", {
  b <- benchmark_gs_boundaries(k = 3, type = "asOF", tolerance = 1e-12)

  expect_false(all(b$agrees))
  expect_identical(unique(b$tolerance), 1e-12)
})

test_that("unequal information rates are passed to both engines", {
  ir <- c(0.4, 0.7, 1)
  b <- benchmark_gs_boundaries(k = 3, type = "OF", information_rates = ir)

  expect_equal(b$information_rate, ir)
  expect_true(all(b$agrees))
  equal_spaced <- benchmark_gs_boundaries(k = 3, type = "OF")
  expect_false(isTRUE(all.equal(b$rpact_z, equal_spaced$rpact_z)))
})

test_that("benchmark_gs_boundaries validates its arguments", {
  expect_error(benchmark_gs_boundaries(k = 1), "`k`")
  expect_error(benchmark_gs_boundaries(k = 3.5), "`k`")
  expect_error(benchmark_gs_boundaries(k = 3, alpha = 0), "alpha")
  expect_error(benchmark_gs_boundaries(k = 3, beta = 1), "beta")
  expect_error(benchmark_gs_boundaries(k = 3, tolerance = 0), "tolerance")
  expect_error(benchmark_gs_boundaries(k = 3, type = "WT"), "type_of_design")
})

test_that("benchmark_report summarises one row per design and all agree", {
  r <- benchmark_report()

  expect_s3_class(r, "data.frame")
  expect_identical(
    names(r),
    c("type", "k", "n_analyses", "max_abs_diff", "tolerance", "all_agree")
  )
  expect_identical(nrow(r), 16L)
  expect_true(all(r$all_agree))
  expect_true(all(r$n_analyses == r$k))
  expect_true(all(r$max_abs_diff < 1e-4))
})

test_that("benchmark_report honours a restricted grid", {
  r <- benchmark_report(k = 2:3, types = c("OF", "P"))

  expect_identical(nrow(r), 4L)
  expect_setequal(r$type, c("OF", "P"))
  expect_setequal(r$k, 2:3)
  expect_true(all(r$all_agree))
})

test_that("benchmark_report validates its arguments", {
  expect_error(benchmark_report(k = 1), "`k`")
  expect_error(benchmark_report(types = character()), "types")
  expect_error(benchmark_report(types = "nonsense"), "type_of_design")
})
