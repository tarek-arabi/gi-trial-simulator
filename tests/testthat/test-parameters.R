temp_pack_dir <- function() {
  dir <- file.path(tempdir(), paste0("packs-", as.integer(stats::runif(1, 1, 1e9))))
  dir.create(dir, recursive = TRUE)
  dir
}

with_pack_options <- function(opts, code) {
  old <- options(opts)
  on.exit(options(old), add = TRUE)
  force(code)
}

test_that("shipped packs load and validate", {
  for (id in c("ercp_acute_cholangitis", "hrs_terlipressin")) {
    pack <- load_pack(id)
    expect_s3_class(pack, "gi_pack")
    expect_silent(validate_pack(pack))
    expect_identical(pack$id, id)
  }
})

test_that("list_packs reports the shipped packs", {
  packs <- list_packs()
  expect_true(all(c("ercp_acute_cholangitis", "hrs_terlipressin") %in% packs$id))
  expect_false(any(packs$restricted))
  expect_true(all(packs$provenance == "published_literature"))
})

test_that("an unknown pack id fails with a message naming what is available", {
  expect_error(load_pack("no_such_pack"), "ercp_acute_cholangitis")
})

test_that("the flagship scenario carries the published anchoring-trial rates", {
  s <- scenario("ercp_acute_cholangitis")
  expect_s3_class(s, "gi_scenario")
  expect_identical(s$endpoint, "mortality_30d")
  expect_equal(s$control_rate, 0.0658)
  expect_equal(s$treatment_rate, 0.0395)
  expect_identical(s$direction, "lower_is_better")
  expect_false(s$overridden)
})

test_that("rates can be overridden for sensitivity analysis and the override is flagged", {
  s <- scenario("ercp_acute_cholangitis", treatment_rate = 0.05)
  expect_equal(s$treatment_rate, 0.05)
  expect_equal(s$control_rate, 0.0658)
  expect_true(s$overridden)
})

test_that("a non-primary endpoint can be selected by name", {
  s <- scenario("ercp_acute_cholangitis", endpoint = "post_ercp_adverse_events")
  expect_equal(s$control_rate, 0.092)
  expect_equal(s$treatment_rate, 0.171)
})

test_that("an unknown endpoint fails with a message listing the real ones", {
  expect_error(
    scenario("ercp_acute_cholangitis", endpoint = "nonexistent"),
    "mortality_30d"
  )
})

test_that("out-of-range rates are rejected", {
  expect_error(scenario("ercp_acute_cholangitis", control_rate = 0), "between 0 and 1")
  expect_error(scenario("ercp_acute_cholangitis", control_rate = 1), "between 0 and 1")
  expect_error(scenario("ercp_acute_cholangitis", treatment_rate = -0.1), "between 0 and 1")
  expect_error(scenario("ercp_acute_cholangitis", control_rate = NA_real_), "between 0 and 1")
})

test_that("validate_pack rejects an out-of-range event rate", {
  pack <- load_pack("ercp_acute_cholangitis")
  pack$endpoints$mortality_30d$control_rate <- 1.5
  expect_error(validate_pack(pack), "strictly between 0 and 1")
})

test_that("validate_pack rejects an endpoint citing an undeclared source", {
  pack <- load_pack("ercp_acute_cholangitis")
  pack$endpoints$mortality_30d$source <- "not_a_declared_source"
  expect_error(validate_pack(pack), "not declared in sources")
})

test_that("validate_pack rejects a missing required field", {
  pack <- load_pack("hrs_terlipressin")
  pack$sources <- NULL
  expect_error(validate_pack(pack), "missing required field 'sources'")
})

test_that("validate_pack rejects an unrecognised provenance", {
  pack <- load_pack("hrs_terlipressin")
  pack$provenance <- "vibes"
  expect_error(validate_pack(pack), "provenance")
})

test_that("every endpoint in a shipped pack declares a verification date", {
  for (id in c("ercp_acute_cholangitis", "hrs_terlipressin")) {
    pack <- load_pack(id)
    for (key in names(pack$endpoints)) {
      expect_false(
        is.null(pack$endpoints[[key]]$verified),
        info = paste0(id, " / ", key, " has no verified date")
      )
    }
  }
})

test_that("packs outside the package take precedence over shipped packs", {
  dir <- temp_pack_dir()
  pack <- yaml::read_yaml(system.file(
    "parameters", "hrs_terlipressin.yaml",
    package = "gitrialsim"
  ))
  pack$title <- "Overridden by an external pack"
  pack$provenance <- "proprietary"
  pack$restricted <- TRUE
  yaml::write_yaml(pack, file.path(dir, "hrs_terlipressin.yaml"))

  with_pack_options(list(gitrialsim.pack_paths = dir), {
    expect_true(dir %in% pack_search_path())
    loaded <- load_pack("hrs_terlipressin")
    expect_identical(loaded$title, "Overridden by an external pack")
    expect_true(isTRUE(loaded$restricted))
    listed <- list_packs()
    expect_equal(sum(listed$id == "hrs_terlipressin"), 1L)
  })

  expect_identical(load_pack("hrs_terlipressin")$provenance, "published_literature")
})

test_that("a pack can be loaded directly from a file path", {
  path <- system.file("parameters", "ercp_acute_cholangitis.yaml", package = "gitrialsim")
  expect_identical(load_pack(path)$id, "ercp_acute_cholangitis")
})
