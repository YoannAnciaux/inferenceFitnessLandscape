test_that("'simulate_fl' of a simple function with multiple parameters and multiple simulations returns the expected matrix",{
  simulation_model <- function(m, sd) {rnorm(n = 12, mean = m, sd = sd)}
  parameter <- cbind(1:10, seq(0.1, 1, 0.1))
  s <- simulate_fl(parameter, simulation_model, ncore = 1)
  sp <- simulate_fl(parameter, simulation_model, ncore = 2)
  expect_equal(unname(s$parameter), parameter)
  expect_equal(unname(sp$parameter), parameter)
  checkmate::expect_matrix(s$simulation, mode = "numeric", any.missing = F, nrows = 10, ncols = 12)
  checkmate::expect_matrix(sp$simulation, mode = "numeric", any.missing = F, nrows = 10, ncols = 12)
})
test_that("'simulate_fl' of a simple function with a single parameter and multiple simulations returns the expected matrix",{
  simulation_model <- function(m) {rnorm(n = 12, mean = m, sd = 0.1)}
  parameter <- matrix(1:10, 10, 1)
  s <- simulate_fl(parameter, simulation_model, ncore = 1)
  sp <- simulate_fl(parameter, simulation_model, ncore = 2)
  expect_equal(unname(s$parameter), parameter)
  expect_equal(unname(sp$parameter), parameter)
  checkmate::expect_matrix(s$simulation, mode = "numeric", any.missing = F, nrows = 10, ncols = 12)
  checkmate::expect_matrix(sp$simulation, mode = "numeric", any.missing = F, nrows = 10, ncols = 12)

})
test_that("'simulate_fl' of a simple function with multiple parameters and a single simulation returns the expected matrix",{
  simulation_model <- function(m, sd) {rnorm(n = 12, mean = m, sd = sd)}
  parameter <- matrix(c(1, 0.1), 1, 2)
  s <- simulate_fl(parameter, simulation_model, ncore = 1)
  sp <- simulate_fl(parameter, simulation_model, ncore = 2)
  expect_equal(unname(s$parameter), parameter)
  expect_equal(unname(sp$parameter), parameter)
  checkmate::expect_matrix(s$simulation, mode = "numeric", any.missing = F, nrows = 1, ncols = 12)
  checkmate::expect_matrix(sp$simulation, mode = "numeric", any.missing = F, nrows = 1, ncols = 12)
})
test_that("'simulate_fl' of a simple function with a single parameter and a single simulation returns the expected matrix",{
  simulation_model <- function(m) {rnorm(n = 12, mean = m, sd = 0.1)}
  parameter <- matrix(1, 1, 1)
  s <- simulate_fl(parameter, simulation_model, ncore = 1)
  sp <- simulate_fl(parameter, simulation_model, ncore = 2)
  expect_equal(unname(s$parameter), parameter)
  expect_equal(unname(sp$parameter), parameter)
  checkmate::expect_matrix(s$simulation, mode = "numeric", any.missing = F, nrows = 1, ncols = 12)
  checkmate::expect_matrix(sp$simulation, mode = "numeric", any.missing = F, nrows = 1, ncols = 12)
})
