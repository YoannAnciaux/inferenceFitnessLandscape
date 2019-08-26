test_that("'ptof_fgm_iso' returns a vector of numeric fitness", {
  checkmate::expect_numeric(ptof_fgm_iso(matrix(1:9, 3, 3), 1), finite = T,
                          any.missing = F, len = 3, null.ok = F)
})
test_that("'ptof_fgm_iso' of positive matrix returns the expected fitness", {
  expect_equal(ptof_fgm_iso(matrix(1:9, 3, 3), 1), c(-32.0, -45.5, -62.0))
})
test_that("'ptof_fgm_iso' of single row matrix returns the expected fitness", {
  expect_equal(ptof_fgm_iso(matrix(1:9, 3, 3)[1, , drop = FALSE], 1), -32)
})
test_that("'ptof_fgm_iso' of single column matrix returns the expected fitness", {
  expect_equal(ptof_fgm_iso(matrix(1:9, 3, 3)[, 1, drop = FALSE], 1), c(0.5, -1.0, -3.5))
})
test_that("'ptof_fgm_iso' of single value matrix returns the expected fitness", {
  expect_equal(ptof_fgm_iso(matrix(1:9, 3, 3)[1, 1, drop = FALSE], 1), 0.5)
})
test_that("'ptof_fgm_iso' of negative matrix returns the expected fitness", {
  expect_equal(ptof_fgm_iso(-matrix(1:9, 3, 3)[1, 1, drop = FALSE], 1), 0.5)
})
test_that("'ptof_fgm_iso' with non-zero position for optimum returns the expected fitness", {
  expect_equal(ptof_fgm_iso(matrix(1:9, 3, 3), 1, pheno_opt = c(-1, 1, 1)), c(-23.5, -36.0, -51.5))
})


