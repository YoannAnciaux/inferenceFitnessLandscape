test_that("'ftop_fgm_iso' returns numeric coordinates", {
  checkmate::expect_numeric(ftop_fgm_iso(fitness = 0, n = 3, maxfitness = 1, alpha = 1/2, Q = 2),
               finite = T, any.missing = F, len = 3, null.ok = F)
})
test_that("'ftop_fgm_iso' of different positive fitness returns expected coordinates", {
  expect_equal(round(ftop_fgm_iso(fitness = 0, n = 3, maxfitness = 1, alpha = 1/2, Q = 2), digits = 6),
               matrix(c(-1.414214, 0.000000, 0.000000), 1, 3))
  expect_equal(ftop_fgm_iso(fitness = 1, n = 3, maxfitness = 1, alpha = 1/2, Q = 2),
               matrix(numeric(3), 1, 3))
  expect_equal(ftop_fgm_iso(fitness = 0.5, n = 3, maxfitness = 1, alpha = 1/2, Q = 2),
               matrix(c(-1, 0, 0), 1, 3))
})
test_that("'ftop_fgm_iso' of non-zero pheno_opt returns expected coordinates", {
  expect_equal(round(ftop_fgm_iso(fitness = 0, n = 3, maxfitness = 1, alpha = 1/2, Q = 2, pheno_opt = c(1, 1, 1)), digits = 7),
               matrix(c(-0.4142136, 1.000000, 1.0000000), 1, 3))
  expect_equal(ftop_fgm_iso(fitness = 1, n = 3, maxfitness = 1, alpha = 1/2, Q = 2, pheno_opt = c(-1, 1, 1)),
               matrix(c(-1, 1, 1), 1, 3))
})
test_that("'ftop_fgm_iso' of negative fitness returns expected coordinates", {
  expect_equal(ftop_fgm_iso(fitness = -1.5, n = 3, maxfitness = -1, alpha = 1/2, Q = 2, pheno_opt = c(0.1, 1, 1)),
               matrix(c(-0.9, 1, 1), 1, 3))
})
test_that("'ftop_fgm_iso' of vector returns a vector/matrix of expected size", {
  expect_equal(round(ftop_fgm_iso(fitness = c(-1.5, 0), n = 3, maxfitness = 1, alpha = 1/2, Q = 2), digits = 6),
               rbind(c(-2.236068, 0, 0),
                     c(-1.414214, 0, 0)))
  expect_equal(round(ftop_fgm_iso(fitness = c(-1.5, 0), n = 1, maxfitness = 1, alpha = 1/2, Q = 2), digits = 6),
               c(-2.236068, -1.414214))
  expect_equal(round(ftop_fgm_iso(fitness = 0, n = 1, maxfitness = 1, alpha = 1/2, Q = 2), digits = 6),
              - 1.414214)
  expect_equal(round(ftop_fgm_iso(fitness = matrix(seq(0.1, 0.4, 0.1), 2, 2), n = 2, maxfitness = 1, alpha = 1/2, Q = 2), digits = 6),
               rbind(c(-1.341641, 0),
                     c(-1.264911, 0),
                     c(-1.183216, 0),
                     c(-1.095445, 0)))
})
