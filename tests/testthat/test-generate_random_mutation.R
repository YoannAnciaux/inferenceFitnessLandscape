test_that("'generate_random_mutation' of multiple mutants in multiple phenotypic dimensions returns the expected matrix", {
  checkmate::expect_matrix(generate_random_mutation(nb_mut = 3, n = 4, lambda = 0.1), mode = "numeric", any.missing = F, nrows = 3, ncols = 4)
})
test_that("'generate_random_mutation' of multiple mutants in a single phenotypic dimension returns the expected matrix", {
  checkmate::expect_matrix(generate_random_mutation(nb_mut = 3, n = 1, lambda = 0.1), mode = "numeric", any.missing = F, nrows = 3, ncols = 1)
})
test_that("'generate_random_mutation' of a single mutant in multiple phenotypic dimensions returns the expected matrix", {
  checkmate::expect_matrix(generate_random_mutation(nb_mut = 1, n = 3, lambda = 0.1), mode = "numeric", any.missing = F, nrows = 1, ncols = 3)
})
test_that("'generate_random_mutation' of a single mutant in a single phenotypic dimension returns the expected matrix", {
  checkmate::expect_matrix(generate_random_mutation(nb_mut = 1, n = 1, lambda = 0.1), mode = "numeric", any.missing = F, nrows = 1, ncols = 1)
})
test_that("'generate_random_mutation' with restricted pleiotropy returns m '0' at random positions per row", {
  set.seed(1)
  mut <- generate_random_mutation(nb_mut = 5, n = 3, lambda = 0.1, m = 2)
  checkmate::expect_matrix(mut, mode = "numeric", any.missing = F, nrows = 5, ncols = 3)
  zero_nb_and_pos <- apply(X = mut, MARGIN = 1, function(r) {p <- which(r==0); c(length(p), p)})
  expect_equal(zero_nb_and_pos[1,], rep(3-2, 5))
  expect_equal(unique(zero_nb_and_pos[2,]), c(2,1,3))
})
