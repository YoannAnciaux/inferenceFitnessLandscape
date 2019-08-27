test_that("'generate_selected_mutation' of multiple mutants in multiple dimensions returns multiple coselected mutant with increasing fitness higher than the wt", {
  nb_mut <- 5
  n <- 3
  pheno_wt <- c(-1, numeric(n-1))
  ps <- generate_coselected_mutation(nb_mut = nb_mut, n = n, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = nb_mut, ncols = n, null.ok = F)
  pheno_opt <-numeric(n)
  pheno_mutant <- array(dim = dim(ps))
  x <- pheno_wt
  for(mut in 1:nb_mut){
    x <- ps[mut, ] + x
    pheno_mutant[mut, ] <- x
  }
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  fitness_mutant <- ptof_fgm_iso(phenotype = pheno_mutant, maxfitness = 1)

  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_equal(sort(dist_to_opt_mutant, decreasing = TRUE), dist_to_opt_mutant)
  expect_true(all(fitness_mutant > ptof_fgm_iso(matrix(pheno_wt, 1, n), maxfitness = 1)))
  expect_equal(sort(fitness_mutant, decreasing = FALSE), fitness_mutant)
})
test_that("'generate_selected_mutation' of a single mutant in multiple dimensions returns a single selected mutant with fitness higher than the wt", {
  nb_mut <- 1
  n <- 3
  pheno_wt <- c(-1, numeric(n-1))
  ps <- generate_coselected_mutation(nb_mut = nb_mut, n = n, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = nb_mut, ncols = n, null.ok = F)
  pheno_opt <-numeric(n)
  pheno_mutant <- array(dim = dim(ps))
  x <- pheno_wt
  for(mut in 1:nb_mut){
    x <- ps[mut, ] + x
    pheno_mutant[mut, ] <- x
  }
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  fitness_mutant <- ptof_fgm_iso(phenotype = pheno_mutant, maxfitness = 1)

  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_equal(sort(dist_to_opt_mutant, decreasing = TRUE), dist_to_opt_mutant)
  expect_true(all(fitness_mutant > ptof_fgm_iso(matrix(pheno_wt, 1, n), maxfitness = 1)))
  expect_equal(sort(fitness_mutant, decreasing = FALSE), fitness_mutant)
})
test_that("'generate_selected_mutation' of multiple mutants in a single dimension returns multiple coselected mutant with increasing fitness higher than the wt", {
  nb_mut <- 5
  n <- 1
  pheno_wt <- c(-1, numeric(n-1))
  ps <- generate_coselected_mutation(nb_mut = nb_mut, n = n, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = nb_mut, ncols = n, null.ok = F)
  pheno_opt <-numeric(n)
  pheno_mutant <- array(dim = dim(ps))
  x <- pheno_wt
  for(mut in 1:nb_mut){
    x <- ps[mut, ] + x
    pheno_mutant[mut, ] <- x
  }
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  fitness_mutant <- ptof_fgm_iso(phenotype = pheno_mutant, maxfitness = 1)

  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_equal(sort(dist_to_opt_mutant, decreasing = TRUE), dist_to_opt_mutant)
  expect_true(all(fitness_mutant > ptof_fgm_iso(matrix(pheno_wt, 1, n), maxfitness = 1)))
  expect_equal(sort(fitness_mutant, decreasing = FALSE), fitness_mutant)
})
test_that("'generate_selected_mutation' of a single mutant in a single dimension returns a single selected mutant with fitness higher than the wt", {
  nb_mut <- 1
  n <- 1
  pheno_wt <- c(-1, numeric(n-1))
  ps <- generate_coselected_mutation(nb_mut = nb_mut, n = n, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = nb_mut, ncols = n, null.ok = F)
  pheno_opt <-numeric(n)
  pheno_mutant <- array(dim = dim(ps))
  x <- pheno_wt
  for(mut in 1:nb_mut){
    x <- ps[mut, ] + x
    pheno_mutant[mut, ] <- x
  }
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  fitness_mutant <- ptof_fgm_iso(phenotype = pheno_mutant, maxfitness = 1)

  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_equal(sort(dist_to_opt_mutant, decreasing = TRUE), dist_to_opt_mutant)
  expect_true(all(fitness_mutant > ptof_fgm_iso(matrix(pheno_wt, 1, n), maxfitness = 1)))
  expect_equal(sort(fitness_mutant, decreasing = FALSE), fitness_mutant)
})
test_that("'generate_selected_mutation' of multiple mutants in multiple dimensions with non default parameters returns multiple coselected mutant with increasing fitness higher than the wt", {
  nb_mut <- 5
  n <- 3
  pheno_wt <- c(-1, numeric(n-1))
  ps <- generate_coselected_mutation(nb_mut = nb_mut, n = n, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt,
                                     alpha = 1, Q = 0.5, m = 2, nb_mut_rand = 10^5)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = nb_mut, ncols = n, null.ok = F)
  pheno_opt <-numeric(n)
  pheno_mutant <- array(dim = dim(ps))
  x <- pheno_wt
  for(mut in 1:nb_mut){
    x <- ps[mut, ] + x
    pheno_mutant[mut, ] <- x
  }
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  fitness_mutant <- ptof_fgm_iso(phenotype = pheno_mutant, maxfitness = 1,
                                 alpha = 1, Q = 0.5)

  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_equal(sort(dist_to_opt_mutant, decreasing = TRUE), dist_to_opt_mutant)
  expect_true(all(fitness_mutant > ptof_fgm_iso(matrix(pheno_wt, 1, n), maxfitness = 1, alpha = 1, Q = 0.5)))
  expect_equal(sort(fitness_mutant, decreasing = FALSE), fitness_mutant)
})
