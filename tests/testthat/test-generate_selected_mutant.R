test_that("'add_mut_to_pheno' of multiple mutations in multiple dimensions returns the expected mutation effects", {
  pheno_mut_effect <- rbind(c(1, 1), c(-1, 2))
  pheno_wt <- c(1,2)
  expect_equal(dim(add_mut_to_pheno(pheno_mut_effect, pheno_wt)),
               dim(pheno_mut_effect))
  expect_equal(add_mut_to_pheno(pheno_mut_effect, pheno_wt),
               rbind(pheno_mut_effect[1,] + pheno_wt,
                     pheno_mut_effect[2,] + pheno_wt))
})
test_that("'add_mut_to_pheno' of multiple mutations in a single dimension returns the expected mutation effects", {
  pheno_mut_effect <- rbind(c(1), c(-1))
  pheno_wt <- c(2)
  expect_equal(dim(add_mut_to_pheno(pheno_mut_effect, pheno_wt)),
               dim(pheno_mut_effect))
  expect_equal(add_mut_to_pheno(pheno_mut_effect, pheno_wt),
               rbind(pheno_mut_effect[1,] + pheno_wt,
                     pheno_mut_effect[2,] + pheno_wt))
})
test_that("'add_mut_to_pheno' of a single mutation in multiple dimensions returns the expected mutation effects", {
  pheno_mut_effect <- matrix(c(1,1), 1, 2)
  pheno_wt <- c(1,2)
  expect_equal(length(add_mut_to_pheno(pheno_mut_effect, pheno_wt)),
               length(pheno_mut_effect))
  expect_equal(add_mut_to_pheno(pheno_mut_effect, pheno_wt),
               pheno_mut_effect + pheno_wt)
})
test_that("'add_mut_to_pheno' of a single mutation in a single dimension returns the expected mutation effect", {
  pheno_mut_effect <- matrix(1, 1, 1)
  pheno_wt <- c(2)
  expect_equal(length(add_mut_to_pheno(pheno_mut_effect, pheno_wt)),
               length(pheno_mut_effect))
  expect_equal(add_mut_to_pheno(pheno_mut_effect, pheno_wt),
               pheno_mut_effect + pheno_wt)
})
test_that("'generate_selected_mutation' of multiple mutants in multiple dimensions returns multipleselected mutants with fitness higher than the wt", {
  pheno_wt <- c(-1, 0, 0)
  ps <- generate_selected_mutation(nb_mut = 5, n = 3, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = 5, ncols = 3, null.ok = F)
  pheno_opt <-numeric(3)
  pheno_mutant <- add_mut_to_pheno(ps, pheno_wt)
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_true(all(ptof_fgm_iso(phenotype = matrix(pheno_wt, nrow = 1, ncol = 3), maxfitness = 1) < ptof_fgm_iso(phenotype = ps, maxfitness = 1)))
})
test_that("'generate_selected_mutation' of a single mutant in multiple dimensions returns a single selected mutant with fitness higher than the wt", {
  pheno_wt <- c(-1, 0, 0)
  ps <- generate_selected_mutation(nb_mut = 1, n = 3, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = 1, ncols = 3, null.ok = F)
  pheno_opt <-numeric(3)
  pheno_mutant <- add_mut_to_pheno(ps, pheno_wt)
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_true(all(ptof_fgm_iso(phenotype = matrix(pheno_wt, nrow = 1, ncol = 3), maxfitness = 1) < ptof_fgm_iso(phenotype = ps, maxfitness = 1)))
})
test_that("'generate_selected_mutation' of multiple mutants in a single dimension returns multiple mutants with fitness higher than the wt", {
  pheno_wt <- c(-1)
  ps <- generate_selected_mutation(nb_mut = 5, n = 1, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = 5, ncols = 1, null.ok = F)
  pheno_opt <-numeric(1)
  pheno_mutant <- add_mut_to_pheno(ps, pheno_wt)
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_true(all(ptof_fgm_iso(phenotype = matrix(pheno_wt, nrow = 1, ncol = 1), maxfitness = 1) < ptof_fgm_iso(phenotype = ps, maxfitness = 1)))
})
test_that("'generate_selected_mutation' of a single mutant in a single dimension returns a single mutant with fitness higher than the wt", {
  pheno_wt <- c(-1)
  ps <- generate_selected_mutation(nb_mut = 1, n = 1, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = 1, ncols = 1, null.ok = F)
  pheno_opt <-numeric(1)
  pheno_mutant <- add_mut_to_pheno(ps, pheno_wt)
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_true(all(ptof_fgm_iso(phenotype = matrix(pheno_wt, nrow = 1, ncol = 1), maxfitness = 1) < ptof_fgm_iso(phenotype = ps, maxfitness = 1)))
})
test_that("'generate_selected_mutation' of multiple mutants in multiple dimensions with non default values returns multiple selected mutants with fitness higher than the wt", {
  pheno_wt <- c(-1, 0, 0)
  ps <- generate_selected_mutation(nb_mut = 5, n = 3, lambda = 0.1, maxfitness = 1, pheno_wt = pheno_wt,
                                   alpha = 1, Q = 0.5, m = 3-1, nb_mut_rand = 10^5)
  checkmate::expect_matrix(ps, mode = "numeric", all.missing = T, nrows = 5, ncols = 3, null.ok = F)
  pheno_opt <-numeric(3)
  pheno_mutant <- add_mut_to_pheno(ps, pheno_wt)
  dist_to_opt_mutant <- apply(X = pheno_mutant,
                              MARGIN = 1,
                              FUN = function(p) {sum((pheno_opt - p)^2)^(1/2)})
  dist_to_opt_wt <- sum((pheno_opt - pheno_wt)^2)^(1/2)
  expect_true(all(dist_to_opt_mutant < dist_to_opt_wt))
  expect_true(all(ptof_fgm_iso(phenotype = matrix(pheno_wt, nrow = 1, ncol = 3), maxfitness = 1, alpha = 1, Q = 0.5) < ptof_fgm_iso(phenotype = ps, maxfitness = 1,alpha = 1, Q = 0.5)))
})
