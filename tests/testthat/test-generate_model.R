test_that("'check_genotype_table' of genotype_table of the 8 combinations of 3 mutants with an error returns an error message", {
  x <- as.matrix(expand.grid(rep(list(0:1), 3)))
  x <- x[sample(1:dim(x)[1]),]
  x[1,1] <- 8
  expect_equal(check_genotype_table(x), "The rows must contains the 8 combinations of the 3 column(s) in a 0/1 format")
})
test_that("'check_genotype_table' of genotype_table of the 8 combinations of 3 mutants returns TRUE", {
  x <- as.matrix(expand.grid(rep(list(0:1), 3)))
  x <- x[sample(1:dim(x)[1]),]
  expect_true(check_genotype_table(x))
})
test_that("'fitness_mutant_genotype' of a 2 mutants genotype table returns the expected genotypes' fitnesses in the expected order", {
  geno_table <- as.matrix(expand.grid(rep(list(0:1), 2)))
  pheno_wt <- ftop_fgm_iso(fitness = 0, n = 3, maxfitness = 1)
  pheno_rand_mut_effect <- rbind(c(0.5,0,0),
                                 c(-0.5,0,0))
  pheno_mutant <- rbind(pheno_wt,
                        pheno_rand_mut_effect[1,] + pheno_wt,
                        pheno_rand_mut_effect[2,] + pheno_wt,
                        pheno_rand_mut_effect[1,] + pheno_rand_mut_effect[2,] + pheno_wt)
  fitnesses_expected <- ptof_fgm_iso(phenotype = pheno_mutant, maxfitness = 1)
  fitnesses_res <- fitness_mutant_genotype(pheno_mut_effect = pheno_rand_mut_effect, pheno_wt = pheno_wt, geno_table = geno_table, maxfitness = 1)
  expect_equal(fitnesses_res, unname(fitnesses_expected))
})
test_that("'generate_model' of 'fgmrmut' returns the expected number of random mutants and the wt",{
  empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), 3))), seq(0.1, 0.8, 0.1)))
  t <- generate_model(empirical_fl = empirical_fl, model_type = "fgmrmut")
  r <- t(3, 0.1, 1, 1/2, 2, 3)
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = dim(empirical_fl)[1], null.ok = F)
  expect_equal(r[1], empirical_fl[1, dim(empirical_fl)[2]])
})
test_that("'generate_model' of 'fgmsmut' returns the expected number of selected mutants and the wt",{
  empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), 3))), seq(0.1, 0.8, 0.1)))
  t <- generate_model(empirical_fl = empirical_fl, model_type = "fgmsmut", fun_args = list(nb_mut_rand = 10^6))
  r <- t(3, 0.1, 1, 1/2, 2, 3)
  checkmate::expect_numeric(round(r, digits = 10), lower = round(empirical_fl[1, dim(empirical_fl)[2]], digits = 10), finite = T, all.missing = T, len = dim(empirical_fl)[1], null.ok = F)
  if (!anyNA(r)) {
    expect_equal(r[1], empirical_fl[1, dim(empirical_fl)[2]])
  }
})
test_that("'generate_model' of 'fgmcsmut' returns the expected number of coselected mutants and the wt",{

  empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), 3))), seq(0.1, 0.8, 0.1)))
  t <- generate_model(empirical_fl = empirical_fl, model_type = "fgmcsmut", fun_args = list(nb_mut_rand = 10^6))
  r <- t(3, 0.1, 1, 1/2, 2, 3)
  checkmate::expect_numeric(r, finite = T, all.missing = T, len = dim(empirical_fl)[1], null.ok = F)
  checkmate::expect_numeric(round(r[c(2, 4, 8)], digits = 10), lower = round(empirical_fl[1, dim(empirical_fl)[2]], digits = 10), finite = T, all.missing = T, sorted = T, null.ok = F)
  if (!anyNA(r)) {
    expect_equal(r[1], empirical_fl[1, dim(empirical_fl)[2]])
  }
})
