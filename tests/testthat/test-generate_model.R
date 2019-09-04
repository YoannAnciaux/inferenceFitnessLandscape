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
  genotype_table <- as.matrix(expand.grid(rep(list(0:1), 2)))
  pheno_wt <- ftop_fgm_iso(fitness = 0, n = 3, maxfitness = 1)
  pheno_rand_mut_effect <- rbind(c(0.5,0,0),
                                 c(-0.5,0,0))
  pheno_mutant <- rbind(pheno_wt,
                        pheno_rand_mut_effect[1,] + pheno_wt,
                        pheno_rand_mut_effect[2,] + pheno_wt,
                        pheno_rand_mut_effect[1,] + pheno_rand_mut_effect[2,] + pheno_wt)
  fitnesses_expected <- ptof_fgm_iso(phenotype = pheno_mutant, maxfitness = 1)
  fitnesses_res <- fitness_mutant_genotype(pheno_mut_effect = pheno_rand_mut_effect, pheno_wt = pheno_wt, genotype_table = genotype_table, maxfitness = 1)
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
test_that("'generate_model' of 'fgmrmut2env' returns the expected fitnesses of random mutants in the new env",{
  #### parameters ####
  empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), 3))), seq(0.1, 0.8, 0.1)))
  genotype_table <- empirical_fl[, 1:3]
  lambda <- 0.1; maxfitness <- 1; alpha <- 1; Q <- 1; theta = pi/4
  fun_args <-  list(fitness_wt_ref = -1, n_ref = 2, lambda_ref = lambda,
                    maxfitness_ref = maxfitness, alpha_ref = alpha, Q_ref = Q, m_ref = 2)
  #### manual compuation of fitness in new env ####
  pheno_wt_ref <- ftop_fgm_iso(fitness = fun_args$fitness_wt_ref,
                               n = fun_args$n_ref, maxfitness = fun_args$maxfitness_ref, alpha = fun_args$alpha_ref, Q = fun_args$Q_ref)
  set.seed(1); pheno_mut_effect <- generate_random_mutation(nb_mut = 3,
                                                            n = fun_args$n_ref,
                                                            lambda = fun_args$lambda_ref,
                                                            m = fun_args$m_ref)
  target <- fitness_mutant_genotype(pheno_mut_effect, pheno_wt_ref, genotype_table, maxfitness, alpha, Q, pheno_opt = c(cos(theta), sin(theta)) * (1 / alpha * (maxfitness-empirical_fl[1,4]))^(1/Q) + pheno_wt_ref)
  #### current ####
  t <- generate_model(empirical_fl = empirical_fl, model_type = "fgmrmut2env",
                      fun_args = fun_args)
  set.seed(1); r <- t(lambda, maxfitness, alpha, Q, theta)
  #### tests ####
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = dim(empirical_fl)[1]*2, null.ok = F)
  expect_equal(r[9], empirical_fl[1, dim(empirical_fl)[2]])
  expect_equal(r[1], fun_args$fitness_wt_ref)
  expect_equal(r[9:16], target)
})
test_that("'generate_model' of 'fgmsmut2env' returns the expected fitnesses of selected mutants in the new env",{
  #### parameters ####
  empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), 3))), seq(0.1, 0.8, 0.1)))
  genotype_table <- empirical_fl[, 1:3]
  lambda <- 0.1; maxfitness <- 1; alpha <- 1; Q <- 1; theta = pi/4
  fun_args <-  list(fitness_wt_ref = -1, n_ref = 2, lambda_ref = lambda,
                    maxfitness_ref = maxfitness, alpha_ref = alpha, Q_ref = Q, m_ref = 2,
                    nb_mut_rand = 10^4)
  #### manual compuation of fitness in new env ####
  pheno_wt_ref <- ftop_fgm_iso(fitness = fun_args$fitness_wt_ref,
                               n = fun_args$n_ref, maxfitness = fun_args$maxfitness_ref, alpha = fun_args$alpha_ref, Q = fun_args$Q_ref)
  set.seed(1); pheno_mut_effect <- do.call(generate_selected_mutation, list(nb_mut = 3,
                                                                                n = fun_args$n_ref,
                                                                                lambda = fun_args$lambda_ref,
                                                                                maxfitness = fun_args$maxfitness_ref,
                                                                                pheno_wt = pheno_wt_ref,
                                                                                alpha = fun_args$alpha_ref,
                                                                                Q = fun_args$Q_ref,
                                                                                m = fun_args$m_ref,
                                                                                nb_mut_rand = fun_args$nb_mut_rand))
  target <- fitness_mutant_genotype(pheno_mut_effect, pheno_wt_ref, genotype_table, maxfitness, alpha, Q, pheno_opt = c(cos(theta), sin(theta)) * (1 / alpha * (maxfitness-empirical_fl[1,4]))^(1/Q) + pheno_wt_ref)
  #### current ####
  t <- generate_model(empirical_fl = empirical_fl, model_type = "fgmsmut2env",
                      fun_args = fun_args)
  set.seed(1); r <- t(lambda, maxfitness, alpha, Q, theta)
  #### tests ####
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = dim(empirical_fl)[1]*2, null.ok = F)
  expect_equal(r[9], empirical_fl[1, dim(empirical_fl)[2]])
  expect_equal(r[1], fun_args$fitness_wt_ref)
  expect_equal(r[9:16], target)
})
test_that("'generate_model' of 'fgmcsmut2env' returns the expected fitnesses of coselected mutants in the new env",{
  #### parameters ####
  empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), 3))), seq(0.1, 0.8, 0.1)))
  genotype_table <- empirical_fl[, 1:3]
  lambda <- 0.1; maxfitness <- 1; alpha <- 1; Q <- 1; theta = pi/4
  fun_args <-  list(fitness_wt_ref = -1, n_ref = 2, lambda_ref = lambda,
                    maxfitness_ref = maxfitness, alpha_ref = alpha, Q_ref = Q, m_ref = 2,
                    nb_mut_rand = 10^4)
  #### manual compuation of fitness in new env ####
  pheno_wt_ref <- ftop_fgm_iso(fitness = fun_args$fitness_wt_ref,
                               n = fun_args$n_ref, maxfitness = fun_args$maxfitness_ref, alpha = fun_args$alpha_ref, Q = fun_args$Q_ref)
  set.seed(1); pheno_mut_effect <- do.call(generate_coselected_mutation, list(nb_mut = 3,
                                                                            n = fun_args$n_ref,
                                                                            lambda = fun_args$lambda_ref,
                                                                            maxfitness = fun_args$maxfitness_ref,
                                                                            pheno_wt = pheno_wt_ref,
                                                                            alpha = fun_args$alpha_ref,
                                                                            Q = fun_args$Q_ref,
                                                                            m = fun_args$m_ref,
                                                                            nb_mut_rand = fun_args$nb_mut_rand))
  target <- fitness_mutant_genotype(pheno_mut_effect, pheno_wt_ref, genotype_table, maxfitness, alpha, Q, pheno_opt = c(cos(theta), sin(theta)) * (1 / alpha * (maxfitness-empirical_fl[1,4]))^(1/Q) + pheno_wt_ref)
  #### current ####
  t <- generate_model(empirical_fl = empirical_fl, model_type = "fgmcsmut2env",
                      fun_args = fun_args)
  set.seed(1); r <- t(lambda, maxfitness, alpha, Q, theta)
  #### tests ####
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = dim(empirical_fl)[1]*2, null.ok = F)
  expect_equal(r[9], empirical_fl[1, dim(empirical_fl)[2]])
  expect_equal(r[1], fun_args$fitness_wt_ref)
  expect_equal(r[9:16], target)
})
test_that("'pos_new_env' of theta pi/4 returns the expected position", {
  pheno_wt_ref <- matrix(c(-1, 0), 1, 2)
  opt_new <- c(cospi(1/4) * 1.5, sinpi(1/4) * 1.5) + pheno_wt_ref
  fitness_wt_new_env <- ptof_fgm_iso(phenotype = pheno_wt_ref, maxfitness = 1, alpha = 1/2, Q = 2, pheno_opt = opt_new)
  expect_equal(matrix(pos_new_env(pheno_wt_ref, fitness_wt_new_env,
                                  maxfitness = 1, alpha = 1/2, Q = 2, theta = pi/4), 1, 2), opt_new)
})
test_that("'pos_new_env' of alpha 2 and Q 1 returns the expected position", {
  pheno_wt_ref <- matrix(c(-1, 0), 1, 2)
  opt_new <- c(cospi(1/4) * 1.5, sinpi(1/4) * 1.5) + pheno_wt_ref
  fitness_wt_new_env <- ptof_fgm_iso(phenotype = pheno_wt_ref, maxfitness = 1, alpha = 2, Q = 1, pheno_opt = opt_new)
  expect_equal(matrix(pos_new_env(pheno_wt_ref, fitness_wt_new_env,
                                  maxfitness = 1, alpha = 2, Q = 1, theta = pi/4), 1, 2), opt_new)
})
test_that("'pos_new_env' of theta pi*5/4 returns the expected position", {
  pheno_wt_ref <- matrix(c(-1, 0), 1, 2)
  opt_new <- c(cospi(5/4) * 1.5, sinpi(5/4) * 1.5) + pheno_wt_ref
  fitness_wt_new_env <- ptof_fgm_iso(phenotype = pheno_wt_ref, maxfitness = 1, alpha = 1/2, Q = 2, pheno_opt = opt_new)
  expect_equal(matrix(pos_new_env(pheno_wt_ref, fitness_wt_new_env,
                                  maxfitness = 1, alpha = 1/2, Q = 2, theta = pi* 5/4), 1, 2), opt_new)
})
test_that("'pos_new_env' of for a single dim returns the expected position in a singke dim", {
  pheno_wt_ref <- matrix(-1, 1, 1)
  opt_new <- cospi(0) * 1 + pheno_wt_ref
  fitness_wt_new_env <- ptof_fgm_iso(phenotype = pheno_wt_ref, maxfitness = 1, alpha = 1/2, Q = 2, pheno_opt = opt_new)
  expect_equal(matrix(pos_new_env(pheno_wt_ref, fitness_wt_new_env,
                                  maxfitness = 1, alpha = 1/2, Q = 2, theta = 0), 1, 1), opt_new)
})
