path2_fl_generate <- system.file("bin", "fl_generate", package = "inferenceFitnessLandscape")
nb_mut <- 3
empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), nb_mut))), seq(0.1, 0.1*2^nb_mut, 0.1)))
genotype_table <- empirical_fl[, -(nb_mut + 1)]
fitness_wt <- empirical_fl[1, nb_mut + 1]

test_that("'generate_model' of 'magellanFix', returns the expected fitnesses", {
  m <- generate_model(empirical_fl, "magellanFix")
  r <- m(f = 1)
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = 2^nb_mut)
  expect_equal(r, rep(fitness_wt, 2^nb_mut))
})
test_that("'generate_model' of 'magellanMult', returns the expected fitnesses", {
  m <- generate_model(empirical_fl, "magellanMult")
  r <- m(s = 0.1, S = 1, d = 1)
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = 2^nb_mut)
  expect_equal(r[1], fitness_wt)
})
test_that("'generate_model' of 'magellanHoC', returns the expected fitnesses", {
  m <- generate_model(empirical_fl, "magellanHoC")
  r <- m(H = 1)
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = 2^nb_mut)
  expect_equal(r[1], fitness_wt)
})
test_that("'generate_model' of 'magellanNK', returns the expected fitnesses", {
  m <- generate_model(empirical_fl, "magellanNK")
  r <- m(K = 2, r = T)
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = 2^nb_mut)
  expect_equal(r[1], fitness_wt)
})
test_that("'generate_model' of 'magellanIsing', returns the expected fitnesses", {
  m <- generate_model(empirical_fl, "magellanIsing")
  r <- m(i = 0.2, I = 1, c = T)
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = 2^nb_mut)
  expect_equal(r[1], fitness_wt)
  expect_equal(r[2^nb_mut], fitness_wt)
})
test_that("'generate_model' of 'magellanEggBox', returns the expected fitnesses", {
  m <- generate_model(empirical_fl, "magellanEggBox")
  r <- m(e = 0.1, E = 1)
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = 2^nb_mut)
  expect_equal(r[1], fitness_wt)
})
test_that("'generate_model' of 'magellanOptimum', returns the expected fitnesses", {
  m <- generate_model(empirical_fl, "magellanOptimum")
  r <- m(o = 0.3, O = 0.1, p = 1, P = 2)
  checkmate::expect_numeric(r, finite = T, any.missing = F, len = 2^nb_mut)
  expect_equal(r[1], fitness_wt)
})
test_that("'fl_generate_magellan' returns a fitness landscape of the expected type and dimensions", {
  command <- paste(path2_fl_generate, "-L", "-f", 1, nb_mut, 2)
  fl <- fl_generate_magellan(nb_mut, genotype_table, command)
  checkmate::expect_matrix(fl, mode = "numeric", nrow = nrow(genotype_table), ncol = nb_mut + 1)
  checkmate::expect_numeric(fl[, -(nb_mut + 1)], lower = 0, upper = 1, any.missing = F)
})
test_that("'equalize_row_order' returns a fitness_landscape in the expected order", {
  rand <- sample(1:nrow(genotype_table))
  genotype_table_rand <- genotype_table[rand,]
  r <- equalize_row_order(genotype_table_rand, empirical_fl)
  expect_equal(r[, -(nb_mut + 1)], genotype_table_rand)
  expect_equal(r[, nb_mut + 1], empirical_fl[rand, nb_mut + 1])
})
test_that("'scale_to_wt_fitness' returns a vector of fitness scaled to the fitness_wt", {
  fit_wt <- 2
  target <- empirical_fl[, nb_mut + 1] - fitness_wt + fit_wt
  r <- scale_to_wt_fitness(nb_mut, genotype_table, fit_wt, empirical_fl)
  expect_equal(r, target)
})
rm(list = c("path2_fl_generate", "nb_mut", "genotype_table", "empirical_fl", "fitness_wt"))
