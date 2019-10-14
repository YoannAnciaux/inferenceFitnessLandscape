#' Closure function which returns a function generating \code{nb_mut} fitnesses
#' of mutant(s) (based on a wild-type of fitness given in the last element of the
#' first row of \code{empirical_fl}). The model used for generating the mutant can
#' chosen with \code{model_type}.
#' The function returned then take a list of varying length as argument depending
#' on \code{model_type}.
#'
#' @param empirical_fl A matrix corresponding to an empirical fitness landscape.
#' The first ncol-1 columns correspond to a genotype table. The rows correspond
#' to different genotypes and the columns to the mutations that are considered for
#' these genotypes. A genotype (at a certain row) has a given mutation when there
#' is a 1 in the corresponding column.
#' The column to of empirical_fl correspond to the fitness of each genotype. Only
#' the fitness of the wt (i.e. the row in the genotype table with only zeros) is
#' used.
#' @param model_type A character corresponding to one of the following implemented
#' models :
#' \describe{
#'   \item{"fgmrmut"}{random mutation in isotropic FGM. See \code{\link{model_fgmrmut}}}
#'   \item{"fgmsmut"}{selected mutation in isotropic FGM. See \code{\link{model_fgmsmut}}}
#'   \item{"fgmcsmut"}{coselected mutation in isotropic FGM. See \code{\link{model_fgmcsmut}}}
#'   \item{"fgmrmut2env"}{random mutation in isotropic FGM with a refence environment and a environment. See \code{\link{model_fgmrmut_2env}}}
#'   \item{"fgmsmut2env"}{selected mutation in isotropic FGM with a refence environment and a environment. See \code{\link{model_fgmsmut_2env}}}
#'   \item{"fgmcsmut2env"}{coselected mutation in isotropic FGM with a refence environment and a environment. See \code{\link{model_fgmcsmut_2env}}}
#'   \item{"magellanFix"}{Fixed fitness for all genotypes. See \code{\link{model_Fix}}}
#'   \item{"magellanMult"}{Additive effect of mutations. The fitnesses of all mutations are independent. See \code{\link{model_Mult}}}
#'   \item{"magellanHoC"}{House of Card model. The fitnesses of all genotypes are iid. See \code{\link{model_HoC}}}
#'   \item{"magellanNK"}{NK model from Kauffman et al. (1988). Each locus interacts with K other loci, that can be its neighbors or can be chosen randomly. Fitness are drawn in uniform [0,1]. See \code{\link{model_NK}}}
#'   \item{"magellanIsing"}{all loci are arranged sequentially, and each locus interacts with its physical neighbors. The last and the first loci will interact only if 'c' is TRUE. For each pair of interacting loci, there is an associated cost if both alleles are not identical (and therefore 'compatible'). See \code{\link{model_Ising}}}
#'   \item{"magellanEggBox"}{Each locus is either high or low fitness, with a systematic change between each neighbor. See \code{\link{model_EggBox}}}
#'   \item{"magellanOptimum"}{Each mutated locus produces a contribution to the fitness according to a production (p,P or uniform) distribution and compared to an optimum with a (o,O or uniform) distribution. See \code{\link{model_EggBox}}}
#' }
#' @param fun_args List of argument for a given \code{model_type}.
#' Argument *_ref are mandatory parameters for the environment of reference in models with two environments
#' \describe{
#'   \item{nb_mut_rand}{optionnal for "fgmsmut", "fgmcsmut", "fgmsmut2env" and "fgmcsmut2env". See \code{\link{model_fgmrmut_2env}}}
#'   \item{fitness_wt_ref}{mandatory for "fgmsmut2env" and "fgmcsmut2env". See \code{\link{model_fgmrmut_2env}}}
#'   \item{n_ref}{mandatory for "fgmrmut2env", "fgmsmut2env" and "fgmcsmut2env". See \code{\link{model_fgmrmut_2env}}}
#'   \item{lambda_ref}{mandatory for "fgmrmut2env", "fgmsmut2env" and "fgmcsmut2env". See \code{\link{model_fgmrmut_2env}}}
#'   \item{maxfitness_ref}{mandatory for "fgmrmut2env", "fgmsmut2env" and "fgmcsmut2env". See \code{\link{model_fgmrmut_2env}}}
#'   \item{alpha_ref}{mandatory for ""fgmrmut2env", fgmsmut2env" and "fgmcsmut2env". See \code{\link{model_fgmrmut_2env}}}
#'   \item{Q_ref}{mandatory for "fgmrmut2env", "fgmsmut2env" and "fgmcsmut2env". See \code{\link{model_fgmrmut_2env}}}
#'   \item{m_ref}{mandatory for "fgmrmut2env", "fgmsmut2env" and "fgmcsmut2env". See \code{\link{model_fgmrmut_2env}}}
#' }
#' @param ... Extra arguments which will be passed to \code{\link[utils]{read.table}} if an empirical
#' fitness landscape is provided as a file in \code{empirical_fl}
#' @return A function for generating mutants following the model from \code{model_type}
#' and the genotype table and the fitness of the wild type given in \code{empirical_fl}.
#' @examples
#' #random mutation in isotropic FGM
#' empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), 3))), seq(0.1, 0.8, 0.1)))
#' model <- generate_model(empirical_fl = empirical_fl, model_type = "fgmrmut")
#' model(3, 0.1, 1, 1/2, 2, 3)
#' #selected mutation in isotropic FGM
#' empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), 3))), seq(0.1, 0.8, 0.1)))
#' model <- generate_model(empirical_fl = empirical_fl, model_type = "fgmsmut",
#'                         fun_args = list(nb_mut_rand = 10^5))
#' model(3, 0.1, 1, 1/2, 2, 3)
#' #coselected mutation in isotropic FGM
#' empirical_fl <- unname(cbind(as.matrix(expand.grid(rep(list(0:1), 3))), seq(0.1, 0.8, 0.1)))
#' model <- generate_model(empirical_fl = empirical_fl, model_type = "fgmcsmut",
#'                         fun_args = list(nb_mut_rand = 10^5))
#' model(3, 0.1, 1, 1/2, 2, 3)
#' @export
generate_model <- function(empirical_fl, model_type, fun_args = list(), ...) {

  #### check arguments ####
  arg_required <- c("empirical_fl", "model_type")
  arg_passed <- names(as.list(match.call())[-1])
  coll1 <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll1$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::reportAssertions(coll1)
  force(empirical_fl)
  force(model_type)
  force(fun_args)
  coll <- checkmate::makeAssertCollection()
  if (is.character(empirical_fl)) {
    checkmate::assert_file_exists(empirical_fl, access = "r", add = coll)
    if (checkmate::test_file_exists(empirical_fl, access = "r")) {
      empirical_fl <- as.matrix(unname(utils::read.table(file = empirical_fl, ...)))
    }
  }
  checkmate::assert_matrix(empirical_fl, mode = "numeric", add = coll)
  if (checkmate::test_matrix(empirical_fl, mode = "numeric")) {
    assert_genotype_table(empirical_fl[, -dim(empirical_fl)[2]],
                          .var.name = paste0(checkmate::vname(empirical_fl), " without the last column (fitness)"),
                          add = coll)
  }
  checkmate::assert_choice(model_type, choices = c("fgmrmut", "fgmsmut", "fgmcsmut",
                                                   "fgmrmut2env", "fgmsmut2env", "fgmcsmut2env",
                                                   "magellanFix", "magellanMult", "magellanHoC", "magellanNK", "magellanIsing", "magellanEggBox", "magellanOptimum"),
                           add = coll)
  #when adding new models check all the values from the two following assertions
  checkmate::assert_list(fun_args, any.missing = F, max.len = 8, names = "unique",
                         null.ok = F, add = coll)

  switch (model_type,
          "fgmsmut" = checkmate::assert_subset(names(fun_args), choices = c("nb_mut_rand"), add = coll),
          "fgmcsmut" = checkmate::assert_subset(names(fun_args), choices = c("nb_mut_rand"), add = coll),
          "fgmrmut2env" = checkmate::assert_subset(names(fun_args),
                                                   choices = c("fitness_wt_ref", "n_ref", "lambda_ref", "maxfitness_ref", "alpha_ref", "Q_ref", "m_ref"),
                                                   add = coll),
          "fgmsmut2env" = checkmate::assert_subset(names(fun_args),
                                                   choices = c("fitness_wt_ref", "n_ref", "lambda_ref", "maxfitness_ref", "alpha_ref", "Q_ref", "m_ref", "nb_mut_rand"),
                                                   add = coll),
          "fgmcsmut2env" = checkmate::assert_subset(names(fun_args),
                                                    choices = c("fitness_wt_ref", "n_ref", "lambda_ref", "maxfitness_ref", "alpha_ref", "Q_ref", "m_ref", "nb_mut_rand"),
                                                    add = coll)
  )
  checkmate::reportAssertions(coll)
  nb_mut <- dim(empirical_fl)[2]-1
  genotype_table <- empirical_fl[, -(nb_mut + 1)]
  fitness_wt <- empirical_fl[which(rowSums(genotype_table) == 0), nb_mut + 1]
  # file path to rhe binary fl_generate (from MAGELLAN) in the package directory
  path2_fl_generate <- system.file("bin", "fl_generate", package = "inferenceFitnessLandscape")

  switch (model_type,
          # FGM random mutation 1 environment
          "fgmrmut" = function(n, lambda, maxfitness, alpha, Q, m) {
            model_fgmrmut(nb_mut = nb_mut, genotype_table = genotype_table,
                          fitness_wt =fitness_wt,
                          n = n, lambda = lambda, maxfitness = maxfitness,
                          alpha = alpha, Q = Q, m = m)
          },
          # FGM selected mutation 1 environment
          "fgmsmut" = function(n, lambda, maxfitness, alpha, Q, m) {
            model_fgmsmut(nb_mut = nb_mut, genotype_table = genotype_table,
                          fitness_wt =fitness_wt,
                          n = n, lambda = lambda, maxfitness = maxfitness,
                          alpha = alpha, Q = Q, m = m, fun_args = fun_args)
          },
          # FGM coselected mutation 1 environment
          "fgmcsmut" = function(n, lambda, maxfitness, alpha, Q, m) {
            model_fgmcsmut(nb_mut = nb_mut, genotype_table = genotype_table,
                           fitness_wt =fitness_wt,
                           n = n, lambda = lambda, maxfitness = maxfitness,
                           alpha = alpha, Q = Q, m = m, fun_args = fun_args)
          },
          # FGM random mutation 2 environments
          "fgmrmut2env" = function(lambda, maxfitness, alpha, Q, theta) {
            model_fgmrmut_2env(nb_mut = nb_mut, genotype_table = genotype_table,
                               fitness_wt_new_env =fitness_wt,
                               lambda = lambda, maxfitness = maxfitness,
                               alpha = alpha, Q = Q, theta = theta, fun_args = fun_args)
          },
          # FGM selected mutation 2 environments
          "fgmsmut2env" = function(lambda, maxfitness, alpha, Q, theta) {
            model_fgmsmut_2env(nb_mut = nb_mut, genotype_table = genotype_table,
                               fitness_wt_new_env =fitness_wt,
                               lambda = lambda, maxfitness = maxfitness,
                               alpha = alpha, Q = Q, theta = theta, fun_args = fun_args)
          },
          # FGM coselected mutation 2 environments
          "fgmcsmut2env" = function(lambda, maxfitness, alpha, Q, theta) {
            model_fgmcsmut_2env(nb_mut = nb_mut, genotype_table = genotype_table,
                                fitness_wt_new_env =fitness_wt,
                                lambda = lambda, maxfitness = maxfitness,
                                alpha = alpha, Q = Q, theta = theta, fun_args = fun_args)
          },
          "magellanFix" = function(f) {
            model_Fix(nb_mut = nb_mut, genotype_table = genotype_table, fitness_wt = fitness_wt,
                      f = f, path2_fl_generate = path2_fl_generate)
          },
          "magellanMult" = function(s, S, d) {
            model_Mult(nb_mut = nb_mut, genotype_table = genotype_table, fitness_wt = fitness_wt,
                       s = s, S = S, d = d, path2_fl_generate = path2_fl_generate)
          },
          "magellanHoC" = function(H) {
            model_HoC(nb_mut = nb_mut, genotype_table = genotype_table, fitness_wt = fitness_wt,
                      H = H, path2_fl_generate = path2_fl_generate)
          },
          "magellanNK" = function(K, r) {
            model_NK(nb_mut = nb_mut, genotype_table = genotype_table, fitness_wt = fitness_wt,
                     K = K, r = r, path2_fl_generate = path2_fl_generate)
          },
          "magellanIsing" = function(i, I, c) {
            model_Ising(nb_mut = nb_mut, genotype_table = genotype_table, fitness_wt = fitness_wt,
                        i = i, I = I, c = c, path2_fl_generate = path2_fl_generate)
          },
          "magellanEggBox" = function(e, E) {
            model_EggBox(nb_mut = nb_mut, genotype_table = genotype_table, fitness_wt = fitness_wt,
                         e = e, E = E, path2_fl_generate = path2_fl_generate)
          },
          "magellanOptimum" = function(o, O, p, P) {
            model_Optimum(nb_mut = nb_mut, genotype_table = genotype_table, fitness_wt = fitness_wt,
                          o = o, O = O, p = p, P = P, path2_fl_generate = path2_fl_generate)
          }

  )
}
#' Checkmate wrapper for argument checks. See \code{\link[checkmate]{makeAssertion}}
#' for more informations.
#'
#' @param x The matrix to check.
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in \code{vname}.
#' @param add Collection to store assertion messages.
assert_genotype_table <- function(x, .var.name = checkmate::vname(x), add = NULL) {
  res <- check_genotype_table(x)
  checkmate::makeAssertion(x, res, .var.name, add)
}
#' Test if the matrix is a genotype_table of all the combinations possible of the
#' number of colums in a 0/1 format.
#' @inheritParams assert_genotype_table
#' @return TRUE or an error message.
check_genotype_table <- function(x) {
  res <- length(which( x != 0 & x != 1)) == 0
  if (res) {return(TRUE)} else {
    return(paste0("The rows must only contains combinations of the ", dim(x)[2], " column(s) in a 0/1 format"))
  }
}

#' Returns the fitnesses for each genotype in \code{genotype_table} based on the phenotype
#' of the wild_type (\code{pheno_wt}), the phenotype to fitness function in
#' \code{fitness_fun} and the phenotypic mutation effects in \code{pheno_mut_effect}
#'
#' @param pheno_mut_effect A matrix of real numbers. Phenotypic effects of mutations
#' by rows. The mutations in rows correspond to the columns of genotype_table.
#' @param pheno_wt A vector of real number. Phenotype of the wild type.
#' @param genotype_table A matrix of 0 and 1. The rows correspond to different genotypes
#' and the columns to the mutations that are considered for these genotypes.
#' A genotype (at a certain row) has a given mutation when there is a 1 in the
#' corresponding column. A row with only zeros correspond to the wild type.
#' @inheritParams ptof_fgm_iso
#' @return A vector of real numbers. Each element is the fitness of the corresponding
#' genotype in \code{genotype_table}. They are in the same order as the rows of
#' \code{genotype_table}.
fitness_mutant_genotype <- function(pheno_mut_effect, pheno_wt, genotype_table, maxfitness, alpha = 1/2, Q = 2, pheno_opt = numeric(dim(pheno_mut_effect)[2])) {

  #### check arguments ####
  arg_required <- c("pheno_mut_effect", "pheno_wt", "genotype_table", "maxfitness")
  arg_passed <- names(as.list(match.call())[-1])
  coll1 <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll1$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::reportAssertions(coll1)
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_matrix(genotype_table, mode = "numeric", add = coll)
  assert_genotype_table(genotype_table, add = coll)
  checkmate::assert_matrix(pheno_mut_effect, mode = "numeric",
                           nrow = dim(genotype_table)[2], add = coll)
  checkmate::assert_numeric(pheno_wt, finite = T, any.missing = F, len = dim(pheno_mut_effect)[2],
                            null.ok = F, add = coll)
  checkmate::reportAssertions(coll)

  #### Compute fitness of all the genotypes from genotype_table ####
  pheno_mut_effect_with_wt <- rbind(pheno_wt, pheno_mut_effect) #add the wt in the first row so that for each genotype the true phenotypic coordinates (i.e. phenotypic_effects + wt_coordinates) are returned and not only the sum of the phenotypic effects.
  genotype_table_with_wt <- cbind(array(1, dim = c(dim(genotype_table)[1], 1), dimnames = list(NULL, "WT")), genotype_table) #add a colum with 1 in every row for the "WT" (see also comment above) because every genotype share the WT background except for the mutation sites
  pheno_genotype <- genotype_table_with_wt %*% pheno_mut_effect_with_wt
  ptof_fgm_iso(phenotype = pheno_genotype, maxfitness = maxfitness, alpha = alpha, Q = Q, pheno_opt = pheno_opt)
}

#' Function generating the fitness of the mutants' genotypes in \code{genotype_table}
#' The mutants are produced by combining \code{nb_mut} random mutations generated
#' from the wt with fitness \code{fitness_wt}. Fitnesses are computed using an
#' isotropic FGM with the parameter \code{n}, \code{lambda}, \code{maxfitness},
#' \code{alpha}, \code{Q}, \code{m}. See \code{\link{generate_random_mutation}}
#' for more information on the random mutations and \code{\link{fitness_mutant_genotype}}
#' for the fitness of the genotypes.
#'
#' @param fitness_wt A real number. Fitness of the wild type in the new environment.
#' @inheritParams generate_random_mutation
#' @inheritParams fitness_mutant_genotype
#' @return A vector of real numbers. Each element is the fitness of the corresponding
#' random mutants in \code{genotype_table}. They are in the same order as the rows of
#' \code{genotype_table}.
model_fgmrmut <- function(nb_mut, genotype_table, fitness_wt, n, lambda, maxfitness, alpha, Q, m) {
  pheno_wt <- ftop_fgm_iso(fitness = fitness_wt, n = n, maxfitness = maxfitness,
                           alpha = alpha, Q = Q)
  pheno_mut_effect <- generate_random_mutation(nb_mut = nb_mut,
                                               n = n,
                                               lambda = lambda,
                                               m = m)
  fitness <- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect,
                                     pheno_wt = pheno_wt,
                                     genotype_table = genotype_table,
                                     maxfitness = maxfitness,
                                     alpha = alpha,
                                     Q = Q)
  fitness
}
#' Function generating the fitness of the mutants' genotypes in \code{genotype_table}
#' The mutants are produced by combining \code{nb_mut} selected mutations generated
#' from the wt with fitness \code{fitness_wt}. Fitnesses are computed using an
#' isotropic FGM with the parameter \code{n}, \code{lambda}, \code{maxfitness},
#' \code{alpha}, \code{Q}, \code{m}. See \code{\link{generate_selected_mutation}}
#' for more information on the selected mutations and \code{\link{fitness_mutant_genotype}}
#' for the fitness of the genotypes.
#'
#' @param fitness_wt A real number. Fitness of the wild type in the new environment.
#' @inheritParams generate_selected_mutation
#' @inheritParams fitness_mutant_genotype
#' @param fun_args List with a single element called nb_mut_rand. For more information
#' see \code{\link{generate_selected_mutation}}
#' @return A vector of real numbers. Each element is the fitness of the corresponding
#' selected mutants in \code{genotype_table}. They are in the same order as the rows of
#' \code{genotype_table}.
model_fgmsmut <- function(nb_mut, genotype_table, fitness_wt, n, lambda, maxfitness, alpha, Q, m, fun_args) {
  pheno_wt <- ftop_fgm_iso(fitness = fitness_wt,
                           n = n, maxfitness = maxfitness, alpha = alpha, Q = Q)
  pheno_mut_effect <- do.call(generate_selected_mutation, append(list(nb_mut = nb_mut,
                                                                      n = n,
                                                                      lambda = lambda,
                                                                      maxfitness = maxfitness,
                                                                      pheno_wt = pheno_wt,
                                                                      alpha = alpha,
                                                                      Q = Q,
                                                                      m = m),
                                                                 fun_args))

  if(!anyNA(pheno_mut_effect)) {
    fitness <- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect,
                                       pheno_wt = pheno_wt,
                                       genotype_table = genotype_table,
                                       maxfitness = maxfitness,
                                       alpha = alpha,
                                       Q = Q)
  } else {
    fitness <- array(dim = dim(genotype_table)[1])
  }
  fitness
}
#' Function generating the fitness of the mutants' genotypes in \code{genotype_table}
#' The mutants are produced by combining \code{nb_mut} coselected mutations generated
#' from the wt with fitness \code{fitness_wt}. Fitnesses are computed using an
#' isotropic FGM with the parameter \code{n}, \code{lambda}, \code{maxfitness},
#' \code{alpha}, \code{Q}, \code{m}. See \code{\link{generate_coselected_mutation}}
#' for more information on the coselected mutations and \code{\link{fitness_mutant_genotype}}
#' for the fitness of the genotypes.
#'
#' @param fitness_wt A real number. Fitness of the wild type in the new environment.
#' @inheritParams generate_coselected_mutation
#' @inheritParams fitness_mutant_genotype
#' @param fun_args List with a single element called nb_mut_rand. For more information
#' see \code{\link{generate_coselected_mutation}}
#' @return A vector of real numbers. Each element is the fitness of the corresponding
#' selected mutants in \code{genotype_table}. They are in the same order as the rows of
#' \code{genotype_table}.
model_fgmcsmut <- function(nb_mut, genotype_table, fitness_wt, n, lambda, maxfitness, alpha, Q, m, fun_args) {
  pheno_wt <- ftop_fgm_iso(fitness = fitness_wt,
                           n = n, maxfitness = maxfitness, alpha = alpha, Q = Q)
  pheno_mut_effect <- do.call(generate_coselected_mutation, append(list(nb_mut = nb_mut,
                                                                        n = n,
                                                                        lambda = lambda,
                                                                        maxfitness = maxfitness,
                                                                        pheno_wt = pheno_wt,
                                                                        alpha = alpha,
                                                                        Q = Q,
                                                                        m = m),
                                                                   fun_args))

  if(!anyNA(pheno_mut_effect)) {
    fitness <- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect,
                                       pheno_wt = pheno_wt,
                                       genotype_table = genotype_table,
                                       maxfitness = maxfitness,
                                       alpha = alpha,
                                       Q = Q)
  } else {
    fitness <- array(dim = dim(genotype_table)[1])
  }
  fitness
}
#' Function generating the fitness of the mutants' genotypes in \code{genotype_table}
#' in both the new and the reference environment. Mutants are produced by
#' combining \code{nb_mut} random mutations generated from the wt with fitness
#' \code{fitness_wt_ref} (in \code{fun_args}) in the reference environment.
#' Fitnesses are computed using an isotropic FGM in two dimensions with the parameters
#' passed to \code{fun_args}. See \code{\link{generate_random_mutation}}
#' for more information on the random mutations. The fitness of the mutants in the
#' new environment are then computed in a two (or one) dimensions isotropic FGM
#' with the dimensions corresponding to the two (or the unique) first dimensions
#' in which the mutants were generated in the environment of reference. For more
#' informations on the positionning of the optimum of the new environment see
#' \code{\link{pos_new_env}}. For more informations on the computation of fitness
#' from phenotypes see \code{\link{fitness_mutant_genotype}}.
#'
#' @param fitness_wt_new_env A real number. Fitness of the wild type in the new environment.
#' @inheritParams generate_random_mutation
#' @inheritParams fitness_mutant_genotype
#' @inheritParams pos_new_env
#' @param fun_args List of parameters used to generate random mutants in the reference environment.
#' Must take the form :
#' fun_args = list(fitness_wt_ref = #, n_ref = #, lambda_ref = #, maxfitness_ref = #,
#' alpha_ref = #, Q_ref = #, m_ref = #). For more information on these parameters, see
#' the parameters in \code{\link{model_fgmrmut}}
#' @return A vector of 2 * \code{nb_mut} real numbers. The elements from 1 to \code{nb_mut}
#' are the fitnesses of the corresponding random mutants in \code{genotype_table}
#' (in the same order as the rows) in the reference environment. The elements from
#' \code{nb_mut} + 1 to 2 * \code{nb_mut} are the fitnesses of the same mutants
#' (in the same order) in the new environment.
model_fgmrmut_2env <- function(nb_mut, genotype_table, fitness_wt_new_env, lambda, maxfitness, alpha, Q, theta, fun_args) {

  #### Generate fitness of random mutant in ref environment ####
  pheno_wt_ref <- ftop_fgm_iso(fitness = fun_args$fitness_wt_ref,
                               n = fun_args$n_ref, maxfitness = fun_args$maxfitness_ref, alpha = fun_args$alpha_ref, Q = fun_args$Q_ref)
  pheno_mut_effect_ref <- generate_random_mutation(nb_mut = nb_mut,
                                                   n = fun_args$n_ref,
                                                   lambda = fun_args$lambda_ref,
                                                   m = fun_args$m_ref)
  fitness_ref <- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect_ref,
                                         pheno_wt = pheno_wt_ref,
                                         genotype_table = genotype_table,
                                         maxfitness = fun_args$maxfitness_ref,
                                         alpha = fun_args$alpha_ref,
                                         Q = fun_args$Q_ref)

  #### Compute position of pheno_opt in new environment ####
  pheno_opt_new <- pos_new_env(pheno_wt_ref = pheno_wt_ref, fitness_wt_new_env = fitness_wt_new_env,
                               maxfitness = maxfitness, alpha = alpha, Q = Q, theta = theta)
  #### Compute fitness of random mutant in new environment ####
  lambda_I_n <- lambda * diag(2)
  pheno_mut_effect_new <- sqrt(lambda / fun_args$lambda) * pheno_mut_effect_ref[, 1:2]
  fitness_new_env<- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect_new,
                                            pheno_wt = pheno_wt_ref[, 1:2],
                                            genotype_table = genotype_table,
                                            maxfitness = maxfitness,
                                            alpha = alpha,
                                            Q = Q,
                                            pheno_opt = pheno_opt_new)

  c(fitness_ref, fitness_new_env)
}
#' Function generating the fitness of the mutants' genotypes in \code{genotype_table}
#' in both the new and the reference environment. Mutants are produced by
#' combining \code{nb_mut} selected mutations generated from the wt with fitness
#' \code{fitness_wt_ref} (in \code{fun_args}) in the reference environment.
#' Fitnesses are computed using an isotropic FGM in two dimensions with the parameters
#' passed to \code{fun_args}. See \code{\link{generate_selected_mutation}}
#' for more information on the selected mutations. The fitness of the mutants in the
#' new environment are then computed in a two (or one) dimensions isotropic FGM
#' with the dimensions corresponding to the two (or the unique) first dimensions
#' in which the mutants were generated in the environment of reference. For more
#' informations on the positionning of the optimum of the new environment see
#' \code{\link{pos_new_env}}. For more informations on the computation of fitness
#' from genotype see \code{\link{fitness_mutant_genotype}}.
#'
#' @param fitness_wt_new_env A real number. Fitness of the wild type in the new environment.
#' @inheritParams generate_selected_mutation
#' @inheritParams fitness_mutant_genotype
#' @inheritParams pos_new_env
#' @param fun_args List of parameters used to generate selected mutants in the reference environment.
#' Must take the form :
#' fun_args = list(fitness_wt_ref = #, n_ref = #, lambda_ref = #, maxfitness_ref = #,
#' alpha_ref = #, Q_ref = #, m_ref = #). Can optionnaly take the element nb_mut_ran.
#' For more information on these parameters, see the parameters in \code{\link{model_fgmsmut}}
#' @return A vector of 2 * \code{nb_mut} real numbers. The elements from 1 to \code{nb_mut}
#' are the fitnesses of the corresponding selected mutants in \code{genotype_table}
#' (in the same order as the rows) in the reference environment. The elements from
#' \code{nb_mut} + 1 to 2 * \code{nb_mut} are the fitnesses of the same mutants
#' (in the same order) in the new environment.
model_fgmsmut_2env <- function(nb_mut, genotype_table, fitness_wt_new_env, lambda, maxfitness, alpha, Q, theta, fun_args) {

  #### Generate fitness of random mutant in ref environment ####
  pheno_wt_ref <- ftop_fgm_iso(fitness = fun_args$fitness_wt_ref,
                               n = fun_args$n_ref, maxfitness = fun_args$maxfitness_ref, alpha = fun_args$alpha_ref, Q = fun_args$Q_ref)
  pheno_mut_effect_ref <- do.call(generate_selected_mutation, list(nb_mut = nb_mut,
                                                                   n = fun_args$n_ref,
                                                                   lambda = fun_args$lambda_ref,
                                                                   maxfitness = fun_args$maxfitness_ref,
                                                                   pheno_wt = pheno_wt_ref,
                                                                   alpha = fun_args$alpha_ref,
                                                                   Q = fun_args$Q_ref,
                                                                   m = fun_args$m_ref,
                                                                   nb_mut_rand = if ("nb_mut_rand" %in% names(fun_args)) {fun_args$nb_mut_rand} else {10^4}))

  if(!anyNA(pheno_mut_effect_ref)) {
    fitness_ref <- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect_ref,
                                           pheno_wt = pheno_wt_ref,
                                           genotype_table = genotype_table,
                                           maxfitness = fun_args$maxfitness_ref,
                                           alpha = fun_args$alpha_ref,
                                           Q = fun_args$Q_ref)

    #### Compute position of pheno_opt in new environment ####
    pheno_opt_new <- pos_new_env(pheno_wt_ref = pheno_wt_ref, fitness_wt_new_env = fitness_wt_new_env,
                                 maxfitness = maxfitness, alpha = alpha, Q = Q, theta = theta)
    #### Compute fitness of random mutant in new environment ####
    lambda_I_n <- lambda * diag(2)
    pheno_mut_effect_new <- sqrt(lambda / fun_args$lambda) * pheno_mut_effect_ref[, 1:2]
    fitness_new_env<- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect_new,
                                              pheno_wt = pheno_wt_ref[, 1:2],
                                              genotype_table = genotype_table,
                                              maxfitness = maxfitness,
                                              alpha = alpha,
                                              Q = Q,
                                              pheno_opt = pheno_opt_new)
    all_fitness <- c(fitness_ref, fitness_new_env)
  } else {
    all_fitness <- array(dim = 2 * dim(genotype_table)[1])
  }
  all_fitness
}
#' Function generating the fitness of the mutants' genotypes in \code{genotype_table}
#' in both the new and the reference environment. Mutants are produced by
#' combining \code{nb_mut} coselected mutations generated from the wt with fitness
#' \code{fitness_wt_ref} (in \code{fun_args}) in the reference environment.
#' Fitnesses are computed using an isotropic FGM in two dimensions with the parameters
#' passed to \code{fun_args}. See \code{\link{generate_coselected_mutation}}
#' for more information on the coselected mutations. The fitness of the mutants in the
#' new environment are then computed in a two (or one) dimensions isotropic FGM
#' with the dimensions corresponding to the two (or the unique) first dimensions
#' in which the mutants were generated in the environment of reference. For more
#' informations on the positionning of the optimum of the new environment see
#' \code{\link{pos_new_env}}. For more informations on the computation of fitness
#' from genotype see \code{\link{fitness_mutant_genotype}}.
#'
#' @param fitness_wt_new_env A real number. Fitness of the wild type in the new environment.
#' @inheritParams generate_selected_mutation
#' @inheritParams fitness_mutant_genotype
#' @inheritParams pos_new_env
#' @param fun_args List of parameters used to generate coselected mutants in the reference environment.
#' Must take the form :
#' fun_args = list(fitness_wt_ref = #, n_ref = #, lambda_ref = #, maxfitness_ref = #,
#' alpha_ref = #, Q_ref = #, m_ref = #). Can optionnaly take the element nb_mut_ran.
#' For more information on these parameters, see the parameters in \code{\link{model_fgmcsmut}}
#' @return A vector of 2 * \code{nb_mut} real numbers. The elements from 1 to \code{nb_mut}
#' are the fitnesses of the corresponding coselected mutants in \code{genotype_table}
#' (in the same order as the rows) in the reference environment. The elements from
#' \code{nb_mut} + 1 to 2 * \code{nb_mut} are the fitnesses of the same mutants
#' (in the same order) in the new environment.
model_fgmcsmut_2env <- function(nb_mut, genotype_table, fitness_wt_new_env, lambda, maxfitness, alpha, Q, theta, fun_args) {

  #### Generate fitness of random mutant in ref environment ####
  pheno_wt_ref <- ftop_fgm_iso(fitness = fun_args$fitness_wt_ref,
                               n = fun_args$n_ref, maxfitness = fun_args$maxfitness_ref, alpha = fun_args$alpha_ref, Q = fun_args$Q_ref)
  pheno_mut_effect_ref <- do.call(generate_coselected_mutation, list(nb_mut = nb_mut,
                                                                     n = fun_args$n_ref,
                                                                     lambda = fun_args$lambda_ref,
                                                                     maxfitness = fun_args$maxfitness_ref,
                                                                     pheno_wt = pheno_wt_ref,
                                                                     alpha = fun_args$alpha_ref,
                                                                     Q = fun_args$Q_ref,
                                                                     m = fun_args$m_ref,
                                                                     nb_mut_rand = if ("nb_mut_rand" %in% names(fun_args)) {fun_args$nb_mut_rand} else {10^4}))
  if(!anyNA(pheno_mut_effect_ref)) {
    fitness_ref <- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect_ref,
                                           pheno_wt = pheno_wt_ref,
                                           genotype_table = genotype_table,
                                           maxfitness = fun_args$maxfitness_ref,
                                           alpha = fun_args$alpha_ref,
                                           Q = fun_args$Q_ref)

    #### Compute position of pheno_opt in new environment ####
    pheno_opt_new <- pos_new_env(pheno_wt_ref = pheno_wt_ref, fitness_wt_new_env = fitness_wt_new_env,
                                 maxfitness = maxfitness, alpha = alpha, Q = Q, theta = theta)
    #### Compute fitness of random mutant in new environment ####
    lambda_I_n <- lambda * diag(2)
    pheno_mut_effect_new <- sqrt(lambda / fun_args$lambda) * pheno_mut_effect_ref[, 1:2]
    fitness_new_env<- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect_new,
                                              pheno_wt = pheno_wt_ref[, 1:2],
                                              genotype_table = genotype_table,
                                              maxfitness = maxfitness,
                                              alpha = alpha,
                                              Q = Q,
                                              pheno_opt = pheno_opt_new)
    all_fitness <- c(fitness_ref, fitness_new_env)
  } else {
    all_fitness <- array(dim = 2 * dim(genotype_table)[1])
  }
  all_fitness
}
#' Returns the coordinates of a new optimum in the plane corresponding to the
#' first two dimensions of \code{pheno_wt_ref}. The position is computed using the
#' inverse of the fitness function of an isotrope FGM on the \code{fitness_wt_new_env}.
#'
#' @param pheno_wt_ref A vector of real number. Phenotype of the wild type.
#' @param fitness_wt_new_env A real number. Fitness of the wild type in the new
#' env.
#' @inheritParams ptof_fgm_iso
#' @param theta A angle in radian between 0 and 2 * pi. Angle between the optimum
#' of the environment of reference (first dimension of \code{pheno_wt_ref}), the
#' \code{pheno_wt_ref} and the optimum of the new environement.
#' @return A pair of coordinates for the position of the optimum of the new environment
#' corresponding to the two first dimensions of \code{pheno_wt_ref}. If \code{pheno_wt_ref}
#' has a single dimension, \code{theta} is ignored and the position of the
#' optimum is given in this single dimension.
pos_new_env <- function(pheno_wt_ref, fitness_wt_new_env, maxfitness, alpha = 1/2, Q = 2, theta) {

  #### check arguments ####
  arg_required <- c("pheno_wt_ref", "fitness_wt_new_env", "maxfitness", "alpha", "Q", "theta")
  arg_passed <- names(as.list(match.call())[-1])
  coll1 <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll1$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::reportAssertions(coll1)
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(pheno_wt_ref, finite = T, any.missing = F,
                            null.ok = F, add = coll)
  checkmate::assert_number(fitness_wt_new_env, na.ok = F, upper = maxfitness,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_number(maxfitness, na.ok = F, finite = T, null.ok = F,
                           add = coll)
  checkmate::assert_number(alpha, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_number(Q, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_number(theta, na.ok = F, lower = 0, upper = 360*pi/180,
                           finite = T, null.ok = F, add = coll)
  checkmate::reportAssertions(coll)

  #### Compute position of pheno_opt in new environment ####
  dist_to_new_opt <- (1 / alpha * (maxfitness - fitness_wt_new_env))^(1/Q)
  pheno_opt_new <- numeric(2)
  if (length(pheno_wt_ref) > 1) {
    pheno_opt_new[1] <- pheno_wt_ref[1] + dist_to_new_opt * cos(theta) # coordinate of new env in dim 1
    pheno_opt_new[2] <- pheno_wt_ref[2] + dist_to_new_opt * sin(theta) # coordinate of new env in dim 2
  } else {
    pheno_opt_new <- pheno_wt_ref + dist_to_new_opt
  }
  pheno_opt_new
}
