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
#'   \item{"fgmrmut"}{random mutation in isotropic FGM. See \code{\link{generate_random_mutation}}}
#'   \item{"fgmsmut"}{selected mutation in isotropic FGM. See \code{\link{generate_selected_mutation}}}
#'   \item{"fgmcsmut"}{coselected mutation in isotropic FGM. See \code{\link{generate_coselected_mutation}}}
#' }
#' @param fun_args List of optionnal argument for certain \code{model_type}:
#' \describe{
#'   \item{nb_mut_rand}{for "fgmsmut" and "fgmcsmut". See \code{\link{generate_selected_mutation}} and \code{\link{generate_coselected_mutation}}}
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
  coll <- checkmate::makeAssertCollection()
  if (is.character(empirical_fl)) {
    checkmate::assert_file_exists(empirical_fl, access = "r", add = coll)
    if (checkmate::test_file_exists(empirical_fl, access = "r")) {
      empirical_fl <- utils::read.table(file = empirical_fl, ...)
    }
  }
  checkmate::assert_matrix(empirical_fl, mode = "numeric", add = coll)
  if (checkmate::test_matrix(empirical_fl, mode = "numeric")) {
    assert_genotype_table(empirical_fl[, -dim(empirical_fl)[2]],
                          .var.name = paste0(checkmate::vname(empirical_fl), " without the last column (fitness)"),
                          add = coll)
  }
  checkmate::assert_choice(model_type, choices = c("fgmrmut", "fgmsmut", "fgmcsmut"),
                           add = coll)
  #when adding new models check all the values from the two following assertions
  checkmate::assert_list(fun_args, type = character(0), any.missing = F,
                         max.len = 1, names = "unique", null.ok = F, add = coll)
  checkmate::assert_subset(names(fun_args), choices = c("nb_mut_rand"), add = coll)
  checkmate::reportAssertions(coll)
  nb_mut <- dim(empirical_fl)[2]-1
  genotype_table <- empirical_fl[, -(nb_mut + 1)]

  switch (model_type,
          # FGM random mutation 1 environment
          "fgmrmut" = function(n, lambda, maxfitness, alpha, Q, m) {
            pheno_wt <- ftop_fgm_iso(fitness = empirical_fl[which(rowSums(genotype_table) == 0), nb_mut + 1],
                                     n = n, maxfitness = maxfitness, alpha = alpha, Q = Q)
            pheno_mut_effect <- generate_random_mutation(nb_mut = nb_mut,
                                                         n = n,
                                                         lambda = lambda,
                                                         m = m)
            fitness <- fitness_mutant_genotype(pheno_mut_effect = pheno_mut_effect,
                                               pheno_wt = pheno_wt,
                                               geno_table = genotype_table,
                                               maxfitness = maxfitness,
                                               alpha = alpha,
                                               Q = Q)
            fitness
          },
          # FGM selected mutation 1 environment
          "fgmsmut" = function(n, lambda, maxfitness, alpha, Q, m) {
            pheno_wt <- ftop_fgm_iso(fitness = empirical_fl[which(rowSums(genotype_table) == 0), nb_mut + 1],
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
                                                 geno_table = genotype_table,
                                                 maxfitness = maxfitness,
                                                 alpha = alpha,
                                                 Q = Q)
            } else {
              fitness <- array(dim = dim(empirical_fl)[1])
            }
            fitness
          },
          # FGM coselected mutation 1 environment
          "fgmcsmut" = function(n, lambda, maxfitness, alpha, Q, m) {
            pheno_wt <- ftop_fgm_iso(fitness = empirical_fl[which(rowSums(genotype_table) == 0), nb_mut + 1],
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
                                                 geno_table = genotype_table,
                                                 maxfitness = maxfitness,
                                                 alpha = alpha,
                                                 Q = Q)
            } else {
              fitness <- array(dim = dim(empirical_fl)[1])
            }
            fitness
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
  x2 <- as.matrix(expand.grid(rep(list(0:1), dim(x)[2])))
  res <- identical(as.logical(duplicated(rbind(x, x2))), c(logical(dim(x)[1]), !logical(dim(x)[1])))
  if (res) {return(TRUE)} else {
    return(paste0("The rows must contains the ", 2^dim(x)[2], " combinations of the ", dim(x)[2], " column(s) in a 0/1 format"))
  }
}

#' Returns the fitnesses for each genotype in \code{geno_table} based on the phenotype
#' of the wild_type (\code{pheno_wt}), the phenotype to fitness function in
#' \code{fitness_fun} and the phenotypic mutation effects in \code{pheno_mut_effect}
#'
#' @param pheno_mut_effect A matrix of real numbers. Phenotypic effects of mutations
#' by rows. The mutations in rows correspond to the columns of geno_table.
#' @param pheno_wt A vector of real number. Phenotype of the wild type.
#' @param geno_table A matrix of 0 and 1. The rows correspond to different genotypes
#' and the columns to the mutations that are considered for these genotypes.
#' A genotype (at a certain row) has a given mutation when there is a 1 in the
#' corresponding column. A row with only zeros correspond to the wild type.
#' @inheritParams ptof_fgm_iso
#' @return A vector of real numbers. Each element is the fitness of the corresponding
#' genotype in \code{geno_table}. They are in the same order as the rows of
#' \code{geno_table}.
#' @examples
#' geno_table <- as.matrix(expand.grid(rep(list(0:1), 3)))
#' pheno_wt <- c(-1, 0)
#' pheno_mut_effect <- generate_random_mutation(nb_mut = 3, n = 2, lambda = 0.1)
#' fitness_mutant_genotype(pheno_mut_effect, pheno_wt, geno_table, maxfitness = 1)
#' @export
fitness_mutant_genotype <- function(pheno_mut_effect, pheno_wt, geno_table, maxfitness, alpha = 1/2, Q = 2) {

  #### check arguments ####
  arg_required <- c("pheno_mut_effect", "pheno_wt", "geno_table", "maxfitness")
  arg_passed <- names(as.list(match.call())[-1])
  coll <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::assert_matrix(geno_table, mode = "numeric", add = coll)
  assert_genotype_table(geno_table, add = coll)
  checkmate::assert_matrix(pheno_mut_effect, mode = "numeric",
                           nrow = dim(geno_table)[2], add = coll)
  checkmate::assert_numeric(pheno_wt, finite = T, any.missing = F, len = dim(pheno_mut_effect)[2],
                            null.ok = F, add = coll)

  #### Compute fitness of all the genotypes from geno_table ####
  pheno_mut_effect_with_wt <- rbind(pheno_wt, pheno_mut_effect) #add the wt in the first row so that for each genotype the true phenotypic coordinates (i.e. phenotypic_effects + wt_coordinates) are returned and not only the sum of the phenotypic effects.
  geno_table_with_wt <- cbind(array(1, dim = c(dim(geno_table)[1], 1), dimnames = list(NULL, "WT")), geno_table) #add a colum with 1 in every row for the "WT" (see also comment above) because every genotype share the WT background except for the mutation sites
  pheno_genotype <- geno_table_with_wt %*% pheno_mut_effect_with_wt
  ptof_fgm_iso(phenotype = pheno_genotype, maxfitness = maxfitness, alpha = alpha, Q = Q)
}
