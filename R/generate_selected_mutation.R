#' Generates \code{nb_rand_mut} mutants and select \code{nb_mut} among them according
#' to their selection coefficient and returns their phenotypic effects.
#'
#' @inheritParams generate_random_mutation
#' @inheritParams ptof_fgm_iso
#' @param pheno_wt A vector of real numbers. Phenotypic coordinates of the wild-type.
#' Its length must be equal to the number of columns of \code{pheno_mut_effect}.
#' @param nb_mut_rand A natural number superior or equal to \code{nb_mut}. Number
#' of random mutations to draw from which the mutant will be selected. This number
#' must be large compared to \code{nb_mut}. Default=10^4.
#' @return A matrix with \code{nb_mut} rows of selected mutation effects in \code{n}
#' columns (phenotypic dimension(s)). Note that if the number of drawn beneficial mutants
#' is inferior to \code{nb_mut} the matrix is filed with NAs. In this case \code{nb_mut_rand}
#' must be increased.
#' @examples
#' generate_selected_mutation(nb_mut = 5, n = 3, lambda = 0.1, maxfitness = 1, pheno_wt = c(-1, 0, 0))
#' generate_selected_mutation(nb_mut = 5, n = 3, lambda = 0.1, maxfitness = 1, pheno_wt = c(-1, 0, 0),
#'                            alpha = 1, Q = 0.5, m = 2, nb_mut_rand = 10^5)
#' @export
generate_selected_mutation <- function(nb_mut, n, lambda, maxfitness, pheno_wt, alpha = 1/2, Q = 2, m = n, nb_mut_rand = 10^4) {

  #### check arguments ####
  arg_required <- c("nb_mut", "n", "lambda", "maxfitness", "pheno_wt")
  arg_passed <- names(as.list(match.call())[-1])
  coll1 <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll1$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::reportAssertions(coll1)
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_count(nb_mut, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::assert_numeric(pheno_wt, finite = T, any.missing = F, len = n,
                            null.ok = F, add = coll)
  checkmate::assert_count(nb_mut_rand, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::reportAssertions(coll)


  #### Generates random mutants ####
  pheno_rand_mut_effect <- generate_random_mutation(nb_mut = nb_mut_rand, n = n, lambda = lambda, m = m)

  #### Computes mutants and wt fitnesses ####
  fitness_wt <- ptof_fgm_iso(phenotype = matrix(pheno_wt, nrow = 1, ncol = n),
                             maxfitness = maxfitness, alpha = alpha, Q = Q)
  pheno_mutant <- add_mut_to_pheno(pheno_mut_effect = pheno_rand_mut_effect, pheno_wt = pheno_wt)
  fitness_mutant <- ptof_fgm_iso(phenotype = pheno_mutant, maxfitness = maxfitness, alpha = alpha, Q = Q)

  #### Select nb_mut mutants according to their selection coefficient ####
  idx_mutant_benef <- which(fitness_mutant > fitness_wt)
  if(length(idx_mutant_benef) >= nb_mut){ #select the mutants among the beneficials
    idx_mutant_benef_selected <- idx_mutant_benef
    if(length(idx_mutant_benef) > 1) {
      idx_mutant_benef_selected <- sample(idx_mutant_benef,
                                          nb_mut,
                                          prob = prob_esc_drift(fitness_mutant[idx_mutant_benef] - fitness_wt))
    }
    pheno_select_mut_effect <- pheno_rand_mut_effect[idx_mutant_benef_selected, , drop=FALSE]
  } else { #produce NA if length(idx_mutant_benef) < nb_sel
    pheno_select_mut_effect <- array(dim = c(nb_mut, dim(pheno_rand_mut_effect)[2]))
  }
  pheno_select_mut_effect
}
#' Add phenotypic effects of mutation(s) from \code{pheno_mut_effect} to a given
#' phenotype \code{pheno_wt}.
#' @param pheno_mut_effect A matrix of real numbers. Each row correspond
#' to the phenotypic effects of given mutation in each dimensions (columns).
#' @param pheno_wt A vector of real numbers. Phenotypic coordinates of the wild-type.
#' Its length must be equal to the number of columns of \code{pheno_mut_effect}.
#' @return A matrix of the same dimensions as \code{pheno_rand_mut_effect} with
#' each row equal to the sum of \code{pheno_wt} with a row of \code{pheno_rand_mut_effect}.
add_mut_to_pheno <- function(pheno_mut_effect, pheno_wt){
  #### Add the pheno_wt to each row of pheno_mut_effect ####
  pheno_mut_effect + matrix(pheno_wt, nrow=dim(pheno_mut_effect)[1], ncol=dim(pheno_mut_effect)[2], byrow = TRUE)
}
#' Computes the probability of escaping drift for a mutation in a clonal population,
#' using its selection coefficient. For more info see Gerrish et al.(1998) Appendix I
#' @param sel_coeff A vector of real numbers corresponfing to selectioncoefficients.
#' @return vetctor of probabilities.
prob_esc_drift <- function(sel_coeff) {4 * sel_coeff / (1 + sel_coeff)^2}
