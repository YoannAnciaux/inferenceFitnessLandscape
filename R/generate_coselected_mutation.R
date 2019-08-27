#' Generates \code{nb_mut} co_selected mutations by calling sequentially
#' \code{\link{generate_selected_mutation}}. So that the second mutation happens
#' on the genetic background of the first selected mutants, the third mutations
#' on the genetic background of the second selected mutants...
#'
#' @inheritParams generate_selected_mutation
#' @return A matrix of \code{nb_mut} rows of coselected mutations effects, ordered
#' by rows from the first mutation to the last mutation in increasing order. Note
#' that if no beneficial mutation is drawn for a given mutant, NA are introduced
#' from the corresponding row in the phenotype matrix to the end of the matrix.
#' In this case \code{nb_mut_rand} must be increased.
#' @examples
#' generate_coselected_mutation(nb_mut = 5, n = 3, lambda = 0.1, maxfitness = 1,
#'                              pheno_wt = c(-1, 0, 0))
#' generate_coselected_mutation(nb_mut = 5, n = 3, lambda = 0.1, maxfitness = 1,
#'                              pheno_wt = c(-1, 0, 0), alpha = 1, Q = 0.5,
#'                              m = 2, nb_mut_rand = 10^5)
#' @export
generate_coselected_mutation <- function(nb_mut, n, lambda, maxfitness, pheno_wt, alpha = 1/2, Q = 2, m = n, nb_mut_rand = 10^4) {

  #### check arguments ####
  arg_required <- c("nb_mut", "n", "lambda", "maxfitness", "pheno_wt")
  arg_passed <- names(as.list(match.call())[-1])
  coll <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::assert_count(nb_mut, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::assert_count(n, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::assert_numeric(pheno_wt, finite = T, any.missing = F, len = n,
                            null.ok = F, add = coll)
  checkmate::reportAssertions(coll)

  pheno_coselect_mut_effect <- array(NA, dim = c(nb_mut, n))
  for(cosel in 1:nb_mut){
    pheno_select_mut_effect <- generate_selected_mutation(nb_mut = 1, n, lambda, maxfitness, pheno_wt,
                                                          alpha, Q, m, nb_mut_rand)
    if(anyNA(pheno_select_mut_effect)) break
    pheno_coselect_mut_effect[cosel, ] <- pheno_select_mut_effect
    pheno_wt <- pheno_select_mut_effect + pheno_wt
  }
  pheno_coselect_mut_effect
}
