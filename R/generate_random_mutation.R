#' Draw \code{nb_mut} mutation(s) effect(s) on phenotype of \code{n} dimensions
#' with a variance of \code{lambda} per phenotypic dimension.
#' Also applies any restricted pleiotropy of level \code{m}
#'
#' @param nb_mut A natural number. Number of mutations in the output.
#' @param n A natural number. Number of dimensions of the phenotypic space in
#' which the random mutations are drawn.
#' @param lambda A positive real number. Variance of the effect of mutations on
#' phenotype per phenotypic dimension.
#' @param m A natural number inferior to \code{n}. Level of restricted pleiotropy
#' which corresponds to the number of dimensions for which each mutation have a
#' non-zero phenotypic effect. i.e. effects are equal to 0 in n-m dimensions.
#' Default=\code{n} which corresponds to full pleiotropy.
#' @return A matrix with \code{nb_mut} rows of random mutation effects in \code{n}
#' columns (phenotypic dimension(s)).
#' @examples
#' #full pleiotropy
#' generate_random_mutation(nb_mut = 3, n = 3, lambda = 0.1)
#' #restricted pleiotropy
#' generate_random_mutation(nb_mut = 3, n = 3, lambda = 0.1, m = 1)
#' @export
generate_random_mutation <- function(nb_mut, n, lambda, m = n) {

  #### check arguments ####
  arg_required <- c("nb_mut", "n", "lambda")
  arg_passed <- names(as.list(match.call())[-1])
  coll1 <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll1$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::reportAssertions(coll1)
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_count(nb_mut, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::assert_count(n, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::assert_number(lambda, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_int(m, na.ok = F, lower = 1, upper = n, null.ok = F,
                        add = coll)
  checkmate::reportAssertions(coll)

  #### Draw random mutation effects ####
  pheno_rand_mut_effect <- matrix(random_mutation(n = n, nb_mut = nb_mut, lambda = lambda),
                                  nrow = nb_mut, ncol = n, byrow = T)
  if(m < n) {
    pheno_rand_mut_effect <- restricted_pleiotropy(mutation = pheno_rand_mut_effect, m = m)
  }
  pheno_rand_mut_effect
}
#' Draw \code{nb_mut} mutation(s) effect(s) on phenotype of \code{n} dimensions
#' @inheritParams generate_random_mutation
#' @return A matrix of \code{nb_mut} rows and \code{n} columns.
random_mutation <- function(n, nb_mut, lambda) {
  #### Draw random mutation effects ####
  lambda_I_n <- lambda * diag(n) #Covariance matrix. Here no covariances and same variance in every dimensions.
  MASS::mvrnorm(n = nb_mut, mu = numeric(n), Sigma = lambda_I_n) #faster than drawing random numbers out of the program and sampling in these number in the program
}
#' Creates restricted pleiotropy (as defined in Chevin et al.(2010)) of a level
#' \code{m} to a matrix of phenotypic effects of mutations. Insert a 0 at \code{m}
#' random positions for each row of the matrix.
#' @param mutation A matrix of real number. The rows are mutations and the columns
#' phenotypic dimensions.
#' @inheritParams generate_random_mutation
#' @return A matrix of the same number of rows and columns as \code{mut} with
#' \code{m} random positions equal to zero in each rows
restricted_pleiotropy <- function(mutation, m) {
  #### Apply restricted pleiotriopy ####
  n <- dim(mutation)[2]
  t(apply(X = mutation, MARGIN = 1, FUN = function(mut) {mut[sample(c(1:n), n-m)] <- 0.0; mut}))
}
