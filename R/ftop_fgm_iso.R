#' Generates the phenotypic coordinates, in a \code{n}-dimension(s) phenotypic space,
#' corresponding to a given fitness using Fisher's Geometric Model (FGM, Fisher (1930)).
#' Use the inverse of the fitness function of an isotrope FGM (with parameters
#' \code{n}, \code{maxfitness}, \code{alpha}, \code{Q}) to compute the euclidian
#' distance to the phenotypic optimum (\code{pheno_opt} -at the origin by default-
#' at which the fitness = \code{maxfitness}). This distance is reported on the
#' first coordinate and all the other \code{n}-1 coordinates are equal to the
#' coordinates of the phenotypic optimum.
#' See also \code{\link{ptof_fgm_iso}} for inverse function.
#'
#' @param fitness A real number. The fitness of the phenotype. Must be lower or
#' equal to \code{maxfitness}.
#' @param n A natural number. Number of dimensions of the phenotypic space.
#' @param maxfitness A real number. The maximum fitness in the landscape. The
#' fitness at the phenotypic optimum (\code{pheno_opt}).
#' @param alpha A strictly positive real number. Scaling factor for the fitness
#' function. Default=1/2 in the cannonical FGM with a quadratic fitness function.
#' @param Q A strictly positive number. "Shape" of the fitness function. Default=2
#' in the cannonical FGM with a quadratic fitness function.
#' @param pheno_opt A vector of coordinates for the position of the phenotypic
#' optimum at which the fitness is equal to \code{maxfitness}. Its length must be
#' equal to \code{n}. Default = numeric(n).
#' @return A vector of \code{n} coordinates in the phenotypic space.
#' @examples
#' ftop_fgm_iso(fitness = 0, n = 3, maxfitness = 1)
#' ftop_fgm_iso(fitness = 0, n = 3, maxfitness = 1, alpha = 1/2, Q = 2, pheno_opt = c(1,1,1))
#' @export
ftop_fgm_iso <- function(fitness, n, maxfitness, alpha = 1/2, Q = 2, pheno_opt = numeric(n)){

  #### check arguments ####
  arg_required <- c("fitness", "n", "maxfitness")
  arg_passed <- names(as.list(match.call())[-1])
  coll1 <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll1$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::reportAssertions(coll1)
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_number(maxfitness, na.ok = F, finite = T, null.ok = F,
                           add = coll)
  checkmate::assert_numeric(fitness, upper = maxfitness, finite = T,
                            any.missing = F, null.ok = F, add = coll)
  checkmate::assert_count(n, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::assert_number(alpha, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_number(Q, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_numeric(pheno_opt, finite = T, any.missing = F, len = n,
                            null.ok = F, add = coll)
  checkmate::reportAssertions(coll)

  #### phenotype from fitness ####
  phenotype <- vapply(X = fitness,
                      FUN = function(f) {pheno_opt - c((1 / alpha * (maxfitness-f))^(1/Q), numeric(n-1))},
                      FUN.VALUE = numeric(n))
  if (n > 1) phenotype <- t(phenotype)
  phenotype
}
