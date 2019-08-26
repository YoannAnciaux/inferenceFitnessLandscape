#' Generates the fitness corresponding to the phenotypic coordinates of \code{phenotype}
#' using Fisher's Geometric Model (FGM, Fisher (1930)). Use the fitness function
#' of an isotrope FGM (with parameters \code{n},  \code{maxfitness}, \code{alpha},
#' \code{Q}) to compute the fitness (malthusian fitness) of phenotype(s)
#' \eqn{log(W) = \code(maxfitness)- \code{alpha} * \code{phenotype}^(\code{Q}/2)}
#' (see Tenaillon et al. (2007))
#' See also \code{\link{ftop_fgm_iso}} for inverse function.
#'
#' @param phenotype A vector, matrix of real number(s) (phenotypic
#' coordinate(s)). For a matrix or a data.frame, the rows are phenotypes and the
#' columns phenotypic dimensions.
#' @param maxfitness A real number. The maximum fitness in the landscape. The
#' fitness at the phenotypic optimum (\code{pheno_opt}).
#' @param alpha A strictly positive real number. Scaling factor for the fitness
#' function. Default=1/2 in the cannonical FGM with a quadratic fitness function.
#' @param Q A strictly positive number. "Shape" of the fitness function. Default=2
#' in the cannonical FGM with a quadratic fitness function.
#' @param pheno_opt A vector of coordinates for the position of the phenotypic
#' optimum at which the fitness is equal to \code{maxfitness}. Its length must be
#' equal to the the number of column of phenotype.
#' @return A vector of fitnesses of length equal to the number of phenotype(s) (row(s)) in \code{phenotype}
#' @examples
#' #' @examples
#' ptof_fgm_iso(phenotype = matrix(1:9, 3, 3, byrow = TRUE), maxfitness = 1)
#' ptof_fgm_iso(phenotype = matrix(1:9, 3, 3, byrow = TRUE), maxfitness = 1,
#' alpha = 1/2, Q = 2, pheno_opt = c(1,1,1))
#' @export
ptof_fgm_iso <- function(phenotype, maxfitness, alpha = 1/2, Q = 2, pheno_opt = numeric(dim(phenotype)[2])) {

  #### check arguments ####
  arg_required <- c("phenotype", "maxfitness")
  arg_passed <- names(as.list(match.call())[-1])
  coll <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::assert_matrix(phenotype, mode = "numeric", any.missing = F, null.ok = F,
                           add = coll)
  checkmate::assert_number(maxfitness, na.ok = F, finite = T, null.ok = F,
                           add = coll)
  checkmate::assert_number(alpha, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_number(Q, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::reportAssertions(coll)


  #### fitness from phenotype ####
  dist_pheno <- as.matrix(apply(X = phenotype,
                                MARGIN = 1,
                                FUN = function(p) {pheno_opt - p}))
  if (dim(phenotype)[2] > 1) dist_pheno <- t(dist_pheno)
  maxfitness - alpha * (rowSums(dist_pheno^2))^(Q/2)
}
