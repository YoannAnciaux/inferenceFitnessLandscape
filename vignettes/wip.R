#' Computes the fitness (malthusian fitness) of phenotype(s) using Fisher's Geometric Model (FGM)
#' (Fisher (1930)).
#' The fitness is log(W) = \code(maxfitness)- \code{alpha} * \code{phenotype}^(Q/2)
#' (see Tenaillon et al. (2007))
#'
#' @param phenotype A vector, matrix of real number(s) (phenotypic
#' coordinate(s)). For a matrix or a data.frame, the rows are phenotypes and the
#' columns phenotypic dimensions.
#' @param maxfitness A real number. Maximum fitness in the landscape reached at the
#' phenotypic optimum.
#' @param alpha A strictly positive number. Scaling factor for the fitness
#' function. Default=1/2 in the cannonical FGM with a quadratic fitness function.
#' @param Q A strictly positive number. "Flatness" of the fitness
#' function. Default=2 in the cannonical FGM with a quadratic fitness function.
#' @return A vector of fitnesses of length equal to the number of phenotypes in \code{phenotype}
#' @examples
#' ## Fitness of 3 single mutants
#' phenotypic_effect_mutations <- random_mutation(n = 3, nb_mut = 5, sigma = 0.1)
#' wild_type_pheno <- rep(0, 3)
#' fgm_iso(phenotype = wild_type_pheno + phenotypic_effect_mutations, maxfitness = 1)
ptof_fgm_iso <- function(phenotype, maxfitness, alpha = 1/2, Q = 2) {

  #checks : phenotype, maxfitness, alpha, Q
  if (missing(phenotype)){
    stop("'phenotype' must be supplied.", call.=FALSE)
  }
  stopifnot(is.numeric(phenotype), length(dim(phenotype)) <=2)
  if (missing(maxfitness)){
    stop("'maxfitness' must be supplied.", call.=FALSE)
  }
  stopifnot(is.numeric(maxfitness), length(maxfitness) == 1)
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0)
  stopifnot(is.numeric(Q), length(Q) == 1, Q > 0)

  UseMethod("fgm_iso")
}

fgm_iso.matrix <- function(phenotype, maxfitness, alpha = 1/2, Q = 2) {
  maxfitness - alpha * (rowSums(phenotype^2))^(Q/2)
}

fgm_iso.default <- function(phenotype, maxfitness, alpha = 1/2, Q = 2) {
  maxfitness - alpha * (sum(phenotype^2))^(Q/2)
}
