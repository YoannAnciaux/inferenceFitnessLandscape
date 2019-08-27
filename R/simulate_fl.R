#' Apply \code{simulation_model} to the rows of the matrix \code{parameter} either
#' sequentially or in parallel (depends on \code{ncore}).
#'
#' @param parameter A matrix of type numeric. Each row correspond to the arguments
#'  in \code{simulation_model}. the order of the column must be the same as the arguments
#'  in \code{simulation_model}.
#' @param simulation_model A function or a function name. The function's arguments
#' must cotrrespond to a row of \code{parameter} with same order.
#' @param ncore A natural number. Number of cores used for the simulations. If
#' \code{ncore} = 1 simulation_model is applied sequentially over the row of
#' \code{parameter}. If \code{ncore} > 1 simulation_model is applied in parallel
#' over the row of \code{parameter}. The parallelization is not avilable for "Windows".
#' @param file_output Full path to the output file for the simulations. If NULL
#' no output file is created.
#' @param multi_file Logical. Applies only if \code{file_output} is not NULL. If
#' TRUE the results of the simulations are saved in two separated files with the
#' suffixes "-simulation" and "-parameter". If FALSE parameters and simulations
#' are saved in the same file with first columns for the parameters and the last
#' for the corresponding simulations.
#' @param ... Extra arguments which will be passed to \code{write.table} when an
#' output is saved in a file. e.g. \code{sep}
#' @return A list of three elements. The first element is the matrix of \code{parameter}.
#' If some parameters are unnamed, all the parameters are renamed to P1, P2...
#' The second element is the matrix of simulations, in which each row is the result
#' of \code{simulation_model} applied on the corresponding row in \code{parameter}. The
#' columns of this matrix are named S1, S2...
#' The third element is the output of \code{system.time} for the total simulation
#' process.
#' @examples
#' simulation_model <- function(m, sd) {rnorm(n = 10, mean = m, sd = sd)}
#' parameter <- cbind(1:10, seq(0.1, 1, 0.1))
#' simulate_fl(parameter, simulation_model, ncore = 1)
#' @export
simulate_fl <- function(parameter, simulation_model, ncore = 1, file_output = NULL, multi_file = TRUE, ...) {

  #### check arguments ####
  arg_required <- c("parameter", "simulation_model")
  arg_passed <- names(as.list(match.call())[-1])
  coll <- checkmate::makeAssertCollection()
  if (!checkmate::test_subset(x = arg_required,
                              choices = names(as.list(match.call())[-1]))) {
    coll$push(paste0("Missing values for ", paste(setdiff(arg_required, arg_passed), collapse=", ")))
  }
  checkmate::assert_matrix(parameter, any.missing = T, all.missing = T, min.rows = 1,
                           add = coll)
  checkmate::assert_function(simulation_model, nargs = dim(parameter)[2], null.ok = F,
                             add = coll)
  if (!checkmate::test_function(simulation_model, nargs = dim(parameter)[2], null.ok = F)){
    warning("See the vignette for more information on requirements for the simulation_model function")
  }
  checkmate::assert_count(ncore, na.ok = F, positive = T, null.ok = F, add = coll)
  if (ncore > 1) {
    checkmate::assert_os(os = c("linux", "mac", "solaris"), add = coll)
    if (!checkmate::test_os(os = c("linux", "mac", "solaris"))) {
      warning("The parallel computation uses 'forking' which is only available for linux, mac and solaris")
    }
  }
  if (!is.null(file_output)) {
    checkmate::assert_path_for_output(file_output, overwrite = TRUE, extension = "csv")
    checkmate::assert_logical(multi_file, any.missing = F, len = 1, null.ok = F,
                              add = coll)
  }
  checkmate::reportAssertions(coll)

  #### Simulation either in parallel (only for linux, mac and solaris) or non-parallel ####
    # be carefull parallel forking only work for linux/macOS for windows must use Psock but need to export the environement
    timing <- system.time(
      simulation <- matrix(unlist(parallel::mclapply(1:nrow(parameter),
                                                     FUN = function(r){do.call(simulation_model, as.list(parameter[r,]))},
                                                     mc.cores = ncore)), nrow = nrow(parameter), byrow = TRUE)
    )

  #### Name columns to be used in other analysis ####
  if(length(colnames(parameter)) < dim(parameter)[2]) {
    warning("Not all parameters have a name. Parameters renamed as P1, P2, ...")
    colnames(parameter) <- sapply(X = 1:dim(parameter)[2], FUN = function(i) {paste0("P", i)})
  }
  colnames(simulation) <- sapply(X = 1:dim(simulation)[2], FUN = function(i) {paste0("S", i)})

  #### Write output either in two files with suffixes or in a single file and return a list####
  if(!is.null(file_output)) {
    if(multi_file) {
      utils::write.table(simulation, file.path(dirname(file_output), paste0(basename(tools::file_path_sans_ext(file_output)), "-simulation.csv")), ...)
      utils::write.table(parameter, file.path(dirname(file_output), paste0(basename(tools::file_path_sans_ext(file_output)), "-parameter.csv")), ...)
    } else {
      utils::write.table(cbind(parameter, simulation), file_output, ...)
    }
  }
  list(parameter = parameter, simulation = simulation, system_time = timing)
}
