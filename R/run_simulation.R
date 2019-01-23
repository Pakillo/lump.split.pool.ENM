#' Run simulation
#'
#' @param delta 'adapt_delta' parameter for brm
#' @inheritParams fit_split
#'
#' @return A data.frame
#' @export
#'
run_simulation <- function(nspp, nsite, delta) {

  simdata <- simul_data(nspp, nsite)

  split <- fit_split(simdata)

  lump <- fit_lump(simdata)

  mixed <- fit_mixed(simdata)

  bayesmix <- fit_bayesmix(simdata, delta)

  pglmm <- fit_pglmm(simdata, delta)

  results <- data.frame(simdata$simdf, split, lump, mixed, bayesmix, pglmm)

  results

}



#' Run a number of simulations with fixed parameters
#'
#' @param nsim Number of replicate simulations to run
#' @inheritParams run_simulation
#'
#' @return A data.frame with nrow = nspp*nsim
#' @export
#'

run_sims <- function(nsim, nspp, nsite, delta) {

  file.remove("pglmm.rds")  # delete file to force model compiling for first run

  reps.list <- replicate(nsim, run_simulation(nspp, nsite, delta), simplify = FALSE)

  reps.df <- do.call("rbind", reps.list)

}


# gives error with brms:
# run_sims_params <- function(nsim, nspp, nsite, delta) {
#
#   sims <- purrr::map2_dfr(as.list(nspp), as.list(nsite), run_sims_fixed, nsim, delta)
#
# }
