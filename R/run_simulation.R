#' Run simulation
#'
#' @param nspp Number of taxa to simulate.
#' @param nsite Number of sites (equal for all taxa)
#' @param delta 'adapt_delta' parameter for brm
#' @param run.pglmm Logical. Run PGLMM? (default is TRUE).
#'
#'
#' @return A data.frame
#' @export
#'
run_simulation <- function(nspp = NULL, nsite = NULL, delta = NULL, run.pglmm = TRUE) {

  simdata <- simul_data(nspp, nsite)

  split <- fit_split(simdata)

  lump <- fit_lump(simdata)

  mixed <- fit_mixed(simdata)

  bayesmix <- fit_bayesmix(simdata, delta)

  results <- data.frame(simdata$simdf, split, lump, mixed, bayesmix)

  if (isTRUE(run.pglmm)) {

    pglmm <- fit_pglmm(simdata, delta)

    results <- data.frame(simdata$simdf, split, lump, mixed, bayesmix, pglmm)
  }


  results

}



#' Run a number of simulations with fixed parameters
#'
#' @param nsim Number of replicate simulations to run
#' @param force.run Logical. If FALSE (the default) simulations will NOT be run if there is a simulation output file available with same parameters (nspp, nsite). If TRUE, simulations will run and the file will be overwritten.
#' @inheritParams run_simulation
#'
#' @return A data.frame with nrow = nspp*nsim
#' @export
#'

run_sims <- function(nsim = 10, nspp = NULL, nsite = NULL, delta = NULL,
                     run.pglmm = TRUE, force.run = FALSE) {

  if (!file.exists(paste0("simulations_v2/", "nsp", nspp, "_nsite", nsite,".rds")) | isTRUE(force.run)) {

    if (file.exists("pglmm.rds")) file.remove("pglmm.rds")  # delete file to force model compiling for first run

    reps.list <- replicate(nsim, run_simulation(nspp, nsite, delta, run.pglmm), simplify = FALSE)

    reps.df <- do.call("rbind", reps.list)

    saveRDS(reps.df, paste0("simulations_v2/", "nsp", nspp, "_nsite", nsite,".rds"))

    reps.df

  }

}


# run_sims_params <- function(nsim = 10, nspp = NULL, nsite = NULL, delta = NULL,
#                             run.pglmm = TRUE, force.run = FALSE) {
#
#   sims <- purrr::map2_dfr(as.list(nspp), as.list(nsite), run_sims, nsim, delta,
#                           run.pglmm, force.run)
#
# }
