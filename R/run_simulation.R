#' Run simulation
#'
#' @param delta adapt_delta parameter for brm
#' @inheritParams fit_split
#'
#' @return A data.frame
#' @export
#'
run_simulation <- function(nspp, nsite, delta = 0.95) {

  simdata <- simul_data(nspp, nsite)

  split <- fit_split(simdata)

  lump <- fit_lump(simdata)

  mixed <- fit_mixed(simdata)

  bayesmix <- fit_bayesmix(simdata, delta)

  pglmm <- fit_pglmm(simdata, delta)

  results <- data.frame(simdata$simdf, split, lump, mixed, bayesmix, pglmm)

  results

}
