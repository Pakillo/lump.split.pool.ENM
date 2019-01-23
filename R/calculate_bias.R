#' Calculate bias after simulations
#'
#' @param simdf Data frame with simulation results, as produced by \code{\link{run_sims}}.
#'
#' @return A tidy data frame with simulation details and resulting bias for each method tested.
#' @export
#'
calculate_bias <- function(simdf) {

  simdf.redux <- dplyr::select(simdf, -contains(".se"))

  # intercept
  simdf.interc <- dplyr::select(simdf.redux, -contains("slope")) %>%
    mutate(interc.split.bias = interc.split - interc,
           interc.lump.bias = interc.lump - interc,
           interc.mixed.bias = interc.mixed - interc,
           interc.bayesmix.bias = interc.bayesmix - interc,
           interc.pglmm.bias = interc.pglmm - interc) %>%
    tidyr::gather("parameter", "sim.value", interc) %>%
    #tidyr::gather("method", "estimate", interc.split:interc.pglmm) %>%
    tidyr::gather("method", "bias", contains(".bias")) %>%
    mutate(method = gsub("interc.", "", .$method)) %>%
    mutate(method = gsub(".bias", "", .$method)) %>%
    dplyr::select(-contains("interc."))


  # slope
  simdf.slope <- dplyr::select(simdf.redux, -contains("interc")) %>%
    mutate(slope.split.bias = slope.split - slope,
           slope.lump.bias = slope.lump - slope,
           slope.mixed.bias = slope.mixed - slope,
           slope.bayesmix.bias = slope.bayesmix - slope,
           slope.pglmm.bias = slope.pglmm - slope) %>%
    tidyr::gather("parameter", "sim.value", slope) %>%
    #tidyr::gather("method", "estimate", slope.split:slope.pglmm) %>%
    tidyr::gather("method", "bias", contains(".bias")) %>%
    mutate(method = gsub("slope.", "", .$method)) %>%
    mutate(method = gsub(".bias", "", .$method)) %>%
    dplyr::select(-contains("slope."))


  ## bind both
  biasdf <- dplyr::bind_rows(simdf.interc, simdf.slope)

  biasdf

}
