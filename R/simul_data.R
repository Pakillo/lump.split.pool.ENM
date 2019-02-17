

#' Simulate dataset
#'
#' @param nspp Number of taxa
#' @param nsite Number of sites
#' @param seed Random seed
#'
#' @return A dataframe
#' @export
#'


simul_data <- function(nspp = NULL, nsite = NULL, seed = NULL) {


  run.id <- gsub(":", "-", as.character(Sys.time()))

  ## random seed
  if (!is.null(seed)) {
    .GlobalEnv$.Random.seed <- seed
  } else {
    seed <- .GlobalEnv$.Random.seed
  }

  dump("seed", file = paste0("analysis/seeds_v2/", run.id, ".R"))

  # see http://www.cookbook-r.com/Numbers/Saving_the_state_of_the_random_number_generator/



  ## Generate environmental site variable
  env <- seq(-10, 10, length.out = nsite)  # mean(env) = 0



  ## Following 'pez' vignette:
  #
  # simulate a phylogenetic tree
  # phy <- ape::rtree(n = nspp)
  # phy <- ape::compute.brlen(phy, method = "Grafen", power = 1)
  #
  # # standardize the phylogenetic covariance matrix to have determinant 1
  # Vphy <- ape::vcv(phy)
  # Vphy.std <- Vphy/(det(Vphy)^(1/nspp))
  #
  # # Perform a Cholesky decomposition of Vphy
  # iD <- t(chol(Vphy.std))
  #
  # ## Set up species-specific regression coefficients as random effects
  # intercept <- iD %*% runif(nspp, -2, 2)
  # slope <- iD %*% runif(nspp, -0.4, 0.1)
  #
  # intercept <- intercept[gtools::mixedorder(rownames(intercept)), ]
  # slope <- slope[gtools::mixedorder(rownames(slope)), ]



  ## Using different approach to simulate continuous trait along phylogeny:

  phy <- ape::rcoal(n = nspp)  # simulate tree
  #plot(phy)

  # simulate intercept values following Brownian Motion
  sim.interc <- ape::rTraitCont(phy, model = "BM", sigma = 1, root.value = runif(1, -1, 1))
  K.interc <- phytools::phylosig(phy, sim.interc)  # phylogenetic signal (Blomberg's K)
  #sim.interc

  # simulate slope values following Brownian Motion
  # (note here intercept and slopes are uncorrelated:
  # check out ape:rTraitMult or mvMORPH::mvSIM for simulating traits correlated evolution)
  sim.slope <- ape::rTraitCont(phy, model = "BM", sigma = 0.2, root.value = runif(1, -0.20, -0.10))
  K.slope <- phytools::phylosig(phy, sim.slope)
  #sim.slope


  intercept <- sim.interc[gtools::mixedorder(names(sim.interc))]
  slope <- sim.slope[gtools::mixedorder(names(sim.slope))]



  ## Calculate suitabilities for each taxa and site
  suitab <- rep(intercept, each = nsite)
  suitab <- suitab + rep(slope, each = nsite) * rep(env, nspp)
  suitab.error <- suitab + rnorm(nspp * nsite, mean = 0, sd = 1) #add some random 'error'
  suitab.invlogit <- arm::invlogit(suitab)
  suitab.error.invlogit <- arm::invlogit(suitab.error)


  # Generate presence-absence data from suitabilities
  pres <- rbinom(length(suitab.error), size = 1, prob = suitab.error.invlogit)


  data2model <- data.frame(taxon = paste("t", sort(rep(1:nspp, nsite)), sep = ""),
                           site = rep(1:nsite, nspp),
                           env = rep(env, nspp),
                           suitab.invlogit = suitab.invlogit,
                           presabs = pres,
                           stringsAsFactors = FALSE)

  simdf <- data.frame(run.id = run.id,
                      nspp = nspp,
                      nsite = nsite,
                      taxon = paste("t", 1:nspp, sep = ""),
                      interc = intercept,
                      K.interc = K.interc,
                      slope = slope,
                      K.slope = K.slope,
                      stringsAsFactors = FALSE)

  simdata <- list(data2model = data2model,
                  simdf = simdf,
                  phylog = phy,
                  seed = seed)

}
