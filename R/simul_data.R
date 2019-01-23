

#' Simulate dataset
#'
#' @param nspp Number of taxa
#' @param nsite Number of sites
#' @param seed Random seed
#'
#' @return A dataframe
#' @export
#'


simul_data <- function(nspp, nsite) {


  run.id <- gsub(":", "-", as.character(Sys.time()))

  ## random seed
  if (is.null(.GlobalEnv$.Random.seed)) {
    stop("Please set seed before.\n")
  } else {
    seed <- .GlobalEnv$.Random.seed
    dump("seed", file = paste0("seeds/", run.id, ".R"))
  }
  # see http://www.cookbook-r.com/Numbers/Saving_the_state_of_the_random_number_generator/




  # simulate a phylogenetic tree
  phy <- ape::rtree(n = nspp)
  phy <- ape::compute.brlen(phy, method = "Grafen", power = 0.5)

  # standardize the phylogenetic covariance matrix to have determinant 1
  Vphy <- ape::vcv(phy)
  Vphy.std <- Vphy/(det(Vphy)^(1/nspp))

  # Perform a Cholesky decomposition of Vphy
  iD <- t(chol(Vphy.std))

  # Generate environmental site variable
  env <- seq(-10, 10, length.out = nsite)  # mean(env) = 0



  ## Set up species-specific regression coefficients as random effects
  intercept <- iD %*% runif(nspp, -1, 2)
  slope <- iD %*% runif(nspp, -0.4, 0.1)

  intercept <- intercept[gtools::mixedorder(rownames(intercept)), ]
  slope <- slope[gtools::mixedorder(rownames(slope)), ]



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
                      slope = slope,
                      stringsAsFactors = FALSE)

  simdata <- list(data2model = data2model,
                  simdf = simdf,
                  phylog = phy,
                  seed = seed)

}
