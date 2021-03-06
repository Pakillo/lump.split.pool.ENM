% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_models.R
\name{fit_mixed}
\alias{fit_mixed}
\title{Fit mixed model (lme4)}
\usage{
fit_mixed(simdata)
}
\arguments{
\item{simdata}{List generated by \code{\link{simul_data}}.}
}
\value{
A data frame with 4 columns (intercept estimate and standard error, and slope estimate and standard error), and as many rows as taxa.
}
\description{
Fit mixed model (lme4)
}
