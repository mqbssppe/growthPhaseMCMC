\name{mcmcSampler}
\alias{mcmcSampler}
\title{
MCMC sampler
}
\description{
This function implements the Metropolis-Hastings MCMC sampler for individual time-series.
}
\usage{
mcmcSampler(myData, nIter, finalIterationPdf, modelVariance, mhPropRange, 
mhSinglePropRange, movesRange, startPoint, postPar, dName, timeScale, 
burn, iterPerPlotPrefix, priorParameters, L = 3, LRange, tau, 
gammaParameter, saveTheta, Prior = "complexity")
}
\arguments{
  \item{myData}{
observed data.
}
  \item{nIter}{
number of mcmc iterations
}
  \item{finalIterationPdf}{
output folder
}
  \item{modelVariance}{
null
}
  \item{mhPropRange}{
positive integer
}
  \item{mhSinglePropRange}{
positive integer
}
  \item{movesRange}{
null
}
  \item{startPoint}{
positive integer
}
  \item{postPar}{
list of emprirically estimated parameters
}
  \item{dName}{
subject ID
}
  \item{timeScale}{
null
}
  \item{burn}{
burn-in period.
}
  \item{iterPerPlotPrefix}{
null
}
  \item{priorParameters}{
prior parameters.
}
  \item{L}{
null
}
  \item{LRange}{
range of possible values of the number of change-points.
}
  \item{tau}{
real.
}
  \item{gammaParameter}{
real.
}
  \item{saveTheta}{
TRUE.
}
  \item{Prior}{
character.
}
}
\value{
MCMC output.
}
\author{
Panagiotis Papastamoulis
}

