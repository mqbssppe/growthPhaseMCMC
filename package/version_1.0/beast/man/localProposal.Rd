\name{localProposal}
\alias{localProposal}
\title{
Move 3.b
}
\description{
Implements Move 3.b of the Metropolis-Hastings MCMC sampler.
}
\usage{
localProposal(cutPoints, nTime, mhPropRange, startPoint)
}
\arguments{
  \item{cutPoints}{
Current configuration of change-points.
}
  \item{nTime}{
Total number of time-points.
}
  \item{mhPropRange}{
Parameter \eqn{d_2}.
}
  \item{startPoint}{
Integer, at least equal to 2.
}
}
\value{
  \item{newState}{Candidate state of the chain.}
  \item{propRatio}{Proposal ratio.}
}
\author{
Panagiotis Papastamoulis
}

