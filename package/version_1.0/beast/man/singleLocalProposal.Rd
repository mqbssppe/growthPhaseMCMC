\name{singleLocalProposal}
\alias{singleLocalProposal}
\title{
Move 3.b
}
\description{
Implement Metropolis-Hastings Move 3.b.
}
\usage{
singleLocalProposal(cutPoints, nTime, mhSinglePropRange, startPoint)
}
\arguments{
  \item{cutPoints}{
Current state
}
  \item{nTime}{
Number of time-points
}
  \item{mhSinglePropRange}{
Prior parameter.
}
  \item{startPoint}{
Optional.
}
}
\value{
  \item{newState }{Candidate state}
  \item{propRatio }{Proposal ratio}
}
\author{
Panagiotis Papastamoulis
}

