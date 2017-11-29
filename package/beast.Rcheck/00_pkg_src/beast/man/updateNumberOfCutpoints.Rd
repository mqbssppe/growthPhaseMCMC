\name{updateNumberOfCutpoints}
\alias{updateNumberOfCutpoints}
\title{
Move 1
}
\description{
Update the number of change-points according to Metropolis-Hastings move 1.
}
\usage{
updateNumberOfCutpoints(cutPoints, nTime, startPoint, LRange, birthProbs)
}
\arguments{
  \item{cutPoints}{
Current configuration
}
  \item{nTime}{
Number of time-points
}
  \item{startPoint}{
Optional integer
}
  \item{LRange}{
Range of possible values
}
  \item{birthProbs}{
Birth probabilities
}
}
\value{
  \item{newState }{Candidate state}
  \item{propRatio }{Proposal ratio}
}
\author{
Panagiotis Papastamoulis
}

