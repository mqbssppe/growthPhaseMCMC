\name{complexityPrior}
\alias{complexityPrior}
\title{
Complexity prior distribution
}
\description{
This function computes the complexity prior distribution on the number of change-points, defined as \eqn{f(\ell) = P(\ell_n = \ell)\propto e^{-\alpha\ell\log(bT/\ell)}, a, b > 0; \ell = 0,1,2,\ldots}. Note that this distribution has exponential decrease (Castillo and van der Vaart, 2012) when \eqn{b>1+e}, so we set \eqn{b=3.72}.
}
\usage{
complexityPrior(Lmax = 20, gammaParameter, nTime)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lmax}{
	maximum number of change-points (default = 20). 
}
  \item{gammaParameter}{
     positive real number, corresponding to \eqn{\alpha}.
}
  \item{nTime}{
	positive integer denoting the total number of time-points.
}
}
\value{
 \item{logPrior}{Prior distribution values in the log-scale.}
}
\references{
Castillo I. and van der Vaart A (2012). Needles and Straw in a Haystack: {P}osterior concentration for possibly sparse sequences. The Annals of Statistics, 40(4), 2069--2101.
}
\author{
Panagiotis Papastamoulis
}

