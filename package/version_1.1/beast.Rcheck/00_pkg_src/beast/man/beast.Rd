\name{beast}
\alias{beast}
\title{
Main function
}
\description{
This is the main function of the package, implementing the MCMC sampler described in Papastamoulis et al (2017). 
}
\usage{
beast(myDataList, burn, nIter, mhPropRange, mhSinglePropRange, startPoint,
    timeScale, savePlots, zeroNormalization, LRange, tau,
    gammaParameter, nu0, alpha0, beta0, subsetIndex, saveTheta, sameVariance,
    Prior
)
}
\arguments{
  \item{myDataList}{
Observed data in the form of a \code{list} with length \eqn{R}, denoting the dimensionality of the multivariate time-series data. For each \eqn{r = 1,\ldots,R}, \code{myDataList[[r]]} should correspond to \eqn{T\times N} array, with \code{myDataList[[r]][t, n]} corresponding to the observed data for time-series \eqn{n = 1,\ldots,N} and time-point \eqn{t = 1,\ldots,T}.
}
  \item{burn}{
Number of iterations that will be discarder as burn-in period. Default value: \code{burn = 20000}.
}
  \item{nIter}{
Number of MCMC iterations. Default value: \code{nIter = 70000}.
}
  \item{mhPropRange}{
Positive integer corresponding to the parameter \eqn{d_1} of MCMC Move 3.a of Papastamoulis et al (2017). Default value: \code{mhPropRange = 1}.
}
  \item{mhSinglePropRange}{
Positive integer denoting the parameter \eqn{d_2} of Papastamoulis et al (2017). Default value: \code{mhPropRange = 40}.
}
  \item{startPoint}{
An (optional) positive integer pointing at the minimum time-point where changes are allowed to occur. Default value: \code{startPoint = 2} (all possible values are taken into account).
}
  \item{timeScale}{
Null.
}
  \item{savePlots}{
Character denoting the name of the folder where various plots will be saved to.
}
  \item{zeroNormalization}{
Logical value denoting whether to normalize to zero all time time-series for \eqn{t=1}. Default: \code{zeroNormalization = FALSE}.
}
  \item{LRange}{
Range of possible values for the number of change-points. Default value: \code{LRange = 0:30}.
}
  \item{tau}{
Positive real number corresponding to parameter \eqn{c} in Move 2 of Papastamoulis et al (2017). Default value: \code{tau = 0.05}.
}
  \item{gammaParameter}{
Positive real number corresponding to parameter \eqn{\alpha} of the exponential prior distribution. Default value: \code{gammaParameter = 2}. 
}
  \item{nu0}{
Positive real number corresponding to prior parameter \eqn{\nu_0} in Papastamoulis et al (2017). Default value: \code{nu0 = 0.1}.
}
  \item{alpha0}{
Positive real number corresponding to prior parameter \eqn{\alpha_0} in Papastamoulis et al (2017). Default value: \code{alpha0 = 1}.
}
  \item{beta0}{
Positive real number corresponding to prior parameter \eqn{\beta_0} in Papastamoulis et al (2017). Default value: \code{beta0 = 1}.
}
  \item{subsetIndex}{
Optional subset of integers corresponding to time-series indexes. If not null, the sampler will be applied only to the specified subset.
}
  \item{saveTheta}{
Logical value indicating whether to save the generated values of the mean per time-point across the MCMC trace. Default: FALSE.
}
  \item{sameVariance}{
Logical value indicating whether to assume the same variance per time-point across time-series. Default value: \code{sameVariance = TRUE}.
}
  \item{Prior}{
Character string specifying the prior distribution of the number of change-points. Allowed values: \code{Prior = "complexity"} (default) or \code{Prior = "Poisson"} (not suggested). 
}
}
\value{
The output of the sampler is returned as a list, with the following features:
\item{Cutpoint_posterior_median}{The estimated medians per change-point, conditionally on the most probable number of change-points per time-series.}
\item{Cutpoint_posterior_variance}{The estimated variances per change-points, conditionally on the most probable number of change-points per time-series.}
\item{NumberOfCutPoints_posterior_distribution}{Posterior distributions of number of change-points per time-series.}
\item{NumberOfCutPoints_MAP}{The most probable number of change-points per time-series.}
\item{Metropolis-Hastings_acceptance_rate}{Acceptance of the MCMC move-types.}
\item{subject_ID}{the identifier of individual time-series.}
\item{Cutpoint_mcmc_trace_map}{The sampled values of each change-point per time series, conditionally on the MAP values.}
\item{theta}{The sampled values of the means per time-series, conditionally on the MAP values.}
\item{nCutPointsTrace}{The sampled values of the number of change-points, per time-series.}
}
\note{
The complexity prior distribution with parameter \code{gammaParameter = 2} is the default prior assumption imposed on the number of change-points. Smaller (larger) values of \code{gammaParameter} will a-priori support larger (respectively: smaller) number of change-points. 

For completeness purposes, the Poisson distribution is also allowed in the \code{Prior}. In this latter case, the \code{gammaParameter} denotes the rate parameter of the Poisson distribution. Note that in this case the interpretation of \code{gammaParameter} is reversed: Smaller (larger) values of \code{gammaParameter} will a-priori support smaller (respectively: larger) number of change-points.
}
\references{
Papastamoulis P., Furukawa T., van Rhijn N., Bromley M., Bignell E. and Rattray M. (2017). Bayesian detection of piecewise linear trends in replicated time-series with application to growth data modelling. arXiv:1709.06111 [stat.AP]
}
\author{
Panagiotis Papastamoulis
}
\examples{
# toy-example (MCMC iterations not enough)
library('beast')	# load package
data("FungalGrowthDataset")	# load dataset
myIndex <- c(392, 62, 3, 117)	# run the sampler only for the 
#                                 specific subset of time-series
set.seed(1)	
# Run MCMC sampler with very small number of iterations (nIter):
run_mcmc <- beast(myDataList = FungalGrowthDataset, subsetIndex = myIndex, 
			zeroNormalization = TRUE, nIter = 40, burn = 20) 
# Print output:
print(run_mcmc)
# Plot output to file: "beast_plot.pdf"
plot(run_mcmc, fileName = "beast_plot_toy.pdf", timeScale=1/6, xlab = "hours", ylab = "growth")
# Run the following commands to obtain convergence:

\dontrun{
# This example illustrates the package using a subset of four 
#      time-series of the fungal dataset. 
library('beast')	# load package
data("FungalGrowthDataset")	# load dataset
myIndex <- c(392, 62, 3, 117)	# run the sampler only for the 
#                                 specific subset of time-series
set.seed(1)		# optional
# Run MCMC sampler with the default number of iterations (nIter =70000):
run_mcmc <- beast(myDataList = FungalGrowthDataset, subsetIndex = myIndex, 
			zeroNormalization = TRUE) 
# Print output:
print(run_mcmc)
# Plot output to file: "beast_plot.pdf"
plot(run_mcmc, fileName = "beast_plot.pdf", timeScale=1/6, xlab = "hours", ylab = "growth")
# NOTE 1: for a complete analysis remove the `subsetIndex = myIndex` argument.
# NOTE 2: `zeroNormalization = TRUE` is an optional argument that forces all 
#	   time-series to start from zero. It is not supposed to be used 
#	   for other applications.
}
}

