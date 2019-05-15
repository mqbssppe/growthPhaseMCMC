\name{beast-package}
\alias{beast-package}
\alias{beast-package}
\docType{package}
\title{
\packageTitle{beast}
}
\description{
\packageDescription{beast}
}
\details{
The \code{beast} package deals with \bold{B}ayesian \bold{e}stimation of ch\bold{a}nge-points in the \bold{s}lope of multivariate \bold{t}ime-series, introduced by Papastamoulis et al (2017). For a given period \eqn{t = 1,\ldots,T} we observe multiple time series which are assumed independent, each one consisting of multiple measurements (replicates). Each time series is assumed to have its own segmentation, which is common among its replicates. Thus, different time series have distinct mean parameters in the underlying normal distribution. The variance, which is assumed known, can be either shared between different time series or not and in practice it is estimated at a pre-processing stage. 

Our model infers the joint posterior distribution of the number and location of change-points by MCMC sampling. For this purpose a Metropolis-Hastings MCMC sampler is used. The main function of the package is \code{\link{beast}}. 

Assume that the observed data consists of \eqn{N} time-series, each one consisting of \eqn{R} variables (or replicates) measured at \eqn{T} consecutive time-points. The input of the main function \code{\link{beast}} should be given in the form of a \code{list} \code{myDataList} with the following attributes: \itemize{
\item \code{length(myDataList)} should be equal to \eqn{R}, that is, the number of variables (or replicates)
\item \code{dim(myDataList)[[1]]}, \eqn{\ldots}, \code{dim(myDataList)[[R]]} should be all \eqn{T\times N}, that is, rows and columns should correspond to time-points and different series, respectively.
}
Then, a basic usage of the package consists of the following commands: \itemize{
\item \code{beastRun <- beast( myDataList = myDataList )}
\item \code{print(beastRun)}
\item \code{plot(beastRun)}
}
which correspond to running the MCMC sampler, printing and plotting output summaries, respectively.
}
\author{
\packageAuthor{beast}

Maintainer: \packageMaintainer{beast}
}
\references{
Papastamoulis P., Furukawa T., van Rhijn N., Bromley M., Bignell E. and Rattray M. (2017). Bayesian detection of piecewise linear trends in replicated time-series with application to growth data modelling. arXiv:\bold{1709.06111} [stat.AP]
}
\keyword{ package }
\seealso{
\code{\link{beast}, \link{print.beast.object}, \link{plot.beast.object}}
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

