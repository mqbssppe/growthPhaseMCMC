#-------------------------------------------------------------------------------------------------------------------------------------
# (1) Download source code:
# 	https://github.com/mqbssppe/growthPhaseMCMC/blob/master/fullMCMC_complexityPrior.R
# (2) Download script for generating simulated data:
#	https://github.com/mqbssppe/growthPhaseMCMC/blob/master/simulate.R
# (3) Make sure that 'fullMCMC_complexityPrior.R' and 'simulate.R' are in the working directory
#-------------------------------------------------------------------------------------------------------------------------------------

# simulate data
cat(paste("Simulating... "))
source('simulate.R')
cat(paste("OK."), "\n")
cat("\n")

# Note: the simulated data is now contained to the list 'myDataList'

source('fullMCMC_complexityPrior.R')
# We will run the MCMC sampler for the 7th time series of the simulated dataset
#	which is shown in Figure 2 of the manuscript
subsetIndex <- 7
# Run sampler using the variance estimator in Equation (3.13) of manuscript
sameVariance=TRUE
set.seed(1)
generalSampler3 <- growthPhaseFullMCMC(myDataList = myDataList, burn = 20000, nIter = 70000, 
                        mhSinglePropRange = 50, savePlots = "sim8_same", zeroNormalization = FALSE,
                        showProgress = FALSE, movesRange = c(2, 3), L = 3, LRange = 0:30, tau = 0.05, 
                        gammaParameter = 2, nu0 = 1/200, alpha0 = 0.1, beta0 = 0.1, saveTheta = TRUE, subsetIndex = subsetIndex, sameVariance=sameVariance)

# Run sampler using the variance estimator in Equation (3.12) of manuscript
sameVariance=FALSE
set.seed(2)
generalSampler4 <- growthPhaseFullMCMC(myDataList = myDataList, burn = 20000, nIter = 70000, 
                        mhSinglePropRange = 50, savePlots = "sim8_diff", zeroNormalization = FALSE,
                        showProgress = FALSE, movesRange = c(2, 3), L = 3, LRange = 0:30, tau = 0.05, 
                        gammaParameter = 2, nu0 = 1/200, alpha0 = 0.1, beta0 = 0.1, saveTheta = TRUE, subsetIndex = subsetIndex, sameVariance=sameVariance)
# clean tmp directories
system("rm -r sim8_same sim8_diff")


# function returning a sequence of points in the line between two points in R^2 

findLine <- function(x, y){
	a <- (y[2] - y[1])/(x[2] - x[1])
	b <- y[1] - a*x[1]
	return(
		a*seq(x[1], x[2], by = 1) + b
	)
}


# Return the variance point estimates per time-point according to both estimators

priorPar <- computeEmpiricalPriorParameters(myDataList = myDataList, nu0 = 1/200, alpha0 = 0.1, beta0 = 0.1)
postPar  <- computePosteriorParametersFree(myDataList = myDataList, priorParameters = priorPar)
postParSame  <- computePosteriorParameters(myDataList = myDataList, priorParameters = priorPar)


# Figure 2 of the manuscript
cat(paste("Producing Figure 2 of the manuscript... "))
pdf(file = "fig2_up.pdf", width = 12, height = 5)
	layout(matrix(c(1,2), ncol = 2), width = c(1.2,3));
	i <- subsetIndex
	plot(generalSampler3$nCutPointsTrace[[i]], type = "l", xlab = "MCMC iteration", ylab = "# change-points")
	m <- dim(generalSampler3$Cutpoint_mcmc_trace_map[[i]])[1]
	Ncuts <- dim(generalSampler3$Cutpoint_mcmc_trace_map[[i]])[2]

	matplot(cbind(myDataList[[1]][,i], myDataList[[2]][,i], myDataList[[3]][,i]), type = "n", col = 2:4, lty = 1, xlab = "t", ylab = "x",ylim = c(-15, 40))
	yMean <- rep(0,10)
	for(iter in floor(seq(1,m,length = 500))){
		x <- c(1, generalSampler3$Cutpoint_mcmc_trace_map[[i]][iter,], nTime)
		y <- generalSampler3$theta[[i]][iter,x]
		yMean <- yMean + y
	}
	yMean <- yMean/500
	xMedian <- c(1, generalSampler3$Cutpoint_posterior_median[[i]], nTime)
	xMean <- floor(c(1, colMeans(generalSampler3$Cutpoint_mcmc_trace_map[[i]]), nTime) + 0.5)
	j <- 2
	myLine <- findLine(x = xMean[c(j-1, j)], y = yMean[c(j-1, j)])
	for (j in 3:10){
		v <- findLine(x = xMean[c(j-1, j)], y = yMean[c(j-1, j)])
		myLine <- c(myLine, v[-1])
	}
	upPoints <- myLine + 2*sqrt( postParSame[[4]]/(postParSame[[3]] - 1) )
	downPoints <- myLine - 2*sqrt( postParSame[[4]]/(postParSame[[3]] - 1) )
	upTPoints <- myLine + sqrt( sigma2[i,] )
	downTPoints <- myLine - sqrt( sigma2[i,] )
	matplot(cbind(upPoints, downPoints), type = "l", col = "gray", lty = 1, add = TRUE)
	polygon(x = c(1:nTime, nTime:1), y = c(upPoints, rev(downPoints)), col = "khaki", border = NA)
	matplot(cbind(myDataList[[1]][,i], myDataList[[2]][,i], myDataList[[3]][,i]), type = "l", col = 2:4, lty = 1, add = TRUE)
	points(xMean, yMean, type = "l", lty = 1, col = "gray", lwd = 2)
	abline(v = generalSampler3$Cutpoint_posterior_median[[i]], lty = 3, col = "gray")
	boxplot(generalSampler3$Cutpoint_mcmc_trace_map[[i]], add = TRUE, horizontal=TRUE, boxwex = 4, col = "gray", xaxt = "n", yaxt = "n", at = rep(-13,8))
dev.off()
	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#
	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#	#
pdf(file = "fig2_down.pdf", width = 12, height = 5)
	layout(matrix(c(1,2), ncol = 2), width = c(1.2,3));
	i <- subsetIndex
	plot(generalSampler4$nCutPointsTrace[[i]], type = "l", xlab = "MCMC iteration", ylab = "# change-points")
	m <- dim(generalSampler4$Cutpoint_mcmc_trace_map[[i]])[1]
	Ncuts <- dim(generalSampler4$Cutpoint_mcmc_trace_map[[i]])[2]

	matplot(cbind(myDataList[[1]][,i], myDataList[[2]][,i], myDataList[[3]][,i]), type = "n", col = 2:4, lty = 1, xlab = "t", ylab = "x",ylim = c(-15, 40))
	yMean <- rep(0,10)
	for(iter in floor(seq(1,m,length = 500))){
		x <- c(1, generalSampler4$Cutpoint_mcmc_trace_map[[i]][iter,], nTime)
		y <- generalSampler4$theta[[i]][iter,x]
		yMean <- yMean + y
	}
	yMean <- yMean/500
	xMedian <- c(1, generalSampler4$Cutpoint_posterior_median[[i]], nTime)
	xMean <- floor(c(1, colMeans(generalSampler4$Cutpoint_mcmc_trace_map[[i]]), nTime) + 0.5)
	j <- 2
	myLine <- findLine(x = xMean[c(j-1, j)], y = yMean[c(j-1, j)])
	for (j in 3:10){
		v <- findLine(x = xMean[c(j-1, j)], y = yMean[c(j-1, j)])
		myLine <- c(myLine, v[-1])
	}
	upPoints <- myLine + 2*sqrt( postPar[[4]][i,]/(postPar[[3]] - 1) )
	downPoints <- myLine - 2*sqrt( postPar[[4]][i,]/(postPar[[3]] - 1) )
	upTPoints <- myLine + sqrt( sigma2[i,] )
	downTPoints <- myLine - sqrt( sigma2[i,] )
	matplot(cbind(upPoints, downPoints), type = "l", col = "gray", lty = 1, add = TRUE)
	polygon(x = c(1:nTime, nTime:1), y = c(upPoints, rev(downPoints)), col = "khaki", border = NA)
	matplot(cbind(myDataList[[1]][,i], myDataList[[2]][,i], myDataList[[3]][,i]), type = "l", col = 2:4, lty = 1, add = TRUE)
	points(xMean, yMean, type = "l", lty = 1, col = "gray", lwd = 2)
	abline(v = generalSampler4$Cutpoint_posterior_median[[i]], lty = 3, col = "gray")
	boxplot(generalSampler4$Cutpoint_mcmc_trace_map[[i]], add = TRUE, horizontal=TRUE, boxwex = 4, col = "gray", xaxt = "n", yaxt = "n", at = rep(-13,8))
dev.off()
cat(paste("done. See `fig2_up.pdf` and `fig2_down.pdf`."), "\n")

