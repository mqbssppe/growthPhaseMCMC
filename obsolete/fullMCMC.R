library(mvoutlier)

simMultiIndNormInvGamma <- function(mu, nu, alpha, beta){
	m <- length(mu)
#	sigma2 <- 1/rgamma(n = m, shape = alpha, rate = beta)
#	tha valw na epistrefei mono to meso tis diasporas wste na thewreitai stathero
	sigma2 <- beta/(alpha - 1)
	mySd <- sqrt( sigma2/nu )
	simMeans <- rnorm(n = m, mean = 0, sd = 1)
	simMeans <- mu + simMeans*mySd
	results <- vector("list", length = 2)
	results[[1]] <- simMeans
	results[[2]] <- sigma2
	names(results) <- c("mean", "variance")
	return(results)
}


proposeTheta <- function(thetaOld, tau, alpha, beta){
	# ta alpha beta xrisimeuoun apla gia na epistrefeis kai tin diaspora pou einai fixed
	m <- length(thetaOld)
#	sigma2 <- 1/rgamma(n = m, shape = alpha, rate = beta)
#	tha valw na epistrefei mono to meso tis diasporas wste na thewreitai stathero
	sigma2 <- beta/(alpha - 1)
	results <- vector("list", length = 2)
	results[[1]] <- rnorm(m, mean = thetaOld, sd = tau)
	results[[2]] <- sigma2
	names(results) <- c("mean", "variance")
	return(results)
}



logLikelihoodFullModel <- function(myData, cutPoints, theta, startPoint){
	nTime <- dim(myData)[1]
	augmentedCutPoints <- c(1, cutPoints, nTime )
	sanityCheck <- diff(augmentedCutPoints)
	if(  min(sanityCheck) <= 0 ){
		logL <- log(0)
	}else{
		goldenBoy <- theta$mean[augmentedCutPoints]
		myMeans <- numeric(nTime)
		for (j in 2:length(augmentedCutPoints)){
			subIndex <- augmentedCutPoints[j-1]:augmentedCutPoints[j]
			myMeans[subIndex] <- (goldenBoy[j] - goldenBoy[j-1])/(augmentedCutPoints[j] - augmentedCutPoints[j-1])*(subIndex - augmentedCutPoints[j-1]) + goldenBoy[j-1]
		}
		myMeanMatrix <- matrix(myMeans, nrow = nTime, ncol = dim(myData)[2] )
		sdPerPoint <- sqrt( theta$var )		
		mySDMatrix <- matrix(sdPerPoint, nrow = nTime, ncol = dim(myData)[2] )
		zData <- (myData - myMeanMatrix)/mySDMatrix
		logL <- sum( apply(zData, 1, function(y){ sum(dnorm(y, mean = 0, sd = 1, log = TRUE)) })   )
	}
	return(logL)	
}


logLikelihoodStandardModel <- function(myData, theta){
	nTime <- dim(myData)[1]
	sdPerPoint <- sqrt( theta$var )		
	myMeanMatrix <- matrix(theta$mean, nrow = nTime, ncol = dim(myData)[2] )
	mySDMatrix <- matrix(sdPerPoint, nrow = nTime, ncol = dim(myData)[2] )
	zData <- (myData - myMeanMatrix)/mySDMatrix
	logL <- sum( apply(zData, 1, function(y){ sum(dnorm(y, mean = 0, sd = 1, log = TRUE)) })   )
	return(logL)	
}


partialLogLikelihoodStandardModel <- function(myData, theta, cutPoints){
	nTime <- dim(myData)[1]
	augmentedCutPoints <- c(1, cutPoints, nTime )
	sdPerPoint <- sqrt( theta$var[augmentedCutPoints] )		
	myMeanMatrix <- matrix(theta$mean[augmentedCutPoints], nrow = length(cutPoints) + 2, ncol = dim(myData)[2] )
	mySDMatrix <- matrix(sdPerPoint, nrow = length(cutPoints) + 2, ncol = dim(myData)[2] )
	zData <- (myData[augmentedCutPoints, ] - myMeanMatrix)/mySDMatrix
	logL <- sum( apply(zData, 1, function(y){ sum(dnorm(y, mean = 0, sd = 1, log = TRUE)) })   )
	return(logL)	
}



# startPoint >= 2 !

logPrior <- function(cutPoints, nTime, startPoint){
	augmentedCutPoints <- c(startPoint - 1, cutPoints, nTime )
	sanityCheck <- diff(augmentedCutPoints)
	if(  min(sanityCheck) <= 0 ){logP <- log(0)}else{
	logP <- -log(nTime-2-startPoint) - log(nTime - cutPoints[1] - 2) - log(nTime - cutPoints[2] - 1)}
	return(logP)
}

simulateFromPrior <- function(nTime, startPoint){
	cutPoints <- numeric(3)
	indexRange <- startPoint:(nTime-3)
	if(length(indexRange) == 1){
		cutPoints[1] <- indexRange}else{
		cutPoints[1] <- sample(indexRange, 1)
	}
	indexRange <- (cutPoints[1] + 1):(nTime-2)
	if(length(indexRange) == 1){
			cutPoints[2] = indexRange
		}else{
			cutPoints[2] <- sample(indexRange, 1)
	}
	indexRange <- (cutPoints[2] + 1):(nTime-1)
	if(length(indexRange) == 1){
			cutPoints[3] = indexRange
		}else{
			cutPoints[3] <- sample(indexRange, 1)
	}
	return(cutPoints)
}

localProposal <- function(cutPoints, nTime, mhPropRange, startPoint){
	if(missing(mhPropRange)){mhPropRange = 1}
	#print("local")
	epsilon <- sample(c(-seq(1:mhPropRange), 0, seq(1:mhPropRange)),3, replace = TRUE)
	newState <- cutPoints + epsilon
	augmentedCutPoints <- c(startPoint - 1, newState, nTime )
	sanityCheck <- diff(augmentedCutPoints)
	if(  min(sanityCheck) <= 0 ){	# not valid state!
		propRatio <- log(0)
	}else{
		propRatio <- 0
	}
	
	results <- vector("list", length = 2)
	results[[1]] <- newState
	results[[2]] <- propRatio
	names(results) <- c("newState", "propRatio")
	return(results)
}


singleLocalProposal <- function(cutPoints, nTime, mhSinglePropRange, startPoint){
	if(missing(mhSinglePropRange)){mhSinglePropRange = 1}
	epsilon <- sample(c(-seq(1:mhSinglePropRange), seq(1:mhSinglePropRange)),1)
	justOneIndex <- sample(3,1)
	newState <- cutPoints
	newState[ justOneIndex ] <- cutPoints[ justOneIndex ] + epsilon
	augmentedCutPoints <- c(startPoint - 1, newState, nTime )
	sanityCheck <- diff(augmentedCutPoints)
	if(  min(sanityCheck) <= 0 ){	# not valid state!
		propRatio <- log(0)
	}else{
		propRatio <- 0
	}
	
	results <- vector("list", length = 2)
	results[[1]] <- newState
	results[[2]] <- propRatio
	names(results) <- c("newState", "propRatio")
	return(results)
}


mcmcSampler <- function(myData, nIter, finalIterationPdf, modelVariance, mhPropRange, mhSinglePropRange, movesRange, startPoint, 
		postPar, dName, timeScale, burn, showProgress, iterPerPlotPrefix, priorParameters){
	if(missing(showProgress)){showProgress = FALSE}
	if(missing(timeScale)){timeScale = 1}
	if(missing(mhPropRange)){mhPropRange = 1}
	if (missing(modelVariance)){modelVariance = FALSE}
	if(missing(startPoint)){ startPoint = 2 }
	if(startPoint < 2){ startPoint = 2 }
	if(missing(burn)){burn = 1}
	if(missing(movesRange)){movesRange = as.character(1:3)}
	if (missing(finalIterationPdf)){ finalIterationPdf = NULL }
	nTime <- dim(myData)[1]
	cutPoints <- array(data = NA, dim = c(nIter, 3))
	mu <- array(data = NA, dim = c(nIter, nTime))
	sigma2 <- array(data = NA, dim = c(nIter, nTime))
	iter <- 1
#	cutPoints[iter, ] <- simulateFromPrior(nTime = nTime, startPoint = startPoint)
	cutPoints[iter, ] <- c(startPoint + 1, startPoint + 10, startPoint + 30)
	overkill <- simMultiIndNormInvGamma(mu = postPar[[1]], nu = postPar[[2]], alpha = postPar[[3]], beta = postPar[[4]])
	theta <- overkill
	mu[iter, ] <- overkill$mean
	sigma2[iter, ] <- overkill$var
	lValues <- lPosterior <- numeric(nIter)
	lValues[iter] <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, ], theta = overkill, startPoint = startPoint)
	lPosterior[iter] <- lValues[iter] + logPrior(cutPoints = cutPoints[iter, ], nTime = nTime, startPoint = startPoint) 
	mhRate0a <- 0
	mhRate0 <- 0
	mhRate0c <- 0
	mhRate <- 0
	mhRate2 <- 0
	mhRate3 <- 0
	myPositions <- floor(seq(0, nTime, length = 7))  #c(0,50,100,150,200,250,300)
	myLabels <- round(myPositions*timeScale,1)
	uchar <- myUnicodeCharacters()
	if(showProgress){pdf(file = paste0(finalIterationPdf,"/",dName,".pdf"),width = 9, height = 6)}
	for(iter in 2:nIter){
		lValues[iter] <- lValues[iter - 1]
		cutPoints[iter, ] <- cutPoints[iter - 1, ]
#		Move 0.a: update theta parameters using the standard posterior as proposal		
		overkill <- simMultiIndNormInvGamma(mu = postPar[[1]], nu = postPar[[2]], alpha = postPar[[3]], beta = postPar[[4]])
		propCurrent <- partialLogLikelihoodStandardModel(myData = myData, theta = theta, cutPoints = cutPoints[iter, ])
		propNew <- partialLogLikelihoodStandardModel(myData = myData, theta = overkill, cutPoints = cutPoints[iter, ])
		propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, ], theta = overkill, startPoint = startPoint)
		mhRatio <- propLogL + propCurrent - lValues[iter] - propCurrent
		if( log(runif(1)) < mhRatio ){
			lValues[iter] <- propLogL
			mhRate0a <- mhRate0a + 1
			theta$mean <- overkill$mean
			lValues[iter] <- propLogL
		}




#		Move 0.c: update all theta parameters using a standard symmetric proposal based on the current values	
		overkill <- proposeTheta(thetaOld = theta$mean, tau = 0.01,  alpha = postPar[[3]], beta = postPar[[4]])
		propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, ], theta = overkill, startPoint = startPoint)
		mhRatio <- propLogL - lValues[iter]
		augCutPoints <- c(1, cutPoints[iter, ], nTime )
		zCurrent <- (theta$mean[augCutPoints] - priorParameters$mu0[augCutPoints])/sqrt(theta$variance[augCutPoints]/priorParameters$nu0)
		zProp <- (overkill$mean[augCutPoints] - priorParameters$mu0[augCutPoints])/sqrt(overkill$variance[augCutPoints]/priorParameters$nu0)
		priorRatio <- sum(dnorm( x = zProp, mean = 0, sd = 1, log = TRUE)) - sum(dnorm( x = zCurrent, mean = 0, sd = 1, log = TRUE))
		mhRatio <- mhRatio + priorRatio
		if( log(runif(1)) < mhRatio ){
			lValues[iter] <- propLogL
			mhRate0c <- mhRate0c + 1
			theta <- overkill
		}

		chooseMove <- sample(movesRange, 1)
		if(chooseMove == "1"){
#			Move a: proposal from prior
			newCutPoints <- simulateFromPrior(nTime = nTime, startPoint = startPoint)
			propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
			mhRatio <- propLogL - lValues[iter - 1]
			u <- runif(1)
			if(log(u) < mhRatio){
				cutPoints[iter, ] <- newCutPoints
				lValues[iter] <- propLogL
				mhRate <- mhRate + length(movesRange)
			}
		}
		if(chooseMove == "2"){
#			Move b: local random walk
			lpProp <- localProposal(cutPoints = cutPoints[iter, ] , nTime = nTime, mhPropRange = mhPropRange, startPoint = startPoint)
			newCutPoints <- lpProp$newState
			propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
			mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
				logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
				logPrior(cutPoints = cutPoints[iter, ], nTime = nTime, startPoint = startPoint)

			u <- runif(1)
			if(log(u) < mhRatio){
				cutPoints[iter, ] <- newCutPoints
				lValues[iter] <- propLogL
				mhRate2 <- mhRate2 + length(movesRange)
			}
		}
		if(chooseMove == "3"){
#			Move c: random walk to just one index
			lpProp <- singleLocalProposal(cutPoints = cutPoints[iter, ], nTime = nTime, 
					mhSinglePropRange = mhSinglePropRange, startPoint = startPoint)
			newCutPoints <- lpProp$newState
			propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
			mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
				logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
				logPrior(cutPoints = cutPoints[iter, ], nTime = nTime, startPoint = startPoint)

			u <- runif(1)
			if(log(u) < mhRatio){
				cutPoints[iter, ] <- newCutPoints
				lValues[iter] <- propLogL
				mhRate3 <- mhRate3 + length(movesRange)
			}
		}
		lPosterior[iter] <- lValues[iter] + logPrior(cutPoints = cutPoints[iter, ], nTime = nTime, startPoint = startPoint)
		if(showProgress){
		myStepSize <- 1
		if(iter > 100){myStepSize <- 5}
		if(iter > 200){myStepSize <- 10}
		if(iter > 1000){myStepSize <- 20}
		if(iter > 2000){myStepSize <- 50}
		if(iter > 20000){myStepSize <- 100}
		if(iter > 30000){myStepSize <- 500}
		if(iter %% myStepSize == 0){
#			par(mfrow = c(1,3))
#			plot(lPosterior[1:iter], type= "l", xlab = "mcmc iteration", ylab = "log-Posterior")
#			matplot(cutPoints[1:iter,]*timeScale, type = "l", col = 2:4, xlab = "mcmc iteration", lty = 1, ylab = "hours")
#		pdf(file = paste0(finalIterationPdf,"/",dName,"-iter-",iter,".pdf"),width = 9, height = 6)
			matplot(myData, type = "l", col = 1, xaxt = "n", xlab = "hours", ylab = "growth")
			xPoints <- c(1, cutPoints[iter,], nTime )
			yPoints <- theta$mean[xPoints]
			points(xPoints[1:2], yPoints[1:2], type = "l", col = "red", lwd = 4)
			points(xPoints[2:3], yPoints[2:3], type = "l", col = "green", lwd = 4)
			points(xPoints[3:4], yPoints[3:4], type = "l", col = "blue", lwd = 4)
			points(xPoints[4:5], yPoints[4:5], type = "l", col = "gray", lwd = 4)
			matplot(myData, type = "l", col = 1, add = TRUE)
			axis(1, at = myPositions, labels = myLabels)
			#abline(v = cutPoints[iter, ], col = c(2,3,4))
			#abline(v = startPoint, lty = 4, col = "gray")
			legend("bottomright", paste0("replicate ",1:3), lty = 1:3)
			text(x = cutPoints[iter, ] - 20, y = 0.8, c("phase 1", "phase 2", "phase 3"), col = c(2,3,4))
			if(missing(dName)  == FALSE){
				title( paste0(dName,", iter = ", iter), outer = TRUE , line=-1)
			}
			cat(paste0("  Iteration: ", iter, ". MH acceptance rates: (M1a) ", 
				round(mhRate0a*100/iter,2), "%, (M1b) ",
				round(mhRate0*100/iter,2), "%, (M1c) ",
				round(mhRate0c*100/iter,2), "%, (M2) ", 
				round(mhRate*100/iter,2), "%, (M3) ", 
				round(mhRate2*100/iter,2) ,"%, (M4) ", 
				round(mhRate3*100/iter,2), "%."),"\n")
		}
		}
		else{
			if(iter %% (nIter/2) == 0){cat(paste0(uchar, " "))}
		}
	}
	if(showProgress == FALSE){
		cat(paste0(" Accepted MH moves: [Move 1] ", round(mhRate0*100/iter,2), "%, [Move 2] ", round(mhRate*100/iter,2), "%, [Move 3] ",round(mhRate2*100/iter,2) ,"%, [Move 4] ", round(mhRate3*100/iter,2), "%."),"\n")
	}else{
		dev.off()
	}
	results <- vector("list", length = 3)
	results[[1]] <- cutPoints
	results[[2]] <- lValues
	names(results) <- c("cutPoints", "logLikelihood", "acceptanceRates")

	if(is.null(finalIterationPdf) == FALSE){
		pdf(file = paste0(finalIterationPdf,"/",dName,"_mcmcTrace.pdf"), width = 12, height = 10)
		split.screen( figs = c( 2, 1 ) )
		split.screen( figs = c( 1, 2 ), screen = 1 )
		split.screen( figs = c( 1, 1 ), screen = 2 )
		screen( 3 )
		plot(lPosterior, ylim = range(lPosterior[-(1:burn)]), type= "l", xlab = "mcmc iteration", ylab = "log-posterior pdf")
		points(lPosterior[1:burn], type = "l", lty = 5, col = "white")
		legend("bottomright",lty = 3, col = 1, "burn-in period")

		screen( 4 )
		ylim = range(cutPoints*timeScale)
		ylim[1] = 0
		matplot(cutPoints*timeScale, ylim = ylim, type = "l", col = 2:4, xlab = "mcmc iteration", lty = 1, ylab = "hours")
		matplot(cutPoints[1:burn,]*timeScale, ylim = ylim, type = "l", col = "white", lty = 5, add = TRUE)
		legend("bottomright", c("cut-point of phase 1", "cut-point of phase 2", "cut-point of phase 3"), lty = 1, col = c(2, 3, 4), bty = "n")

		screen( 5 )
		matplot(myData, type = "l", col = 1, xaxt = "n", xlab = "hours", ylab = "growth")
		axis(1, at = myPositions, labels = myLabels)
		abline(v = apply(cutPoints[-(1:burn), ], 2, median), col = c(2,3,4), lwd = 2)
		legend("bottomright", paste0("replicate ",1:3), lty = 1:3, bty = "n")
		legend("topleft", c("posterior median 1", "posterior median 2", "posterior median 3"), lty = c(1,1,1), col = c(2, 3, 4), bty = "n")
		mcmcVar <- apply(cutPoints[-(1:burn),],2,var)
		if(mcmcVar[1] > 0){
			points(density(cutPoints[-(1:burn),1], adjust = 7), type = "l", col = 2, lwd = 4)
			points(density(cutPoints[-(1:burn),1], adjust = 7), type = "l", col = "gray", lwd = 2)
		}
		if(mcmcVar[2] > 0){
			points(density(cutPoints[-(1:burn),2], adjust = 7), type = "l", col = 3, lwd = 4)
			points(density(cutPoints[-(1:burn),2], adjust = 7), type = "l", col = "gray", lwd = 2)
		}
		if(mcmcVar[3] > 0){
			points(density(cutPoints[-(1:burn),3], adjust = 7), type = "l", col = 4, lwd = 4)
			points(density(cutPoints[-(1:burn),3], adjust = 7), type = "l", col = "gray", lwd = 2)
		}
		text(x = cutPoints[iter, ] - 20, y = 0.8, c("phase 1", "phase 2", "phase 3"), col = c(2,3,4))
		if(missing(dName)  == FALSE){
			title( dName, outer = TRUE , line=-1)
		}
		close.screen( all = TRUE )
		dev.off()
	}

	return(results)
}




areaPerPhase <- function(cutPoints, myData, timeScale){
	m <- dim(cutPoints)[1]
	nCutPoints <- dim(cutPoints)[2]
	auc <- array(data = NA, dim = c(m, nCutPoints))
	myMeans <- rowMeans(myData)
	if(missing(timeScale)){timeScale = 1}

	for (i in 1:m){
		xValues <- c(1, cutPoints[i,])
		yValues <- myMeans[xValues]
		auc[i, ] <- diff(yValues)*diff(xValues)*timeScale/2 + c(0, (cutPoints[i,2] - cutPoints[i,1])*yValues[2] , (cutPoints[i,3] - cutPoints[i,3])*yValues[3])*timeScale
	}
	Sresults <- array(data = NA, dim = c(3, 2))
	colnames(Sresults) <- c("mean", "var")
	rownames(Sresults) <- c("area_phase_1", "area_phase_2", "area_phase_3")
	Sresults[,1] <- apply(auc, 2, mean)
	Sresults[,2] <- apply(auc, 2, var)
	results <- vector("list", length=2)
	names(results) <- c("PosteriorSummary", "mcmc_auc")
	results[[1]] <- Sresults
	results[[2]] <- auc
	return(results)
}



ratePerPhase <- function(cutPoints, myData, timeScale){
	m <- dim(cutPoints)[1]
	nCutPoints <- dim(cutPoints)[2]
	myRates <- array(data = NA, dim = c(m, nCutPoints))
	myMeans <- rowMeans(myData)
	if(missing(timeScale)){timeScale = 1}

	for (i in 1:m){
		xValues <- c(1, cutPoints[i,])
		yValues <- myMeans[xValues]
		myRates[i, ] <- diff(yValues)/(diff(xValues)*timeScale)
	}
	Sresults <- array(data = NA, dim = c(3, 2))
	colnames(Sresults) <- c("mean", "var")
	rownames(Sresults) <- c("rate_phase_1", "rate_phase_2", "rate_phase_3")
	Sresults[,1] <- apply(myRates, 2, mean)
	Sresults[,2] <- apply(myRates, 2, var)
	results <- vector("list", length=2)
	names(results) <- c("PosteriorSummary", "mcmc_rates")
	results[[1]] <- Sresults
	results[[2]] <- myRates
	return(results)
}





getVariance <- function(myDataList, blankThreshold){
	if(missing(blankThreshold)){blankThreshold = 0.02}
	n <- dim(myDataList[[1]])[2]
	nTime <- dim(myDataList[[1]])[1]
	nReps <- length(myDataList)
	sdPerPoint <- numeric(nTime)
	j <- 0
	for(i in 1:n){
		myData <- myDataList[[1]][,i]
		for(k in 2:nReps){
			myData <- rbind(myData, myDataList[[k]][,i])
		}
		if(max(myData) > blankThreshold){
			sdPerPoint  <- sdPerPoint  + apply(myData, 2, var)
			j <- j + 1
		}
	}
	sdPerPoint <- sdPerPoint/j
	smallIndex <- which(sdPerPoint < 10^{-5})
	if(length(smallIndex) > 0){
		sdPerPoint[smallIndex] <- rep(10^{-5}, length(smallIndex))
	}
	sdPerPoint <- sqrt(sdPerPoint)
	return(sdPerPoint)
}


growthPhaseFullMCMC <- function(myDataList, burn, nIter, mhPropRange, mhSinglePropRange, movesRange, startPoint, getSDvalues, timeScale, blankThreshold, savePlots, showProgress, zeroNormalization){
#	burn = 2000, nIter = 5000,mhPropRange = 1, mhSinglePropRange = 50, getSDvalues = T, startPoint=54, timeScale = 1/6,
	myColNames <- colnames(myDataList[[1]])
	if(missing(timeScale)){timeScale = 1/6}
	if(missing(showProgress)){showProgress = FALSE}
	if(missing(blankThreshold)){blankThreshold = 0.02}
	if(missing(burn)){burn = 2000}
	if(missing(nIter)){nIter = 5000}
	if(missing(zeroNormalization)){zeroNormalization = TRUE}
	if(burn > nIter - 1){stop("`burn` should be smaller than `nIter`.")}
	if(missing(savePlots)){
		savePlots = NULL; finalIterationPdf = FALSE
	}else{ 
		if(dir.exists(savePlots)){
			myPrompt <- readline(paste0("*  [WARNING] Directory `",getwd(),"/", savePlots, "` exists. Overwrite? 1 = YES, 0 = ABORT "))
			if(myPrompt != 1){stop("Process killed by the user.")}
		}else{
			dir.create(savePlots)
			cat(paste0("*  Plots will be saved to directory: `",getwd(),"/",savePlots,"`"),"\n")
		} 
	}
	if(missing(startPoint)){ startPoint = 48 }
	if( startPoint < 2){ startPoint = 2 }
	if(missing(mhPropRange)){mhPropRange = 1}
	if(missing(mhSinglePropRange)){mhSinglePropRange = 50}
	if(missing(movesRange)){movesRange = as.character(1:3)}
	if(missing(getSDvalues)){getSDvalues = TRUE}
	if(getSDvalues == FALSE){sdValues = NULL}
	if(zeroNormalization){
		cat(paste0("*  Normalizing at time zero... "))	
		myDataList <- normalizeTime0(myDataList = myDataList)
		cat(paste0(" done.","\n"))
	}
	if(getSDvalues == TRUE){
		cat(paste0("*  Computing posterior parameters... "))
		priorParameters <- computeEmpiricalPriorParameters(myDataList = myDataList)
		posteriorParameters <- computePosteriorParameters(myDataList = myDataList, priorParameters = priorParameters)
		cat(paste0(" done.","\n"))
	}
	nReps <- length(myDataList)
	n <- dim(myDataList[[1]])[2]
	cutPoints <- array(data = NA, dim = c(n, 3))
	cutPointsVar <- array(data = NA, dim = c(n, 3))
	areaMeanPerPhase <- array(data=NA, dim = c(n, 3))
	areaVarPerPhase <- array(data=NA, dim = c(n, 3))
	rateMeanPerPhase <- array(data=NA, dim = c(n, 3))
	rateVarPerPhase <- array(data=NA, dim = c(n, 3))
	colnames(areaMeanPerPhase) <- colnames(areaVarPerPhase) <- c("phase_1", "phase_2", "phase_3")
	if(nReps < 2){stop("no replicates")}
	cat(paste0("*  MCMC sampler parameters: nIter = ", nIter, ", burn = ", burn, ", startPoint = ", startPoint ),"\n")
	cat(paste0("*  Running MCMC for ", n, " subjects..."), "\n")	
	NAindex <- c()
	for (i in 1:n){
		cat(paste0("*    i = ", i, ", name: ", myColNames[i], " " ))
		myData <- myDataList[[1]][ , i]
		if( max( myData ) < blankThreshold ){
			cat(paste0("                [SKIPPED] this looks like a blank well.","\n"))
			NAindex <- c(NAindex, i)
		}else{
			for (j in 2:nReps){
				myData <- cbind(myData, myDataList[[j]][ , i])
			}
			posteriorParametersIter <- posteriorParameters
			posteriorParametersIter[[1]] <- posteriorParameters[[1]][i,]
			mhRunForSubject <- mcmcSampler(myData = myData, nIter = nIter, mhPropRange = mhPropRange, dName = myColNames[i], burn = burn,
							mhSinglePropRange = mhSinglePropRange, movesRange = movesRange, finalIterationPdf = savePlots,
							startPoint = startPoint, postPar = posteriorParametersIter, timeScale = timeScale, 
							showProgress = showProgress, iterPerPlotPrefix = paste0("plot-", i), priorParameters = priorParameters)
			cutPoints[i, ] <- apply(mhRunForSubject$cutPoints[-(1:burn), ], 2, median)
			cutPointsVar[i, ] <- apply(mhRunForSubject$cutPoints[-(1:burn), ], 2, var)
			getArea <- areaPerPhase(cutPoints = mhRunForSubject$cutPoints[-(1:burn), ], myData = myData, timeScale = timeScale)$PosteriorSummary
			getRate <- ratePerPhase(cutPoints = mhRunForSubject$cutPoints[-(1:burn), ], myData = myData, timeScale = timeScale)$PosteriorSummary
			rateMeanPerPhase[i, ] <- getRate[,1]
			rateVarPerPhase[i, ] <- getRate[,2]
			areaMeanPerPhase[i, ] <- getArea[,1]
			areaVarPerPhase[i, ] <- getArea[,2]
#			stop("stop")
		}
	}


	results <- vector("list", length = 7)
	results[[1]] <- cutPoints*timeScale
	results[[2]] <- cutPointsVar*(timeScale^2)
	results[[3]] <- rateMeanPerPhase
	results[[4]] <- rateVarPerPhase
	results[[5]] <- areaMeanPerPhase
	results[[6]] <- areaVarPerPhase
	if(length(NAindex) > 0){
		myDF <- results[[1]][-NAindex,]
		rownames(myDF) <- myColNames[-NAindex]
	}else{
		myDF <- results[[1]]
		rownames(myDF) <- myColNames
	}
	cat(paste0("*  Outlier detection at the 0.01 level... "),"\n")
	cat(paste0("*  "),"\n")
	if(missing(savePlots) == FALSE){
		pdf(file = paste0(savePlots,"/outliers_projection.pdf"), width = 18, height = 12)
			mvOut <- aq.plot(myDF, alpha=0.01)
		dev.off()
	}else{
		mvOut <- aq.plot(myDF, alpha=0.01)
	}
	results[[7]] <- mvOut$outliers
	cat(paste0("*  "),"\n")
	cat(paste0("*                                     ... done."),"\n")
	rownames(results[[1]]) <- rownames(results[[2]]) <- rownames(results[[3]]) <- rownames(results[[4]]) <- rownames(results[[5]])<- rownames(results[[6]]) <- myColNames
	names(results) <- c("posteriorMedian", "posteriorVar", "rateMean", "rateVar", "areaMean", "areaVar", "outliers")
	cat(paste0("*  ALL DONE."),"\n")
	if(missing(savePlots) == FALSE){ 
		cat(paste0("*  See produced *.pdf files in: `",getwd(),"/",savePlots,"`"),"\n")
	}



	return(results)

}

myUnicodeCharacters <- function(){
	mySymbols <- c("\U0000A4", "\U0000A3", "\U0000A5", "\U000285","\U0009F8", "\U000E5B","\U001405","\U001518","\U001620","\U0018F2","\U00204D","\U0021AF","\U00220E","\U00261D","\U00262F",
"\U00266C","\U00267B","\U002687","\U002713","\U002730","\U00272A", "\U0027C6","\U002FC2","\U00265E","\U00269D","\U002A12", "\U002605", "\U0029EB", "\U002300", "\U002301", "\U002302", "\U002303", "\U002304", "\U002305", "\U002306", "\U002307", "\U002308", "\U002309", "\U0023F0", "\U0023ED", "\U0023E7", "\U0025F4", "\U0025CD", "\U0025B5", "\U002615", "\U002660", "\U0026C7", "\U002667", "\U002706", "\U00270D", "\U0026F7")
	return(mySymbols[floor(length(mySymbols)*runif(1)+1)])
}


normalizeTime0 <- function(myDataList){

	normalizedData <- myDataList
	nReps <- length(myDataList)
	n <- dim(myDataList[[1]])[2]
	for(i in 1:n){
		for(j in 1:nReps){
			valueAtZero <- myDataList[[j]][1,i]
			normalizedData[[j]][ , i] <- myDataList[[j]][ , i] - valueAtZero
		}
	}
	return(normalizedData)
}


computePosteriorParameters <- function(myDataList, priorParameters, blankThreshold){
	if(missing(blankThreshold)){blankThreshold = 0.02}
	n <- dim(myDataList[[1]])[2]
	nTime <- dim(myDataList[[1]])[1]
	nReps <- length(myDataList)
	mu0 <- priorParameters$mu0
	nu0 <- priorParameters$nu0
	alpha0 <- priorParameters$alpha0
	beta0 <- priorParameters$beta0
	j <- 0
	alive <- c()
	for(i in 1:n){
		myData <- myDataList[[1]][,i]
		for(k in 2:nReps){
			myData <- rbind(myData, myDataList[[k]][,i])
		}
		if(max(myData) > blankThreshold){
			j <- j + 1
			alive <- c(alive, i)
#			xMean <- xMean + colMeans(myData)
		}
	}
	nAlive <- j

	meanParameters1 <- array(data = 0, dim = c(n, nTime))
	meanParameters2 <- nReps + nu0
	varianceParameters1 <- alpha0 + nAlive*nReps/2	
	varianceParameters2 <- numeric(nTime)	
	sumOverN <- numeric(nTime)
	for(i in alive){
		myData <- myDataList[[1]][,i]
		for(k in 2:nReps){
			myData <- rbind(myData, myDataList[[k]][,i])
		}
		cSum <- colSums(myData)
		meanParameters1[i, ] <- (cSum + nu0*mu0)/meanParameters2
		sumOverN <- sumOverN + meanParameters2*colSums(myData^2) - cSum^2 - 2*nu0*mu0*cSum
	}
	varianceParameters2 <- beta0 + 0.5*(nAlive*nReps*nu0*mu0^2 + sumOverN)/meanParameters2
	results <- vector("list", length = 4)
	results[[1]] <- meanParameters1
	results[[2]] <- meanParameters2
	results[[3]] <- varianceParameters1
	results[[4]] <- varianceParameters2
	names(results) <- c("meanPar1", "meanPar2", "varPar1", "varPar2")
	return(results)
	
}

# ara prepei na valw opwsdipote tin empirical prior ston meso.


computeEmpiricalPriorParameters <- function(myDataList, blankThreshold){
	priorParameters <- vector("list",length = 4)
	if(missing(blankThreshold)){blankThreshold = 0.02}
	n <- dim(myDataList[[1]])[2]
	nTime <- dim(myDataList[[1]])[1]
	nReps <- length(myDataList)
	j <- 0
	alive <- c()
	xMean <- numeric(nTime)
	for(i in 1:n){
		myData <- myDataList[[1]][,i]
		for(k in 2:nReps){
			myData <- rbind(myData, myDataList[[k]][,i])
		}
		if(max(myData) > blankThreshold){
			j <- j + 1
			alive <- c(alive, i)
			xMean <- xMean + colMeans(myData)
		}
	}
	nAlive <- j
	xMean <- xMean/nAlive
	
	names(priorParameters) <- c("mu0", "nu0", "alpha0", "beta0")
	priorParameters$mu0 <- xMean
	priorParameters$nu0 <- 0.5
	priorParameters$alpha0 <- 1
	priorParameters$beta0 <- 1
	return(priorParameters)

}


