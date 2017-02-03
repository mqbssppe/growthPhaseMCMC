library(mvoutlier)

logLikelihood <- function(myData, cutPoints, sdPerPoint, modelVariance, startPoint){
	if(missing(modelVariance)){modelVariance = FALSE}
	if(missing(sdPerPoint)){modelVariance = TRUE}
	nTime <- dim(myData)[1]
	augmentedCutPoints <- c(startPoint - 1, cutPoints, nTime )
	sanityCheck <- diff(augmentedCutPoints)
	if(  min(sanityCheck) <= 0 ){
		logL <- log(0)
	}else{
		goldenBoy <- numeric(length(cutPoints) + 2)
		goldenBoy[1] <- mean(myData[1,])
		goldenBoy[length(augmentedCutPoints)] <- mean(myData[nTime,])
		j <- 1
		for (t in cutPoints){
			j <- j + 1
			goldenBoy[j] <- mean(myData[t,])
		}
		myMeans <- numeric(nTime)
		for (j in 2:length(augmentedCutPoints)){
			subIndex <- augmentedCutPoints[j-1]:augmentedCutPoints[j]
			myMeans[subIndex] <- (goldenBoy[j] - goldenBoy[j-1])/(augmentedCutPoints[j] - augmentedCutPoints[j-1])*(subIndex - augmentedCutPoints[j-1]) + goldenBoy[j-1]
		}
		myMeanMatrix <- matrix(myMeans, nrow = nTime, ncol = dim(myData)[2] )
		if(modelVariance == TRUE){
			sdPerPoint <- sqrt(rowMeans( (myData - myMeanMatrix)^2 ))
			sdPerPoint <- (sdPerPoint < 10^{-7})*10^{-7} + sdPerPoint
		}
		mySDMatrix <- matrix(sdPerPoint, nrow = nTime, ncol = dim(myData)[2] )
		zData <- (myData - myMeanMatrix)/mySDMatrix
		logL <- sum( apply(zData, 1, function(y){ sum(dnorm(y, mean = 0, sd = 1, log = TRUE)) })   )
	}
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


mhSampler <- function(myData, nIter, finalIterationPdf, modelVariance, mhPropRange, mhSinglePropRange, movesRange, startPoint, sdValues, dName, timeScale, burn, showProgress){
	if(missing(showProgress)){showProgress = FALSE}
	if(missing(timeScale)){timeScale = 1}
	if(missing(mhPropRange)){mhPropRange = 1}
	if (missing(modelVariance)){modelVariance = FALSE}
	if(missing(startPoint)){ startPoint = 2 }
	if(startPoint < 2){ startPoint = 2 }
	if(missing(burn)){burn = 1}
	if(missing(movesRange)){movesRange = as.character(1:3)}
	if(is.null(sdValues)){
		cat(paste0("WARNING: sdValues are not provided."),"\n")
		sdPerPoint <- sqrt(apply(myData, 1, function(y){max(var(y),10^{-7})}))
	}else{
		sdPerPoint <- sdValues		
	}
	if (missing(finalIterationPdf)){ finalIterationPdf = NULL }
	nTime <- dim(myData)[1]
	cutPoints <- array(data = NA, dim = c(nIter, 3))
	iter <- 1
	cutPoints[iter, ] <- simulateFromPrior(nTime = nTime, startPoint = startPoint)
	lValues <- lPosterior <- numeric(nIter)
	lValues[iter] <- logLikelihood(myData = myData, cutPoints = cutPoints[iter, ], sdPerPoint = sdPerPoint, modelVariance = modelVariance, startPoint = startPoint)
	lPosterior[iter] <- lValues[iter] + logPrior(cutPoints = cutPoints[iter, ], nTime = nTime, startPoint = startPoint) 
	mhRate <- 0
	mhRate2 <- 0
	mhRate3 <- 0
	myPositions <- floor(seq(0, nTime, length = 7))  #c(0,50,100,150,200,250,300)
	myLabels <- round(myPositions*timeScale,1)
	uchar <- myUnicodeCharacters()
	for(iter in 2:nIter){
		cutPoints[iter, ] <- cutPoints[iter - 1, ]
		lValues[iter] <- lValues[iter - 1]
		chooseMove <- sample(movesRange, 1)
		if(chooseMove == "1"){
#			Move a: proposal from prior
			newCutPoints <- simulateFromPrior(nTime = nTime, startPoint = startPoint)
			propLogL <- logLikelihood(myData = myData, cutPoints = newCutPoints,
					sdPerPoint = sdPerPoint, modelVariance = modelVariance, startPoint = startPoint)
			mhRatio <- propLogL - lValues[iter - 1]
			u <- runif(1)
			if(log(u) < mhRatio){
				cutPoints[iter, ] <- newCutPoints
				lValues[iter] <- propLogL
				mhRate <- mhRate + 3
			}
		}
		if(chooseMove == "2"){
#			Move b: local random walk
			lpProp <- localProposal(cutPoints = cutPoints[iter, ] , nTime = nTime, mhPropRange = mhPropRange, startPoint = startPoint)
			newCutPoints <- lpProp$newState
			propLogL <- logLikelihood(myData = myData, cutPoints = newCutPoints,
					sdPerPoint = sdPerPoint, modelVariance = modelVariance, startPoint = startPoint)
			mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
				logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
				logPrior(cutPoints = cutPoints[iter, ], nTime = nTime, startPoint = startPoint)

			u <- runif(1)
			if(log(u) < mhRatio){
				cutPoints[iter, ] <- newCutPoints
				lValues[iter] <- propLogL
				mhRate2 <- mhRate2 + 3
			}
		}
		if(chooseMove == "3"){
#			Move c: random walk to just one index
			lpProp <- singleLocalProposal(cutPoints = cutPoints[iter, ], nTime = nTime, 
					mhSinglePropRange = mhSinglePropRange, startPoint = startPoint)
			newCutPoints <- lpProp$newState
			propLogL <- logLikelihood(myData = myData, cutPoints = newCutPoints,
						sdPerPoint = sdPerPoint, modelVariance = modelVariance, startPoint = startPoint)
			mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
				logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
				logPrior(cutPoints = cutPoints[iter, ], nTime = nTime, startPoint = startPoint)

			u <- runif(1)
			if(log(u) < mhRatio){
				cutPoints[iter, ] <- newCutPoints
				lValues[iter] <- propLogL
				mhRate3 <- mhRate3 + 3
			}
		}
		lPosterior[iter] <- lValues[iter] + logPrior(cutPoints = cutPoints[iter, ], nTime = nTime, startPoint = startPoint)
		if(showProgress){
		if(iter %% 1000 == 0){
			par(mfrow = c(1,3))
			plot(lPosterior[1:iter], type= "l", xlab = "mcmc iteration", ylab = "log-Posterior")
			matplot(cutPoints[1:iter,]*timeScale, type = "l", col = 2:4, xlab = "mcmc iteration", lty = 1, ylab = "hours")
			matplot(myData, type = "l", col = 1, xaxt = "n", xlab = "hours", ylab = "growth")
			axis(1, at = myPositions, labels = myLabels)
			abline(v = cutPoints[iter, ], col = c(2,3,4))
			abline(v = startPoint, lty = 4, col = "gray")
			legend("bottomright", paste0("replicate ",1:3), lty = 1:3)
			text(x = cutPoints[iter, ] - 20, y = 0.8, c("phase 1", "phase 2", "phase 3"), col = c(2,3,4))
			if(missing(dName)  == FALSE){
				title( dName, outer = TRUE , line=-1)
			}
			cat(paste0("  Iteration: ", iter, ". MH acceptance rates: (M1) ", 
				round(mhRate*100/iter,2), "%, (M2) ", 
				round(mhRate2*100/iter,2) ,"%, (M3) ", round(mhRate3*100/iter,2), "%."),"\n")
		}
		}
		else{
			if(iter %% (nIter/2) == 0){cat(paste0(uchar, " "))}
		}
	}
	if(showProgress == FALSE){
		cat(paste0(" Accepted MH moves: [Move 1] ", round(mhRate*100/iter,2), "%, [Move 2] ",round(mhRate2*100/iter,2) ,"%, [Move 3] ", round(mhRate3*100/iter,2), "%."),"\n")
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
		abline(v = apply(cutPoints, 2, median), col = c(2,3,4), lwd = 2)
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
	sdPerPoint <- numeric(nTime)
	j <- 0
	for(i in 1:n){
		myData <- rbind(myDataList[[1]][,i], myDataList[[2]][,i], myDataList[[3]][,i])
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


growthPhaseMCMC <- function(myDataList, burn, nIter, mhPropRange, mhSinglePropRange, movesRange, startPoint, getSDvalues, timeScale, blankThreshold, savePlots, showProgress, zeroNormalization){
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
		cat(paste0("*  Estimating variances... "))
		sdValues <- getVariance(myDataList = myDataList)
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
			mhRunForSubject <- mhSampler(myData = myData, nIter = nIter, mhPropRange = mhPropRange, dName = myColNames[i], burn = burn,
							mhSinglePropRange = mhSinglePropRange, movesRange = movesRange, finalIterationPdf = savePlots,
							startPoint = startPoint, sdValues = sdValues, timeScale = timeScale, showProgress = showProgress)
			cutPoints[i, ] <- apply(mhRunForSubject$cutPoints[-(1:burn), ], 2, median)
			cutPointsVar[i, ] <- apply(mhRunForSubject$cutPoints[-(1:burn), ], 2, var)
			getArea <- areaPerPhase(cutPoints = mhRunForSubject$cutPoints[-(1:burn), ], myData = myData, timeScale = timeScale)$PosteriorSummary
			getRate <- ratePerPhase(cutPoints = mhRunForSubject$cutPoints[-(1:burn), ], myData = myData, timeScale = timeScale)$PosteriorSummary
			rateMeanPerPhase[i, ] <- getRate[,1]
			rateVarPerPhase[i, ] <- getRate[,2]
			areaMeanPerPhase[i, ] <- getArea[,1]
			areaVarPerPhase[i, ] <- getArea[,2]
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







