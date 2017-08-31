#library("mvoutlier")
library("truncnorm")

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



replaceIndexWithMean <- function(mu, nu, alpha, beta, myIndex){
	m <- length(mu)
	sigma2 <- beta/(alpha - 1)
	mySd <- sqrt( sigma2/nu )
	simMeans <- rnorm(n = m, mean = 0, sd = 1)
	simMeans <- mu + simMeans*mySd
	simMeans[myIndex] <- mu[myIndex]
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
#	results[[1]] <- rnorm(m, mean = thetaOld, sd = tau)
	results[[1]] <- rnorm(m, mean = 0, sd = 1)
	results[[1]] <- thetaOld + tau*results[[1]]
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
	L <- length(cutPoints)
	augmentedCutPoints <- c(startPoint - 1, cutPoints, nTime )
	sanityCheck <- diff(augmentedCutPoints)
	if(  min(sanityCheck) <= 0 ){
		logP <- log(0)
	}else{
		logP <- 0
		if(L > 0){
			logP <- -log(nTime - L - startPoint + 1)
			if(L > 1){
				for(l in 2:L){
					logP <- logP - log(nTime - L + l - 1 - cutPoints[l-1]) 
				}
			}
		}
	}
	return(logP)
}

simulateFromPrior <- function(nTime, startPoint, L = 3){
	cutPoints <- numeric(L)
	indexRange <- startPoint:(nTime-L)
	if(L > 0){
		if(length(indexRange) == 1){
			cutPoints[1] <- indexRange}else{
			cutPoints[1] <- sample(indexRange, 1)
		}
		if(L > 1){
			for(l in 2:L){
				indexRange <- (cutPoints[l-1] + 1):(nTime - L + l - 1)
				if(length(indexRange) == 1){
					cutPoints[l] <- indexRange}
				else{			
					cutPoints[l] <- sample(indexRange, 1)
				}
			}
		}
	}
	return(cutPoints)
}

localProposal <- function(cutPoints, nTime, mhPropRange, startPoint){
	L <- length(cutPoints)
	if(L < 1){stop("L = 0.")}
	if(missing(mhPropRange)){mhPropRange = 1}
	epsilon <- sample(c(-seq(1:mhPropRange), 0, seq(1:mhPropRange)),L, replace = TRUE)
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
	L <- length(cutPoints)
	if(L < 1){stop("L = 0.")}
	if(missing(mhSinglePropRange)){mhSinglePropRange = 1}
	epsilon <- sample(c(-seq(1:mhSinglePropRange), seq(1:mhSinglePropRange)),1)
	justOneIndex <- sample(L,1)
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


birthProbs <- function(LRange){
	Lmin = 0
	Lmax <- max(LRange)
	probs <- numeric(Lmax + 1)
	probs[1] = 1
	probs[Lmax + 1] = 0
	if(Lmax - Lmin + 1 > 2) {
		probs[(Lmin+2):(Lmax)] = 0.5
	}
	names(probs) <- 0:Lmax
	return(probs)
}



truncatedPoisson <- function(Lmax = 20, gammaParameter){
	# default: gammaParameter = 6
        logPrior <- log(numeric(Lmax))
	Lmin = 0
        logPrior1 <- exp(-(0:Lmax)*gammaParameter)
        logPrior1 <- logPrior1/sum(logPrior1)
        logPrior <- logPrior1
	names(logPrior) <- 0:Lmax
	return(logPrior)

}


complexityPrior <- function(Lmax = 20, gammaParameter, nTime){
	Lmin = 0
        logPrior1 <- exp(-(0:Lmax)*gammaParameter*log(3*nTime/(0:Lmax)))
	logPrior1[1] <- 1
        logPrior1 <- logPrior1/sum(logPrior1)
        logPrior <- log(logPrior1)
	names(logPrior) <- 0:Lmax
	return(logPrior)

}


priorNumberOfCutpoints <- function(epsilon, LRange, secondComponent){
	if(epsilon > 1){stop("epsilon should be less than 1.")}
	if(epsilon < 0){stop("epsilon should be larger than 0.")}
	
	Lmax <- max(LRange)
	Lmin <- 0
	# truncatedPoisson (this is not logged)
        logPrior <- secondComponent

	diracPrior <- rep(0,Lmax+1)
	diracPrior[1] <- 1
	diracPrior <- diracPrior/sum(diracPrior)
	finalPrior <- (1-epsilon)*diracPrior + epsilon*logPrior
	logPrior <- log(finalPrior)
	names(logPrior) <- 0:Lmax
	return(logPrior)
	
}



updateEpsilon <- function(epsilon_current, beta_prior_parameters, proposal_sd, secondComponent, L, Lmax){
	# proposal via the truncated normal distribution in the (0,1) interval
	epsilon_new <- rtruncnorm(n = 1, a = 0, b = 1, mean = epsilon_current, sd = proposal_sd)
	dens_epsilon_new <- log(dtruncnorm(x = epsilon_new, a = 0, b = 1, mean = epsilon_current, sd = proposal_sd))
	dens_epsilon_current <- log(dtruncnorm(x = epsilon_current, a = 0, b = 1, mean = epsilon_new, sd = proposal_sd))
	prop_dens_ratio <- dens_epsilon_current - dens_epsilon_new
	target_dens_new <- (beta_prior_parameters[1] - 1)*log(epsilon_new) + (beta_prior_parameters[2] - 1)*log(1-epsilon_new) + 
				priorNumberOfCutpoints(epsilon = epsilon_new, LRange = 0:Lmax, secondComponent = secondComponent)[L+1]
	target_dens_current <- (beta_prior_parameters[1] - 1)*log(epsilon_current) + (beta_prior_parameters[2] - 1)*log(1-epsilon_current) + 
				priorNumberOfCutpoints(epsilon = epsilon_current, LRange = 0:Lmax, secondComponent = secondComponent)[L+1]
	target_dens_ratio <- target_dens_new - target_dens_current
	mh_ratio <- target_dens_ratio + prop_dens_ratio
	myResult <- vector("list", length = 2)
	myResult[[1]] <- mh_ratio
	myResult[[2]] <- epsilon_new
	names(myResult) <- c("mh_ratio", "epsilon_new")
	return(myResult)
}




updateNumberOfCutpoints <- function(cutPoints, nTime, startPoint, LRange, birthProbs){
	results <- vector("list", length = 4)
	L <- length(cutPoints)
	augmentedCutPoints <- c(startPoint - 1, cutPoints, nTime )
	#choose birth or death according to birthProbs
	u = runif(1)
	if( u < birthProbs[as.character(L)]){
		#	birth move
		nIntervals <- L + 1
		intervalIndex <- sample(L+1, 1)	
		myIndex <- (augmentedCutPoints[intervalIndex]+1):(augmentedCutPoints[intervalIndex+1]-1)
		lengthInterval <- augmentedCutPoints[intervalIndex+1] - augmentedCutPoints[intervalIndex] - 1
		if( lengthInterval < 1 ){
			propRatio <- log(0)
			moveType = "notValidBirth"
			newState = cutPoints
		}else{
			moveType = "birth"	
			if(lengthInterval > 1){
				newCutPoint <- sample(myIndex, 1)
			}else{
				newCutPoint <- myIndex
			}
			results[[4]] <- newCutPoint
			newState <- numeric(L+1)
			newState[intervalIndex] <- newCutPoint
			if( intervalIndex > 1 ){
				newState[1:(intervalIndex-1)] <- cutPoints[1:(intervalIndex-1)]
			}
			if(intervalIndex < L+1){
				newState[(intervalIndex+1):(L+1)] <- cutPoints[intervalIndex:L]  
			}
			propRatio <- log(lengthInterval) + log(1 - birthProbs[as.character(L+1)]) - log(birthProbs[as.character(L)])
		}
	}else{
		moveType = "death"
		#	death move
		myIndex <- sample(L, 1)
		results[[4]] <- cutPoints[myIndex]
		newState <- cutPoints[-myIndex]
		lengthInterval <- augmentedCutPoints[myIndex+2] - augmentedCutPoints[myIndex] - 1
		propRatio <- - log(lengthInterval) - log(1 - birthProbs[as.character(L)]) + log(birthProbs[as.character(L-1)])
	}

	results[[1]] <- newState
	results[[2]] <- propRatio
	results[[3]] <- moveType
	names(results) <- c("newState", "propRatio", "moveType", "timePoint")
	return(results)
}


mcmcSampler <- function(myData, nIter, finalIterationPdf, modelVariance, mhPropRange, mhSinglePropRange, movesRange, startPoint, 
		postPar, dName, timeScale, burn, showProgress, iterPerPlotPrefix, priorParameters, L = 3, LRange, tau, 
		gammaParameter, saveTheta){
	# tau is the sd of the MH proposal
	if(missing(showProgress)){showProgress = FALSE}
	if(missing(timeScale)){timeScale = 1}
	if(missing(mhPropRange)){mhPropRange = 1}
	if (missing(modelVariance)){modelVariance = FALSE}
	if(missing(startPoint)){ startPoint = 2 }
	if(startPoint < 2){ startPoint = 2 }
	if(missing(burn)){burn = 1}
	if(missing(movesRange)){movesRange = as.character(1:3)}
	if (missing(finalIterationPdf)){ finalIterationPdf = NULL }
	if (missing(LRange)){LRange = L}
	Lmax <- max(LRange)
#	truncPois <- truncatedPoisson(Lmax = Lmax, gammaParameter = gammaParameter)
	nTime <- dim(myData)[1]
	truncPois <- complexityPrior(Lmax = Lmax, gammaParameter = gammaParameter, nTime = nTime)
	logPriorNPoints <- truncPois
	cutPoints <- array(data = NA, dim = c(nIter, Lmax))
	if(saveTheta == TRUE){
		thetaValues <- array(data = NA, dim = c(nIter, nTime))
	}
	Lvalues <- numeric(nIter)
	mu <- array(data = NA, dim = c(nIter, nTime))
	sigma2 <- array(data = NA, dim = c(nIter, nTime))
	iter <- 1
#	L <- Lvalues[iter] <- floor(Lmax/2)
	L <- Lvalues[iter] <- 1
	cutPoints[iter, 1:L] <- startPoint + (1:L)*floor((nTime - startPoint)/(L+1))

#	if( L == 3){
#		cutPoints[iter, 1:L] <- c(startPoint + 1, startPoint + 10, startPoint + 30)
#	}else{
#		cutPoints[iter, 1] <- startPoint + 1
#		if(L>1){
#			for(l in 2:L){
#				cutPoints[iter, l] <- startPoint + 10*(l-1)
#			}
#		}
#	}

	myBirthProbs <- birthProbs(LRange)

	overkill <- simMultiIndNormInvGamma(mu = postPar[[1]], nu = postPar[[2]], alpha = postPar[[3]], beta = postPar[[4]])
	theta <- overkill
	mu[iter, ] <- overkill$mean
	sigma2[iter, ] <- overkill$var
	lValues <- lPosterior <- numeric(nIter)
	lValues[iter] <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, 0:L], theta = overkill, startPoint = startPoint)
	lPosterior[iter] <- lValues[iter] + logPrior(cutPoints = cutPoints[iter, 0:L], nTime = nTime, startPoint = startPoint) 
	mhRate0a <- 0
	mhRate0 <- 0
	mhRate0c <- 0
	mhRate <- 0
	mhRate2 <- 0
	mhRate3 <- 0
	mhRatePoints <- 0
	arEpsilon <- 0
	myPositions <- floor(seq(0, nTime, length = 10))  #c(0,50,100,150,200,250,300)
	myLabels <- round(myPositions*timeScale,1)
	uchar <- myUnicodeCharacters()
#	epsilonValues <- numeric(nIter)
#	epsilonValues[1] <- rbeta(n = 1, shape1 = beta_prior_parameters[1], shape2 = beta_prior_parameters[2])
	sampleSD <- sqrt(postPar[[4]]/(postPar[[3]] - 1))
	#plot(sampleSD)
	myFl <- floor((nIter/4))
	if(showProgress){pdf(file = paste0(finalIterationPdf,"/",dName,".pdf"),width = 9, height = 6)}
	for(iter in 2:nIter){
		lValues[iter] <- lValues[iter - 1]
		cutPoints[iter, ] <- cutPoints[iter - 1, ]
#-------------------------------------------------------------------------------------------------------------------------------------
#		update number of cutpoints		
#		L <- sample(LRange, 1)
		lValues[iter] <- lValues[iter - 1]
		if(L > 0){
			cutPoints[iter, 1:L] = cutPoints[iter - 1, 1:L]
		}
		myProposal <- updateNumberOfCutpoints(cutPoints = cutPoints[iter - 1, 0:L], nTime = nTime, 
				startPoint = startPoint, LRange = LRange, birthProbs = myBirthProbs)
#		overkillProposed <- replaceIndexWithMean(mu = postPar[[1]], nu = postPar[[2]], alpha = postPar[[3]], beta = postPar[[4]], myIndex = myProposal$timePoint)
		overkillProposed <- theta
		propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = myProposal$newState, theta = overkillProposed, startPoint = startPoint)
		if(myProposal$moveType != "notValidBirth"){
			if(myProposal$moveType == "birth"){Lprop = L + 1}else{Lprop = L - 1}
#			logPriorNPoints <- priorNumberOfCutpoints(epsilon = epsilonValues[iter - 1], LRange = LRange, secondComponent = truncPois)
			mhRatio <- propLogL + logPriorNPoints[as.character(Lprop)] - lValues[iter] - logPriorNPoints[as.character(L)] + myProposal$propRatio
			if(L > 0){
				augCutPoints <- c(1, cutPoints[iter, 1:L], nTime )
			}else{augCutPoints <- c(1, nTime )}
			augCutPointsProp <- c(1, myProposal$newState, nTime )
			#zCurrent <- (theta$mean[augCutPoints] - priorParameters$mu0[augCutPoints])/sqrt(theta$variance[augCutPoints]/priorParameters$nu0)
			#zProp <- (overkillProposed$mean[augCutPointsProp] - priorParameters$mu0[augCutPointsProp])/sqrt(overkillProposed$variance[augCutPointsProp]/priorParameters$nu0)
			#priorRatio <- sum(dnorm( x = zProp, mean = 0, sd = 1, log = TRUE)) - sum(dnorm( x = zCurrent, mean = 0, sd = 1, log = TRUE))
			priorRatio <- 0
			mhRatio <- mhRatio + priorRatio
			if( log(runif(1)) < mhRatio ){
				lValues[iter] <- propLogL
				mhRatePoints <- mhRatePoints + 1
				theta$mean <- overkillProposed$mean
				L = Lprop
				cutPoints[iter, ] <- rep(0, Lmax)
				if(L > 0){
					cutPoints[iter, 1:L] <- myProposal$newState
				}
				lValues[iter] <- propLogL
			}
#			cutPoints[iter, 1:L] <- simulateFromPrior(nTime = nTime, startPoint = startPoint, L = L)
		}
		Lvalues[iter] <- L
#-------------------------------------------------------------------------------------------------------------------------------------
#		Move 0.a: update theta parameters using the standard posterior as proposal
	#	overkill <- simMultiIndNormInvGamma(mu = postPar[[1]], nu = postPar[[2]], alpha = postPar[[3]], beta = postPar[[4]])
	#	propCurrent <- partialLogLikelihoodStandardModel(myData = myData, theta = theta, cutPoints = cutPoints[iter, 0:L])
	#	propNew <- partialLogLikelihoodStandardModel(myData = myData, theta = overkill, cutPoints = cutPoints[iter, 0:L])
	#	propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, 0:L], theta = overkill, startPoint = startPoint)
	#	mhRatio <- propLogL + propCurrent - lValues[iter] - propCurrent
	#	if( log(runif(1)) < mhRatio ){
	#		lValues[iter] <- propLogL
	#		mhRate0a <- mhRate0a + 1
	#		theta$mean <- overkill$mean
	#		lValues[iter] <- propLogL
	#	}


#		Move 0.c: update all theta parameters using a standard symmetric proposal based on the current values	
		overkill <- proposeTheta(thetaOld = theta$mean, tau = tau*sampleSD,  alpha = postPar[[3]], beta = postPar[[4]])
		propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, 0:L], theta = overkill, startPoint = startPoint)
		mhRatio <- propLogL - lValues[iter]
		augCutPoints <- c(1, cutPoints[iter, 0:L], nTime )
#		zCurrent <- (theta$mean[augCutPoints] - priorParameters$mu0[augCutPoints])/sqrt(theta$variance[augCutPoints]/priorParameters$nu0)
#		zProp <- (overkill$mean[augCutPoints] - priorParameters$mu0[augCutPoints])/sqrt(overkill$variance[augCutPoints]/priorParameters$nu0)
		zCurrent <- (theta$mean - priorParameters$mu0)/sqrt(theta$variance/priorParameters$nu0)
		zProp <- (overkill$mean - priorParameters$mu0)/sqrt(overkill$variance/priorParameters$nu0)
		priorRatio <- sum(dnorm( x = zProp, mean = 0, sd = 1, log = TRUE)) - sum(dnorm( x = zCurrent, mean = 0, sd = 1, log = TRUE))
		mhRatio <- mhRatio + priorRatio
		if( log(runif(1)) < mhRatio ){
			lValues[iter] <- propLogL
			mhRate0c <- mhRate0c + 1
			theta <- overkill
		}
###########################################################################################################
#		Gibss move: update all theta[-augmentedCutPoints] from the prior.
		if(iter > burn){
			overkill <-  simMultiIndNormInvGamma(mu = priorParameters$mu0, 
								nu = priorParameters$nu0, 
								alpha = postPar[[3]], 
								beta = postPar[[4]]
							)
			if(L < nTime - 1){
				theta$mean[-augCutPoints] <- overkill$mean[-augCutPoints]
			}
		}
###########################################################################################################
		chooseMove <- sample(movesRange, 1)
		if(chooseMove == "1"){
			if(L > 0){
#			Move a: proposal from prior
				newCutPoints <- simulateFromPrior(nTime = nTime, startPoint = startPoint, L = L)
				propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
				mhRatio <- propLogL - lValues[iter - 1]
				u <- runif(1)
				if(log(u) < mhRatio){
					cutPoints[iter, 1:L] <- newCutPoints
					lValues[iter] <- propLogL
					mhRate <- mhRate + length(movesRange)
				}
			}
		}
		if(chooseMove == "2"){
#			Move b: local random walk
			if(L > 0){
				lpProp <- localProposal(cutPoints = cutPoints[iter, 1:L] , nTime = nTime, mhPropRange = mhPropRange, startPoint = startPoint)
				newCutPoints <- lpProp$newState
				propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
				mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
					logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
					logPrior(cutPoints = cutPoints[iter, 1:L], nTime = nTime, startPoint = startPoint)

				u <- runif(1)
				if(log(u) < mhRatio){
					cutPoints[iter, 1:L] <- newCutPoints
					lValues[iter] <- propLogL
					mhRate2 <- mhRate2 + length(movesRange)
				}
			}
		}
		if(chooseMove == "3"){
#			Move c: random walk to just one index
			if( L > 0 ){
				lpProp <- singleLocalProposal(cutPoints = cutPoints[iter, 1:L], nTime = nTime, 
						mhSinglePropRange = mhSinglePropRange, startPoint = startPoint)
				newCutPoints <- lpProp$newState
				propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
				mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
					logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
					logPrior(cutPoints = cutPoints[iter, 1:L], nTime = nTime, startPoint = startPoint)

				u <- runif(1)
				if(log(u) < mhRatio){
					cutPoints[iter, 1:L] <- newCutPoints
					lValues[iter] <- propLogL
					mhRate3 <- mhRate3 + length(movesRange)
				}
			}
		}
		if(chooseMove == "4"){
#			Move c: random walk to just one index
			if(L > 0){
				lpProp <- singleLocalProposal(cutPoints = cutPoints[iter, 1:L], nTime = nTime, 
						mhSinglePropRange = 2, startPoint = startPoint)
				newCutPoints <- lpProp$newState
				propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
				mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
					logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
					logPrior(cutPoints = cutPoints[iter, 1:L], nTime = nTime, startPoint = startPoint)

				u <- runif(1)
				if(log(u) < mhRatio){
					cutPoints[iter, 1:L] <- newCutPoints
					lValues[iter] <- propLogL
					mhRate3 <- mhRate3 + length(movesRange)
				}
			}
		}
		lPosterior[iter] <- lValues[iter] + logPrior(cutPoints = cutPoints[iter, 0:L], nTime = nTime, startPoint = startPoint)
		if(saveTheta == TRUE){
			thetaValues[iter, ] <- theta$mean
		}

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
			xPoints <- c(1, cutPoints[iter,1:L], nTime )
			yPoints <- theta$mean[xPoints]
			for(l in 1:L){
				points(xPoints[l:(l+1)], yPoints[l:(l+1)], type = "l", col = l+1, lwd = 4)
			}
			points(xPoints[(L+1):(L+2)], yPoints[(L+1):(L+2)], type = "l", col = "gray", lwd = 4)
			matplot(myData, type = "l", col = 1, add = TRUE)
			axis(1, at = myPositions, labels = myLabels)
			#abline(v = cutPoints[iter, ], col = c(2,3,4))
			#abline(v = startPoint, lty = 4, col = "gray")
			legend("bottomright", paste0("replicate ",1:3), lty = 1:3)
			text(x = cutPoints[iter, 1:L] - 20, y = 0.8, paste0("phase ", 1:L), col = 2:(L+1))
			if(missing(dName)  == FALSE){
				title( paste0(dName,", iter = ", iter), outer = TRUE , line=-1)
			}
			cat(paste0("  Iteration: ", iter, ". MH acceptance rates: (M1a) ", 
				round(mhRate0a*100/iter,2), "%, (M1b) ",
				round(mhRate0*100/iter,2), "%, (M1c) ",
				round(mhRate0c*100/iter,2), "%, (M2) ", 
				round(mhRate*100/iter,2), "%, (M3) ", 
				round(mhRate2*100/iter,2) ,"%, (M4) ", 
				round(mhRate3*100/iter,2), "%, (L) ",
				round(mhRatePoints*100/iter,2),"%, ",
				"L = ", L, ", epsilon = ", epsilonValues[iter]),"\n")
		}
		}
		else{
			if(iter %% myFl == 0){cat(paste0(uchar, " "))}
		}
	}
	acceptanceRates <- c( mhRatePoints*100/iter,  mhRate0c*100/iter, mhRate0a*100/iter, mhRate2*100/iter, mhRate3*100/iter )
	names(acceptanceRates) <- c("move_l_n", "move_theta_randomWalk", "move_theta_post", "move_tau_all", "move_tau_j")
	if(showProgress == FALSE){
		cat(paste0(" Accepted MH moves: [l_n] ", round(mhRatePoints*100/iter,2), "%, [","\U03B8","+","\U03B5",'] ',round(mhRate0c*100/iter,2) ,"%, [","\U03B8","|x] ", round(mhRate0a*100/iter,2),"%, [","\U03C4","] ", round(mhRate2*100/iter,2), "%, [","\U03C4","_j] ", round(mhRate3*100/iter,2), "%."),"\n")
	}else{
		dev.off()
	}
	results <- vector("list", length = 7)
	results[[2]] <- lValues
	results[[3]] <- acceptanceRates
	results[[4]] <- Lvalues
	LVtab <- table(Lvalues[-(1:burn)])/(length(Lvalues[-(1:burn)]))
	results[[6]] <- LVtab
	names(results) <- c("cutPoints", "logLikelihood", "acceptanceRates", "nCutPoints", "nCutPointsMAP", "nCutPointsPosterior", "theta")
	#posterior mode of number of cutpoints
	L <- as.numeric(names(table(Lvalues[-(1:burn)]))[order(table(Lvalues[-(1:burn)]), decreasing=T)[1]])
	cat(paste0( "     Most probable number of cutpoints equals ", L, " with: P(l_n = ",L,"|x) = ",as.numeric(-sort(-LVtab)[1])),"\n")
	results[[5]] <- L
	Lindex <- which(Lvalues == L)
	LindexFull <- which(Lvalues == L)
	burnIndex <- which(Lindex < burn + 1) 
	if(length(burnIndex) > 0){
		Lindex <- Lindex[-burnIndex]
	}

	results[[1]] = cutPoints[Lindex,1:L]
	if(saveTheta == TRUE){
		results[[7]] <- thetaValues[Lindex, ]
	}
	if(is.null(finalIterationPdf) == FALSE){
		pdf(file = paste0(finalIterationPdf,"/",dName,"_mcmcTrace.pdf"), width = 20, height = 10)
		split.screen( figs = c( 2, 1 ) )
		split.screen( figs = c( 1, 4 ), screen = 1 )
		split.screen( figs = c( 1, 1 ), screen = 2 )
		screen( 3 )
		plot(Lvalues, type = "l", xlab = "mcmc iteration", ylab = "number of cutpoints")
		#print('OK1')
		screen( 4 )
		barplot(table(Lvalues[-(1:burn)]), xlab = "number of cutpoints", ylab = "frequency")
		#print('OK2')
		screen( 5 )
		plot(lPosterior, ylim = range(lPosterior[-(1:burn)]), type= "l", xlab = "mcmc iteration", ylab = "log-posterior pdf")
		points(lPosterior[1:burn], type = "l", lty = 5, col = "white")
		legend("bottomright",lty = 3, col = 1, "burn-in period")
		#print('OK3')
		screen( 6 )
		if(L > 0){
			ylim = range(cutPoints[LindexFull,1:L]*timeScale)
			ylim[1] = 0
			matplot(cutPoints[LindexFull,1:L]*timeScale, ylim = ylim, type = "l", col = 2:(L+1), xlab = "mcmc iteration", lty = 1, ylab = "hours")
			matplot(cutPoints[Lindex,1:L]*timeScale, ylim = ylim, type = "l", col = "white", lty = 5, add = TRUE)
			legend("bottomright", paste0("cut-point of phase ", 1:L), lty = 1, col = 2:(L+1), bty = "n")
		}else{plot(c(1,1))}
		screen( 7 )
		matplot(myData, type = "l", col = 1, xaxt = "n", xlab = "hours", ylab = "growth")
		axis(1, at = myPositions, labels = myLabels)
		if(L > 0){
			if(L > 1){
				abline(v = apply(cutPoints[Lindex, 1:L], 2, median), col = 2:(L+1), lwd = 2)
			}else{
				abline(v = median(cutPoints[Lindex, 1:L]), col = 2:(L+1), lwd = 2)
			}
			legend("bottomright", paste0("replicate ",1:3), lty = 1:3, bty = "n")
			legend("topleft", paste0("posterior median ",1:L), lty = c(1,1,1), col = 2:(L+1), bty = "n")
			if(L > 1){mcmcVar <- apply(cutPoints[Lindex, 1:L],2,var)}else{
				mcmcVar <- var(cutPoints[Lindex, 1:L])
			}
			for(l in 1:L){
				if(mcmcVar[l] > 0){
					points(density(cutPoints[Lindex,l], adjust = 7), type = "l", col = l+1, lwd = 4)
					points(density(cutPoints[Lindex,l], adjust = 7), type = "l", col = "gray", lwd = 2)
				}
			}
			text(x = cutPoints[iter,1:L ] - 20, y = 0.8, paste0("phase ", 1:L), col = 2:(L+1))
		}
		if(missing(dName)  == FALSE){
			title( dName, outer = TRUE , line=-1)
		}
		close.screen( all = TRUE )
		dev.off()
	}

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


growthPhaseFullMCMC <- function(myDataList, burn, nIter, mhPropRange = 1, mhSinglePropRange, movesRange = c(2,3), startPoint = 2, 
		getSDvalues, timeScale, blankThreshold, savePlots, showProgress = FALSE, zeroNormalization = FALSE, L = 3, LRange = 0:30, tau = 0.05,
		gammaParameter = 6, nu0 = 1, alpha0 = 1, beta0 = 1, subsetIndex, saveTheta = FALSE, sameVariance = TRUE){
#	burn = 2000, nIter = 5000,mhPropRange = 1, mhSinglePropRange = 50, getSDvalues = T, startPoint=54, timeScale = 1/6,
	myColNames <- colnames(myDataList[[1]])
	if(missing(timeScale)){timeScale = 1/1}
	if(missing(showProgress)){showProgress = FALSE}
#	if(missing(blankThreshold)){blankThreshold = 0.02}
	if(missing(blankThreshold)){blankThreshold = -10^9}
	if(missing(burn)){burn = 2000}
	if(missing(nIter)){nIter = 5000}
	if(missing(zeroNormalization)){zeroNormalization = TRUE}
	if(missing(LRange)){LRange = L}
	Lmax <- max(LRange)
	LRange = 0:Lmax
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
	if(missing(startPoint)){ startPoint = 2 }
	if( startPoint < 2){ startPoint = 2 }
	if(missing(mhPropRange)){mhPropRange = 1}
	if(missing(mhSinglePropRange)){ mhSinglePropRange = dim(myDataList[[1]])[1]/10  }
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
		priorParameters <- computeEmpiricalPriorParameters(myDataList = myDataList, nu0 = nu0, alpha0 = alpha0, beta0 = beta0)
		if(sameVariance == TRUE){
			cat(paste0(" using the same variance per time series... "))
			posteriorParameters <- computePosteriorParameters(myDataList = myDataList, priorParameters = priorParameters)
		}else{
			cat(paste0(" using different variance per time series... "))
			posteriorParameters <- computePosteriorParametersFree(myDataList = myDataList, priorParameters = priorParameters)
		}
		cat(paste0(" done.","\n"))
	}
	zz <- file(paste0(savePlots,"/numberCutPointsMAP.txt"), open = "w")
	nReps <- length(myDataList)
	n <- dim(myDataList[[1]])[2]
	cutPoints <- vector("list", length = n)
	cutPointsVar <- vector("list", length = n)
#	cutPoints <- array(data = NA, dim = c(n, L))
#	cutPointsVar <- array(data = NA, dim = c(n, L))
	areaMeanPerPhase <- array(data=NA, dim = c(n, L))
	areaVarPerPhase <- array(data=NA, dim = c(n, L))
	rateMeanPerPhase <- array(data=NA, dim = c(n, L))
	rateVarPerPhase <- array(data=NA, dim = c(n, L))
	colnames(areaMeanPerPhase) <- colnames(areaVarPerPhase) <- paste0("phase_", 1:L)
	if(nReps < 2){stop("no replicates")}
	cat(paste0("*  MCMC sampler parameters: nIter = ", nIter, ", burn = ", burn, ", startPoint = ", startPoint ),"\n")
	cat(paste0("*  Running MCMC for ", n, " subjects..."), "\n")	
	NAindex <- c()
	results <- vector("list", length = 9)
	results[[3]] <- vector("list", length = n)
	results[[5]] <- vector("list", length = n)
	results[[4]] <- numeric(n)
	results[[6]] <- myColNames
	results[[7]] <- vector("list", length = n)
	results[[8]] <- vector("list", length = n)
	results[[9]] <- vector("list", length = n)
	checkCondition <- -99
	if( missing(subsetIndex) == FALSE ){
		checkCondition <- (1:n)[-subsetIndex]
		cat(paste0("*                [NOTE] skipping ", n - length(subsetIndex), " subjects"), "\n")	
	}
	for (i in 1:n){
		myData <- myDataList[[1]][ , i]
		if( i %in% checkCondition ){
			#cat(paste0("                [SKIPPED] subsetIndex enabled.","\n"))
			NAindex <- c(NAindex, i)
		}else{
			cat(paste0("*    i = ", i, ", name: ", myColNames[i], " " ))
			for (j in 2:nReps){
				myData <- cbind(myData, myDataList[[j]][ , i])
			}
			posteriorParametersIter <- posteriorParameters
			posteriorParametersIter[[1]] <- posteriorParameters[[1]][i,]
			if(sameVariance == FALSE){
				posteriorParametersIter[[4]] <- posteriorParameters[[4]][i,]
			}
			mhRunForSubject <- mcmcSampler(myData = myData, nIter = nIter, mhPropRange = mhPropRange, dName = myColNames[i], burn = burn,
							mhSinglePropRange = mhSinglePropRange, movesRange = movesRange, finalIterationPdf = savePlots,
							startPoint = startPoint, postPar = posteriorParametersIter, timeScale = timeScale, 
							showProgress = showProgress, iterPerPlotPrefix = paste0("plot-", i), priorParameters = priorParameters, 
							L = L, LRange = LRange, tau = tau, 
							gammaParameter = gammaParameter, saveTheta = saveTheta)
			results[[3]][[i]] <- mhRunForSubject$nCutPointsPosterior
			results[[4]][i] <- mhRunForSubject$nCutPointsMAP
			results[[5]][[i]] <- mhRunForSubject$acceptanceRates
			results[[7]][[i]] <- mhRunForSubject$cutPoints
			results[[8]][[i]] <- mhRunForSubject$theta
			results[[9]][[i]] <- mhRunForSubject$nCutPoints
			if(is.array(mhRunForSubject$cutPoints) == TRUE){
				cutPoints[[i]] <- apply(mhRunForSubject$cutPoints, 2, median)
				cutPointsVar[[i]] <- apply(mhRunForSubject$cutPoints, 2, var)
			}else{
				cutPoints[[i]] <- median(mhRunForSubject$cutPoints)
				cutPointsVar[[i]] <- var(mhRunForSubject$cutPoints)
			}
			cat(paste0(myColNames[i], " ", mhRunForSubject$nCutPointsMAP), file = zz, "\n")
		}
	}
	close(zz)
	results[[1]] <- lapply(cutPoints, function(x){x*timeScale})
	results[[2]] <- lapply(cutPointsVar, function(x){x*(timeScale^2)})
	names(results) <- c("Cutpoint_posterior_median", 
				"Cutpoint_posterior_variance", 
				"NumberOfCutPoints_posterior_distribution", 
				"NumberOfCutPoints_MAP", 
				"Metropolis-Hastings_acceptance_rate", 
				"subject_ID",
				"Cutpoint_mcmc_trace_map",
				"theta",
				"nCutPointsTrace")
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
	if(missing(blankThreshold)){blankThreshold = -10^9}
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



# different variance per time series
computePosteriorParametersFree <- function(myDataList, priorParameters, blankThreshold){
	if(missing(blankThreshold)){blankThreshold = -10^9}
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
	varianceParameters1 <- alpha0 + nReps/2	
	varianceParameters2 <- array(data = 0, dim = c(n, nTime))	
	for(i in alive){
		myData <- myDataList[[1]][,i]
		for(k in 2:nReps){
			myData <- rbind(myData, myDataList[[k]][,i])
		}
		cSum <- colSums(myData)
		meanParameters1[i, ] <- (cSum + nu0*mu0)/meanParameters2
		sumOverN <- meanParameters2*colSums(myData^2) - cSum^2 - 2*nu0*mu0*cSum
		varianceParameters2[i, ] <- beta0 + 0.5*(nReps*nu0*mu0^2 + sumOverN)/meanParameters2
	}

	results <- vector("list", length = 4)
	results[[1]] <- meanParameters1
	results[[2]] <- meanParameters2
	results[[3]] <- varianceParameters1
	results[[4]] <- varianceParameters2

	names(results) <- c("meanPar1", "meanPar2", "varPar1", "varPar2")
	return(results)
	
}






# ara prepei na valw opwsdipote tin empirical prior ston meso.


computeEmpiricalPriorParameters <- function(myDataList, blankThreshold, nu0 = 1, alpha0 = 1, beta0 = 1){
	priorParameters <- vector("list",length = 4)
	if(missing(blankThreshold)){blankThreshold = -10^9}
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
	priorParameters$nu0 <- nu0
	priorParameters$alpha0 <- alpha0
	priorParameters$beta0 <- beta0
	return(priorParameters)

}


