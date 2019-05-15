
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
		x_minus_mean2 <- logL
	}else{
		goldenBoy <- theta$mean[augmentedCutPoints]
		myMeans <- numeric(nTime)
		for (j in 2:length(augmentedCutPoints)){
			subIndex <- augmentedCutPoints[j-1]:augmentedCutPoints[j]
			myMeans[subIndex] <- (goldenBoy[j] - goldenBoy[j-1])/(augmentedCutPoints[j] - augmentedCutPoints[j-1])*(subIndex - augmentedCutPoints[j-1]) + goldenBoy[j-1]
		}
		myMeanMatrix <- matrix(myMeans, nrow = nTime, ncol = dim(myData)[2] )
		sdPerPoint <- sqrt( theta$variance )		
		mySDMatrix <- matrix(sdPerPoint, nrow = nTime, ncol = dim(myData)[2] )
		x_minus_mean2 <- myData - myMeanMatrix 
		zData <- x_minus_mean2/mySDMatrix
		logL <- sum( apply(zData, 1, function(y){ sum(dnorm(y, mean = 0, sd = 1, log = TRUE)) })   )
		x_minus_mean2 <- rowSums(x_minus_mean2^2)
	}
	myResult <- vector("list", length = 2)
	myResult[[1]] <- logL
	myResult[[2]] <- x_minus_mean2
	names(myResult) <- c("logL", "x_minus_mean2")
	return(myResult)	
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



truncatedPoisson <- function(Lmax = 20, gammaParameter = 1){
        logPrior <- log(numeric(Lmax))
        Lmin = 0
        logPrior1 <- dpois(0:Lmax, lambda=gammaParameter)
        logPrior1 <- logPrior1/sum(logPrior1)
        logPrior <- log(logPrior1)
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
		postPar, dName, timeScale, burn,  iterPerPlotPrefix, priorParameters, L = 3, LRange, tau, 
		gammaParameter, saveTheta, Prior = "complexity"){
	# tau is the sd of the MH proposal
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
	nTime <- dim(myData)[1]
	if(Prior == "complexity"){
		logPriorNPoints <- complexityPrior(Lmax = Lmax, gammaParameter = gammaParameter, nTime = nTime)
	}else{
		logPriorNPoints <- truncatedPoisson(Lmax = Lmax, gammaParameter = gammaParameter)
	}
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

	myBirthProbs <- birthProbs(LRange)

	overkill <- simMultiIndNormInvGamma(mu = postPar[[1]], nu = postPar[[2]], alpha = postPar[[3]], beta = postPar[[4]])
	theta <- overkill
	mu[iter, ] <- overkill$mean
	sigma2[iter, ] <- overkill$variance
	lValues <- lPosterior <- numeric(nIter)
	lValues[iter] <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, 0:L], theta = overkill, startPoint = startPoint)$logL
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
		overkillProposed <- theta
		propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = myProposal$newState, theta = overkillProposed, startPoint = startPoint)$logL
		if(myProposal$moveType != "notValidBirth"){
			if(myProposal$moveType == "birth"){Lprop = L + 1}else{Lprop = L - 1}
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

#		Move 0.c: update all theta parameters using a standard symmetric proposal based on the current values	
		overkill <- proposeTheta(thetaOld = theta$mean, tau = tau*sampleSD,  alpha = postPar[[3]], beta = postPar[[4]])
		propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, 0:L], theta = overkill, startPoint = startPoint)$logL
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
		if(Prior == "complexity"){
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
		}else{
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
				propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)$logL
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
				propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)$logL
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
				propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)$logL
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
				propLogL <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)$logL
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
		if(iter %% myFl == 0){cat(paste0(uchar, " "))}
	}
	acceptanceRates <- c( mhRatePoints*100/iter,  mhRate0c*100/iter, mhRate2*100/iter, mhRate3*100/iter )
	names(acceptanceRates) <- c("move_l_n", "move_theta_randomWalk", "move_tau_all", "move_tau_j")
	cat(paste0(" Accepted MH moves: [l_n] ", round(mhRatePoints*100/iter,2), "%, [","\U03B8","+","\U03B5",'] ',round(mhRate0c*100/iter,2) ,"%, [","\U03C4","] ", round(mhRate2*100/iter,2), "%, [","\U03C4","_j] ", round(mhRate3*100/iter,2), "%."),"\n")
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
		close.screen( all.screens = TRUE )
		dev.off()
	}

	return(results)
}



mcmcSamplerUnknownVariance <- function(myData, nIter, finalIterationPdf, modelVariance, mhPropRange, mhSinglePropRange, movesRange, startPoint, 
		postPar, dName, timeScale, burn,  iterPerPlotPrefix, priorParameters, L = 3, LRange, tau, 
		gammaParameter, saveTheta, Prior = "complexity"){
	# tau is the sd of the MH proposal
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
	nTime <- dim(myData)[1]
	nReps <- dim(myData)[2]
	if(Prior == "complexity"){
		logPriorNPoints <- complexityPrior(Lmax = Lmax, gammaParameter = gammaParameter, nTime = nTime)
	}else{
		logPriorNPoints <- truncatedPoisson(Lmax = Lmax, gammaParameter = gammaParameter)
	}
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

	myBirthProbs <- birthProbs(LRange)

	overkill <- simMultiIndNormInvGamma(mu = postPar[[1]], nu = postPar[[2]], alpha = postPar[[3]], beta = postPar[[4]])
	theta <- overkill
	mu[iter, ] <- overkill$mean
	sigma2[iter, ] <- overkill$variance
	lValues <- lPosterior <- numeric(nIter)
	voivod <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, 0:L], theta = overkill, startPoint = startPoint)
	lValues[iter] <- voivod$logL
	x_minus_mean2 <- voivod$x_minus_mean2
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
		overkillProposed <- theta
		voivod <- logLikelihoodFullModel(myData = myData, cutPoints = myProposal$newState, theta = overkillProposed, startPoint = startPoint)
		propLogL <- voivod$logL
		if(myProposal$moveType != "notValidBirth"){
			if(myProposal$moveType == "birth"){Lprop = L + 1}else{Lprop = L - 1}
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
				x_minus_mean2 <- voivod$x_minus_mean2
			}
#			cutPoints[iter, 1:L] <- simulateFromPrior(nTime = nTime, startPoint = startPoint, L = L)
		}
		Lvalues[iter] <- L
#-------------------------------------------------------------------------------------------------------------------------------------

#		Move 0.c: update all theta parameters using a standard symmetric proposal based on the current values	
		overkill <- proposeTheta(thetaOld = theta$mean, tau = tau*sampleSD,  alpha = 2, beta = sigma2[iter - 1, ])
		voivod <- logLikelihoodFullModel(myData = myData, cutPoints = cutPoints[iter, 0:L], theta = overkill, startPoint = startPoint)
		propLogL <- voivod$logL
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
			x_minus_mean2 <- voivod$x_minus_mean2
		}
###########################################################################################################
#		Gibss move: update all theta[-augmentedCutPoints] from the prior.
		if(Prior == "complexity"){
		if(iter > burn){
			overkill <-  simMultiIndNormInvGamma(mu = priorParameters$mu0, 
								nu = priorParameters$nu0, 
								alpha = 2, 
								beta = sigma2[iter-1, ]
							)
			if(L < nTime - 1){
				theta$mean[-augCutPoints] <- overkill$mean[-augCutPoints]
			}
		}
		}else{
			overkill <-  simMultiIndNormInvGamma(mu = priorParameters$mu0, 
								nu = priorParameters$nu0, 
								alpha = 2, 
								beta = sigma2[iter-1, ]
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
				voivod <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
				propLogL <- voivod$logL
				mhRatio <- propLogL - lValues[iter - 1]
				u <- runif(1)
				if(log(u) < mhRatio){
					cutPoints[iter, 1:L] <- newCutPoints
					lValues[iter] <- propLogL
					mhRate <- mhRate + length(movesRange)
					x_minus_mean2 <- voivod$x_minus_mean2
				}
			}
		}
		if(chooseMove == "2"){
#			Move b: local random walk
			if(L > 0){
				lpProp <- localProposal(cutPoints = cutPoints[iter, 1:L] , nTime = nTime, mhPropRange = mhPropRange, startPoint = startPoint)
				newCutPoints <- lpProp$newState
				voivod <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
				propLogL <- voivod$logL
				mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
					logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
					logPrior(cutPoints = cutPoints[iter, 1:L], nTime = nTime, startPoint = startPoint)

				u <- runif(1)
				if(log(u) < mhRatio){
					cutPoints[iter, 1:L] <- newCutPoints
					lValues[iter] <- propLogL
					mhRate2 <- mhRate2 + length(movesRange)
					x_minus_mean2 <- voivod$x_minus_mean2
				}
			}
		}
		if(chooseMove == "3"){
#			Move c: random walk to just one index
			if( L > 0 ){
				lpProp <- singleLocalProposal(cutPoints = cutPoints[iter, 1:L], nTime = nTime, 
						mhSinglePropRange = mhSinglePropRange, startPoint = startPoint)
				newCutPoints <- lpProp$newState
				voivod <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
				propLogL <- voivod$logL
				mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
					logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
					logPrior(cutPoints = cutPoints[iter, 1:L], nTime = nTime, startPoint = startPoint)

				u <- runif(1)
				if(log(u) < mhRatio){
					cutPoints[iter, 1:L] <- newCutPoints
					lValues[iter] <- propLogL
					mhRate3 <- mhRate3 + length(movesRange)
					x_minus_mean2 <- voivod$x_minus_mean2
				}
			}
		}
		if(chooseMove == "4"){
#			Move c: random walk to just one index
			if(L > 0){
				lpProp <- singleLocalProposal(cutPoints = cutPoints[iter, 1:L], nTime = nTime, 
						mhSinglePropRange = 2, startPoint = startPoint)
				newCutPoints <- lpProp$newState
				voivod <- logLikelihoodFullModel(myData = myData, cutPoints = newCutPoints, theta = theta, startPoint = startPoint)
				propLogL <- voivod$logL
				mhRatio <- propLogL - lValues[iter] + lpProp$propRatio + 
					logPrior(cutPoints = newCutPoints, nTime = nTime, startPoint = startPoint) - 
					logPrior(cutPoints = cutPoints[iter, 1:L], nTime = nTime, startPoint = startPoint)

				u <- runif(1)
				if(log(u) < mhRatio){
					cutPoints[iter, 1:L] <- newCutPoints
					lValues[iter] <- propLogL
					mhRate3 <- mhRate3 + length(movesRange)
					x_minus_mean2 <- voivod$x_minus_mean2
				}
			}
		}
		########################################################################################
#		Update variance  via gibbs sampling:
		if(iter > 30000){
		sINV_rate = 0.5*(x_minus_mean2 + priorParameters$nu0*(theta$mean - priorParameters$mu0)^2 + 2*priorParameters$beta0)
		sINV_shape = (nReps + 1)/2 + priorParameters$alpha0
		sigma2[iter, ] = 1/rgamma(n = nTime, shape = sINV_shape, rate = sINV_rate )
		theta$variance <- sigma2[iter, ]
		}else{
			sigma2[iter, ] = theta$variance
		}
		########################################################################################

		lPosterior[iter] <- lValues[iter] + logPrior(cutPoints = cutPoints[iter, 0:L], nTime = nTime, startPoint = startPoint)
		if(saveTheta == TRUE){
			thetaValues[iter, ] <- theta$mean
		}
		if(iter %% myFl == 0){cat(paste0(uchar, " "))}
	}
	acceptanceRates <- c( mhRatePoints*100/iter,  mhRate0c*100/iter, mhRate2*100/iter, mhRate3*100/iter )
	names(acceptanceRates) <- c("move_l_n", "move_theta_randomWalk", "move_tau_all", "move_tau_j")
	cat(paste0(" Accepted MH moves: [l_n] ", round(mhRatePoints*100/iter,2), "%, [","\U03B8","+","\U03B5",'] ',round(mhRate0c*100/iter,2) ,"%, [","\U03C4","] ", round(mhRate2*100/iter,2), "%, [","\U03C4","_j] ", round(mhRate3*100/iter,2), "%."),"\n")
	results <- vector("list", length = 8)
	results[[2]] <- lValues
	results[[3]] <- acceptanceRates
	results[[4]] <- Lvalues
	LVtab <- table(Lvalues[-(1:burn)])/(length(Lvalues[-(1:burn)]))
	results[[6]] <- LVtab
	names(results) <- c("cutPoints", "logLikelihood", "acceptanceRates", "nCutPoints", "nCutPointsMAP", "nCutPointsPosterior", "theta", "variance")
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
		results[[8]] <- sigma2
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
		close.screen( all.screens = TRUE )
		dev.off()
	}

	return(results)
}



beast <- function(myDataList, burn = 20000, nIter = 70000, mhPropRange = 1, mhSinglePropRange=40, startPoint = 2, 
		timeScale, savePlots, zeroNormalization = FALSE, LRange = 0:30, tau = 0.05,
		gammaParameter = 2, nu0 = 0.1, alpha0 = 1, beta0 = 1, subsetIndex, saveTheta = TRUE, sameVariance = TRUE, Prior = "complexity"){
#	burn = 2000, nIter = 5000,mhPropRange = 1, mhSinglePropRange = 50, getSDvalues = T, startPoint=54, timeScale = 1/6,
	cat("\n")
	myColNames <- colnames(myDataList[[1]])
	if(missing(timeScale)){timeScale = 1/1}
#	if(missing(zeroNormalization)){zeroNormalization = TRUE}
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
	getSDvalues = TRUE
	if(Prior == "complexity"){
		cat(paste0("*  Assuming a complexity prior distribution on the number of change-points with `gammaParameter = ", gammaParameter,"`."), "\n")	
	}else{
		cat(paste0("*  Assuming a Poisson prior distribution on the number of change-points with rate `gammaParameter = ", gammaParameter,"`."), "\n")
		cat(paste0("*      [WARNING]: The Poisson distribution is NOT suggested!"), "\n")		
	}
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
	L = 1	# not really used anywhere
#	zz <- file(paste0(savePlots,"/numberCutPointsMAP.txt"), open = "w")
	nReps <- length(myDataList)
	n <- dim(myDataList[[1]])[2]
	cutPoints <- vector("list", length = n)
	cutPointsVar <- vector("list", length = n)
	if(nReps < 2){stop("no replicates")}
	cat(paste0("*  MCMC sampler parameters: nIter = ", nIter, ", burn = ", burn, ", startPoint = ", startPoint ),"\n")
	cat(paste0("*  Running MCMC for ", n, " subjects..."), "\n")	
	NAindex <- c()
	results <- vector("list", length = 11)
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
	}else{
		subsetIndex <- 1:n
	}
	movesRange = c(2,3)
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
							iterPerPlotPrefix = paste0("plot-", i), priorParameters = priorParameters, 
							L = L, LRange = LRange, tau = tau, 
							gammaParameter = gammaParameter, saveTheta = saveTheta, Prior = Prior)
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
#			cat(paste0(myColNames[i], " ", mhRunForSubject$nCutPointsMAP), file = zz, "\n")
		}
	}
#	close(zz)
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
				"nCutPointsTrace",
				"subsetIndex",
				"data")
	results$subsetIndex = subsetIndex
	results$data = myDataList
	cat(paste0("*  ALL DONE."),"\n")
	if(is.null(savePlots) == FALSE){ 
		cat(paste0("*  See produced *.pdf files in: `",getwd(),"/",savePlots,"`"),"\n")
	}

        class(results) <- c('list', 'beast.object')
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


computePosteriorParameters <- function(myDataList, priorParameters){
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
		j <- j + 1
		alive <- c(alive, i)
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
computePosteriorParametersFree <- function(myDataList, priorParameters){
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
		j <- j + 1
		alive <- c(alive, i)
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

computeEmpiricalPriorParameters <- function(myDataList, nu0 = 1, alpha0 = 1, beta0 = 1){
	priorParameters <- vector("list",length = 4)
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
		j <- j + 1
		alive <- c(alive, i)
		xMean <- xMean + colMeans(myData)
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

#' @export
print.beast.object <- function(x, ...){
        if( 'beast.object' %in% class(x) ){
                cat("\n")
		subsetIndex <- x$subsetIndex
		cat(paste0("Table of frequencies of the most probable number of change-points:"),'\n')
		print(table(names(unlist(lapply(x$NumberOfCutPoints_MAP, function(y){table(y)})))[subsetIndex]))
		cat("\n")
		cat(paste0("Summary per time-series:"),'\n')
		for(i in subsetIndex){
			cat(paste0(x$subject_ID[[i]],": "), 
				paste0("P(",x$NumberOfCutPoints_MAP[[i]], " change-points|data) ="), 
				paste0(round(as.numeric(x$NumberOfCutPoints_posterior_distribution[[i]][as.character(x$NumberOfCutPoints_MAP[[i]])]),2)), 
				paste0("located at"),paste0( x$Cutpoint_posterior_median[[i]]),
			 "\n")
		}
                cat('\n')
        }else{
                cat(paste("    The input is not in class `beast.object`"),'\n')
        }
}


#' @export
plot.beast.object <- function(x, fileName, width = 9, height = 6, pointsize = 12, ylab = "x", xlab = "t", timeScale = 1, myPal, ...){
        if( 'beast.object' %in% class(x) ){
                cat("\n")
		if(timeScale < 0){stop("timeScale should be positive.")}
		subsetIndex <- x$subsetIndex
		mutantNames <- colnames(x$data[[1]])[subsetIndex]	
		nTime <- dim(x$data[[1]])[1]
		n <- dim(x$data[[1]])[2]
		nReps <- length(x$data)
		xMax <- max(unlist(x$data))
		myPositions <- floor(seq(0, nTime, length = 10))  
		myLabels <- round(myPositions*timeScale,1)
		pdf(fileName, width = width, height = height, pointsize = pointsize)
#		1st: summary plot
		nColors <- length(table(names(unlist(lapply(x$NumberOfCutPoints_MAP, function(y){table(y)})))[subsetIndex]))
		if(nColors < 10){
			myPal <- brewer.pal(n = nColors, name = "Set1")
		}else{
			stop("the range of possible numbers of change-points is larger than 9, please manually define the colors to `myPal` argument")
		}
		myMeanTimeSeries <- array(data = 0, dim = c(n, nTime))
		for (i in subsetIndex){
			for (j in 1:nReps){
				myMeanTimeSeries[i,] <- myMeanTimeSeries[i,] + x$data[[j]][ ,i]
			}
			myMeanTimeSeries[i,] <- myMeanTimeSeries[i,]/nReps
		}
		if(timeScale == 1){
			matplot(t(myMeanTimeSeries[subsetIndex,]), type = "l", col = myPal[as.numeric(as.factor(unlist(x$NumberOfCutPoints_MAP)[subsetIndex]))], lty = 1, xlab = "t", ylab = paste0("average ", ylab))
		}else{
			matplot(t(myMeanTimeSeries[subsetIndex,]), xaxt = "n", type = "l", col = myPal[as.numeric(as.factor(unlist(x$NumberOfCutPoints_MAP)[subsetIndex]))], lty = 1, xlab = "t", ylab = paste0("average ", ylab))
			axis(1, at=myPositions,labels=myLabels)
		}
		legend("topleft", col = c(0,myPal), lty = 1, c("MAP number of change-points:", paste0(names(table(unlist(x$NumberOfCutPoints_MAP)[subsetIndex]))," (",as.numeric(table(unlist(x$NumberOfCutPoints_MAP)[subsetIndex]))," time-series)")))
		
#		then, individual plots
			iter <- 0
			for (i in subsetIndex){
				iter <- iter+1
				tmp <- x$data[[1]][,i] 
				if(nReps > 1){
					for(j in 2:nReps){
						tmp <- cbind(tmp, x$data[[j]][,i])
					}
				}
				if(timeScale == 1){
					matplot(tmp, 
							type = "l", lwd = 2, lty = 1, col = 2:(nReps+1), ylim = c(0,xMax), 
							ylab = ylab, xlab = xlab, main = paste0( "",mutantNames[iter],""))
				}else{
					matplot(tmp, xaxt = "n",
							type = "l", lwd = 2, lty = 1, col = 2:(nReps+1), ylim = c(0,xMax), 
							ylab = ylab, xlab = xlab, main = paste0( "",mutantNames[iter],""))
					axis(1, at=myPositions,labels=myLabels)
				}
				abline(v = x$Cutpoint_posterior_median[[i]], lty = 3)
				l <- x$NumberOfCutPoints_MAP[[i]]
				yRange <- c(0, max(tmp))
				yPoints <- seq(from = 0, to = xMax,length=l+2)[-c(1,l+2)]
				if(l == 1){
					boxplot(x$Cutpoint_mcmc_trace_map[[i]], horizontal = TRUE, add = TRUE, at = yPoints, boxwex = 0.3, col = "gray", xaxt = "n")
				}else{
					for (j in 1:l){ 
						boxplot(x$Cutpoint_mcmc_trace_map[[i]][,j], horizontal = TRUE, add = TRUE, at = yPoints[j], boxwex = 0.3, col = "gray", xaxt = "n")
					}
				}
			}
		dev.off()
		cat(paste0("See produced file: ", fileName),"\n")
        }else{
                cat(paste("    The input is not in class `beast.object`"),'\n')
        }
}




