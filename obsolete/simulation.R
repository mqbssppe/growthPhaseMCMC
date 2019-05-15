set.seed(10)
nTime <- 500
nReps <- 3
n <- 100
sigma2 <- rgamma(n = nTime, shape = 1, rate = seq(1,0.1, length = nTime))
myDataList <- vector("list", length = nReps)

meanLinearFunction <- function(theta1, theta2, t1, t2){
	myMean <- theta1 + (t1:t2 - t1)*(theta2 - theta1)/(t2 - t1)
	return(myMean)
}


myMean <- vector("list", length = n)
Ncuts <- numeric(n)
for(i in 1:n){
	myMean[[i]] <- array(data = 0,dim = c(nReps,nTime))
	nCutPoints <- sample(0:10, 1)
	Ncuts[i] <- nCutPoints
	cutPoints <- numeric(nCutPoints)
	if(nCutPoints > 0){
		cutPoints <- sample(10:(nTime - 10), Ncuts[i], replace = FALSE)
	}
	cutPoints <- c(sort(cutPoints),nTime) #just in case
	cutValues <- numeric(nCutPoints + 1)
	if(nCutPoints > 0){
		cutValues[1] <- -0.01 + 0.02*runif(1)
		for(l in 2:(nCutPoints+1)){
#			easy scenario
			if(l%%2 == 0){
				cutValues[l] <- rnorm(1, mean = 15 + l, sd = 5)
			}else{
				cutValues[l] <- rnorm(1, mean = -15 - l, sd = 5)
			}
			if(cutPoints[l] - cutPoints[l-1] < 15){
				cutValues[l] <- - cutValues[l]*(cutValues[l-1]/abs(cutValues[l-1]))
			}


#			ligo pio diskolo
#			mySlope <- rnorm(1, mean = 0, sd = 0.3) 
#			cutValues[l] <- cutValues[l-1] + mySlope*(cutPoints[l] - cutPoints[l-1])
		}
	}
	nCutPoints <- nCutPoints + 1	
	for(k in 1:nReps){
		cutPointsReplicate <- cutPoints
		if(nCutPoints > 1){
			cutPointsReplicate[1:(nCutPoints-1)] <- cutPoints[1:(nCutPoints-1)] + sample(seq(-1,1),nCutPoints - 1, replace = TRUE)*rpois(n = nCutPoints - 1, lambda = 1)
			if( min(diff(cutPointsReplicate)) < 1){cutPointsReplicate <- cutPoints}
			if( min(cutPointsReplicate) < 3){cutPointsReplicate <- cutPoints}
			for(j in 2:nCutPoints){
#				myMean[[i]][ k, (cutPointsReplicate[j-1]+1):(cutPointsReplicate[j])] <- seq(cutValues[j-1] + rnorm(n = 1, mean = 0, sd = 0.3), cutValues[j] + rnorm(n = 1, mean = 0, sd = 0.3), length = cutPointsReplicate[j] - cutPointsReplicate[j-1]  )
			myMean[[i]][ k, (cutPointsReplicate[j-1]):(cutPointsReplicate[j])] <- meanLinearFunction(cutValues[j-1] + rnorm(n = 1, mean = 0, sd = 0.3), cutValues[j] + rnorm(n = 1, mean = 0, sd = 0.3), cutPointsReplicate[j-1], cutPointsReplicate[j])
			}
		}else{
			myMean[[i]][ k, ] <- rep(0,nTime)
		}

	}
}
write.table(Ncuts, file = "NcutPoints_true.txt")
for(k in 1:nReps){
	myDataList[[k]] <- array( data = 0, dim = c(nTime, n))
	for(i in 1:n){
		myDataList[[k]][ ,i] <- myMean[[i]][k, ] + rnorm(n = nTime)*sqrt(sigma2)
	}
	colnames(myDataList[[k]]) <- paste0("item_",1:n)
}

source('~/Dropbox/all_time_course_data/growthPhaseMCMC/growthPhaseMCMC/fullMCMC_complexityPrior.R')
gP <- 0.2
generalSampler <- growthPhaseFullMCMC(myDataList = myDataList, burn = 50000, nIter = 200000, 
			mhSinglePropRange = 5, savePlots = paste0("gamma_",gP), zeroNormalization = FALSE,
			showProgress = FALSE, movesRange = c(2, 3), L = 3, LRange = 0:30, tau = 0.1, 
			gammaParameter = gP, nu0 = 1, alpha0 = 1, beta0 = 1)
save.image(paste0("gamma_",gP,".RData"))

