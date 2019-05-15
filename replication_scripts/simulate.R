##########################################################################################################################################
# SIMULATE DATA
##########################################################################################################################################
# Simulation parameters
p <- 0.075       # probability of same sign with previous                                ####
nReps <- 3      # number of replicates                                                  ####
############################################################################################
#                                                                                       ####
############################################################################################

set.seed(3)
nTime <- 1000
n <- 1
myDataList <- vector("list", length = nReps)
cutPointsList <- vector("list", length = n)


switchSignMatrix <- array(data = NA, dim = c(2,2))
colnames(switchSignMatrix) <- rownames(switchSignMatrix) <- c("-","+")
switchSignMatrix[1, ] <- c(p, 1-p)
switchSignMatrix[2, ] <- c(1-p, p)
sigma2 <- array(data = NA, dim = c(n, nTime))

myMean <- vector("list", length = n)
Ncuts <- numeric(n)

#ss2 <- rgamma(n = nTime, shape = 1, rate = seq(1,0.1, length = nTime))

for(i in 1:n){
	set.seed(75+i)
	sigma2[i,] <- rgamma(n = nTime, shape = 1, rate = seq(1,0.1, length = nTime))
#	sigma2[i,] <- ss2
	myMean[[i]] <- array(data = 0,dim = c(nReps,nTime))
	nCutPoints <- sample(0:10, 1)
	Ncuts[i] <- nCutPoints
	cutPoints <- numeric(nCutPoints)
	if(nCutPoints > 0){
		centerPoints <- floor(seq(0,nTime,by = nTime/(Ncuts[i]+1)))[-c(1,Ncuts[i]+2)]	
		mySeq0 <- ((50-20):(50+20))
		myBinProbs <- dbinom(x= mySeq0, size = 2*50, prob = 0.5)
		for(l in 1:Ncuts[i]){
			cutProbs_given_Ncuts <- rep(0,nTime)
			mySeq <- ((centerPoints[l]-20):(centerPoints[l]+20))
			cutProbs_given_Ncuts[mySeq] <- myBinProbs
			cutPoints[l] <- sample(nTime, 1, prob = cutProbs_given_Ncuts)
		}
		cutPointsList[[i]] <- cutPoints
	}
	cutPoints <- c(sort(cutPoints),nTime) #just in case
        cutValues <- numeric(nCutPoints + 1)
        if(nCutPoints > 0){
                cutValues[1] <- -0.01 + 0.02*runif(1)
		precSign <- "-"
		if(runif(1) < 0.5){precSign <- "+"}
                for(l in 2:(nCutPoints+1)){
			u <- runif(1)
			if (u < switchSignMatrix[precSign, "-"]){
				currentSign <- "-"
				mySlope <- - abs(rnorm(1, mean = 0, sd = 0.3))
			}else{
				currentSign <- "+"
				mySlope <- abs(rnorm(1, mean = 0, sd = 0.3))
			}
			cutValues[l] <- cutValues[l-1] + mySlope*(cutPoints[l] - cutPoints[l-1])
			precSign <- currentSign
                }
        }
	nCutPoints <- nCutPoints + 1	
	for(k in 1:nReps){
		cutPointsReplicate <- cutPoints
		if(nCutPoints > 1){
			cutPointsReplicate[1:(nCutPoints-1)] <- cutPoints[1:(nCutPoints-1)] + sample(seq(-1,1),nCutPoints - 1, replace = TRUE)*rpois(n = nCutPoints - 1, lambda = 2)
			if( min(diff(cutPointsReplicate)) < 1){cutPointsReplicate <- cutPoints}
			if( min(cutPointsReplicate) < 3){cutPointsReplicate <- cutPoints}
			for(j in 2:nCutPoints){
				myMean[[i]][ k, (cutPointsReplicate[j-1]+1):(cutPointsReplicate[j])] <- seq(cutValues[j-1] + rnorm(n = 1, mean = 0, sd = 1), cutValues[j] + rnorm(n = 1, mean = 0, sd = 1), length = cutPointsReplicate[j] - cutPointsReplicate[j-1]  )
			}
		}else{
			myMean[[i]][ k, ] <- rep(0,nTime)
		}

	}
}
#write.table(Ncuts, file = "NcutPoints_true.txt")
for(k in 1:nReps){
	myDataList[[k]] <- array( data = 0, dim = c(nTime, n))
	for(i in 1:n){
		myDataList[[k]][ ,i] <- myMean[[i]][k, ] + rnorm(n = nTime)*sqrt(sigma2[i,])
	}
	colnames(myDataList[[k]]) <- paste0("item_",1:n)
}

# Produce supplementary figure 1
##########################################################################################################################################
# #	end of data simulation
##########################################################################################################################################

