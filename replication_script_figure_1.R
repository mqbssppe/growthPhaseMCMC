#-------------------------------------------------------------------------------------------------------------------------------------
# (1) Download source code:
# 	https://github.com/mqbssppe/growthPhaseMCMC/blob/master/fullMCMC_complexityPrior.R
# (2) Download and extract the data:
#	https://github.com/mqbssppe/growthPhaseMCMC/blob/master/fungal_time_series_dataset.tar.gz
# (3) Make sure that 'fullMCMC_complexityPrior.R' and the extracted folder 'fungal_time_series_dataset' are in the working directory
#-------------------------------------------------------------------------------------------------------------------------------------



source('fullMCMC_complexityPrior.R')

# Read data:
cat("Reading input data", "\n")
data_folder <- "fungal_time_series_dataset" # contains one *txt file per replicate
data_files <- list.files(data_folder)
nReps <- length(data_files)
myDataList <- vector("list", length = nReps)
i <- 0
for (myFile in data_files){
	i <- i + 1
	f <- paste0(data_folder,"/",myFile)
	cat(paste0("  Reading file ", f),"\n")
	myDataList[[i]] <- read.table(f)
}
cat("Done.", "\n")
cat(paste0("Number of different time-series: ", dim(myDataList[[1]])[2]), "\n")
cat(paste0("Length of each time-series: ", dim(myDataList[[1]])[1]), "\n")
cat(paste0("Number of replicates: ", nReps), "\n")

# Run the MCMC sampler for a subset of four time-series:
myIndex <- c(433, 64, 3, 124)
# Specify the main parameters of model and sampler:
gP <- 2		# $\alpha$ 
nu0 <- 0.1	# $\nu_0$
alpha0 <- 1	# $\alpha_0$
beta0 <- 1	# $\beta_0$
burn <- 20000	# burn-in period
nIter <- 70000	# total MCMC iterations
LRange = 0:30	# range of possible number of change-points
tmp_folder <- "temp_folder"	# temporary folder
run_mcmc <- growthPhaseFullMCMC(myDataList = myDataList, burn = burn, nIter = nIter, blankThreshold = 0.02, 
                        mhSinglePropRange = 40, savePlots = tmp_folder, zeroNormalization = TRUE,
                        showProgress = FALSE, movesRange = c(2, 3), L = 3, LRange = LRange, tau = 0.05, 
                        gammaParameter = gP, nu0 = nu0, alpha0 = alpha0, beta0 = beta0, subsetIndex = myIndex)

mutantNames <- colnames(myDataList[[1]])[myIndex]
mutantNames <- unlist(lapply(strsplit(mutantNames, split = "_"), function(y){paste(y, collapse = " ")}))
system(paste0("rm -r ", tmp_folder))
# Manuscript Figure 1:
cat(paste0("Producing Figure 1 of the manuscript: fig1.pdf"), "\n")
pdf(file = "fig1.pdf", width = 11, height = 8)
	par(mfrow = c(2,2), mar = c(4,4,3,0.6))
	nTime <- dim(myDataList[[1]])[1]
	n <- dim(myDataList[[1]])[2]
	myMeanTimeSeries <- array(data = NA, dim = c(n, nTime))
	for (i in 1:n){
		myMeanTimeSeries[i,] <- (myDataList[[1]][ ,i] + myDataList[[2]][ ,i] + myDataList[[3]][ ,i])/nReps
	}

	iter <- 0
	for (i in myIndex){
		iter<-iter+1
		matplot(cbind(myDataList[[1]][,i], myDataList[[2]][,i], myDataList[[3]][,i]), 
				type = "l", lwd = 2, lty = 1, col = 2:4, ylim = c(0,1.5), 
				ylab = "growth level", xlab = "t", main = paste0( "",mutantNames[iter],""))
		abline(v = run_mcmc$Cutpoint_posterior_median[[i]], lty = 3)
		l <- run_mcmc$NumberOfCutPoints_MAP[[i]]
		yPoints <- myMeanTimeSeries[i,run_mcmc$Cutpoint_posterior_median[[i]]] + 0.2
		if(l == 1){
		        boxplot(run_mcmc$Cutpoint_mcmc_trace_map[[i]], horizontal = TRUE, add = TRUE, at = yPoints, boxwex = 0.3, col = "gray", xaxt = "n")
		}else{
		        for (j in 1:l){ 
		                boxplot(run_mcmc$Cutpoint_mcmc_trace_map[[i]][,j], horizontal = TRUE, add = TRUE, at = yPoints[j], boxwex = 0.3, col = "gray", xaxt = "n")
		        }
		}
		myProb <- paste0("P(", run_mcmc$NumberOfCutPoints_MAP[[i]], " cutpoints | data)", " = ",
				round(as.numeric(run_mcmc$NumberOfCutPoints_posterior_distribution[[i]][ 
							as.character(run_mcmc$NumberOfCutPoints_MAP[[i]]
							)
						]), 
					2))
	}
dev.off()
cat("Done.", "\n")
