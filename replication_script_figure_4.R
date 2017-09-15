#-------------------------------------------------------------------------------------------------------------------------------------
# In a terminal use 
# (1) > git clone https://github.com/mqbssppe/growthPhaseMCMC  
# (2) > cd growthPhaseMCMC
# (3) > tar -zxvf fungal_time_series_dataset.tar.gz
# (4) > R CMD BATCH replication_script_figure_4.R
#............................................................................................................
# or inside R> source('replication_script_figure_4.R')
#-------------------------------------------------------------------------------------------------------------------------------------



source('fullMCMC_complexityPrior.R')

# Read data:
cat("Reading input data", "\n")
data_folder <- "fungal_time_series_dataset" # contains one *txt file per replicate
data_files <- list.files(data_folder, pattern = "replicate")
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


# Specify the main parameters of model and sampler:
gP <- 2		# $\alpha$ 
nu0 <- 0.1	# $\nu_0$
alpha0 <- 1	# $\alpha_0$
beta0 <- 1	# $\beta_0$
burn <- 20000	# burn-in period
nIter <- 70000	# total MCMC iterations
LRange = 0:30	# range of possible number of change-points
tmp_folder <- "individual_plots"	#  folder
startPoint = 40 	# apply method imposing the prior restriction that the change-point of the first phase is at least equal to 40
#	this will take ~ 10 hours to complete
generalSampler <- growthPhaseFullMCMC(myDataList = myDataList, burn = burn, nIter = nIter, startPoint = startPoint,
                        mhSinglePropRange = 40, savePlots = tmp_folder, zeroNormalization = TRUE,
                        showProgress = FALSE, movesRange = c(2, 3), L = 3, LRange = LRange, tau = 0.05, 
                        gammaParameter = gP, nu0 = nu0, alpha0 = alpha0, beta0 = beta0)


# Manuscript Figure 4:
cat(paste0("Producing Figure 4 of the manuscript: fig4.pdf"), "\n")
nTime <- dim(myDataList[[1]])[1]
n <- dim(myDataList[[1]])[2]
myMeanTimeSeries <- array(data = 0, dim = c(n, nTime))
for (i in 1:n){
	for (j in 1:nReps){
	        myMeanTimeSeries[i,] <- myMeanTimeSeries[i,] + myDataList[[j]][ ,i]
	}
	myMeanTimeSeries[i,] <- myMeanTimeSeries[i,]/nReps
}
library(RColorBrewer)
myPal <- brewer.pal(n = 4, name = "Set1")
pdf(file = "fig4.pdf", width = 9, height = 6)
        matplot(t(myMeanTimeSeries), type = "l", col = myPal[generalSampler$NumberOfCutPoints_MAP], lty = 1, xlab = "t", ylab = "average growth level")
        legend("topleft", col = c(0,myPal), lty = 1, c("MAP number of change-points:", paste0(1:4," (",as.numeric(table(generalSampler$NumberOfCutPoints_MAP))," mutants)")))
dev.off()


myDF <- data.frame(read.table("fungal_time_series_dataset/AFUB_GENE_ID.txt", header=TRUE), NUMBER_OF_CHANGEPOINTS = generalSampler$NumberOfCutPoints_MAP)

cat("Done.", "\n")
