library('beast')        # load package
library('not')


#**************************************************************************
# Reproduce Figure 1 of manuscript
#**************************************************************************

data("FungalGrowthDataset")     # load dataset
myIndex <- c(392, 62, 3, 117)   # run the sampler only for the 
#                                 specific subset of time-series

set.seed(1)             # optional
# Run MCMC sampler with the default number of iterations (nIter =70000):
run_mcmc <- beast(myDataList = FungalGrowthDataset, subsetIndex = myIndex, 
                     zeroNormalization = TRUE) 
# Print output:
print(run_mcmc)
# Plot output to file: "beast_plot.pdf"
plot(run_mcmc, fileName = "beast_plot_fungal_data.pdf", timeScale=1/6, xlab = "hours", ylab = "growth")



#**************************************************************************
# Reproduce Figure D.10 of Appendix 
# (comparison with the `not` package of Baranowski et al, 2016)
#**************************************************************************
# (a) We consider one example of the Fungal Growth dataset.
#	In order to apply `not` we have to average the data across the three replicates
i <- myIndex[2]
myData1 <- cbind(FungalGrowthDataset[[1]][,i],
		FungalGrowthDataset[[2]][,i],
		FungalGrowthDataset[[3]][,i])

#-----	apply `not`
notData1 <- rowMeans(myData1)
w <- not(notData1, contrast = "pcwsLinContMean") 
fo1 <- features(w)

#	comparison of `not` vs `beast` (as shown in 'real data' of Figure D.10 in the appendix):
par(mfrow=c(1,2))
matplot(myData1, col = 2:4, type = 'l', lty = 1, xlab = 't', ylab = 'x', main = 'proposed method (real data)')
abline(v = run_mcmc$Cutpoint_posterior_median[[62]],  lty = 3)

plot(notData1, type ='l', lty = 1, xlab = 't', ylab = 'average of x', main = '"not" package (real data)')
abline(v = fo1$cpt,  lty = 3, xlab = 't', ylab = 'x')


# (a) Simulated data example
#	Simulate N = 1 time-series with 3 replicates and 3 change-points

source('simulate.R')

#-----	apply proposed method
set.seed(1)
run_mcmc_sim <- beast(myDataList = myDataList, 
                             zeroNormalization = TRUE, nIter = 40000, burn = 20000) 

#-----	apply `not` (on averaged data across replicates)
i <- 1
myData2 = myDataList[[1]][,i] 
for(j in 2:nReps){
        myData2 <- cbind(myData2, myDataList[[j]][,i])
}
notData2 <- rowMeans(myData2)
w <- not(notData2, contrast = "pcwsLinContMean") 
fo2 <- features(w)

#	comparison of `not` vs `beast` (as shown in 'simulated data' of Figure D.10 in the appendix):
par(mfrow=c(1,2))
matplot(myData2, col = 2:4, type = 'l', lty = 1, xlab = 't', ylab = 'x', main = 'proposed method (synthetic data)')
abline(v = run_mcmc_sim$Cutpoint_posterior_median[[1]],  lty = 3)

plot(notData2, type ='l', lty = 1, xlab = 't', ylab = 'average of x', main = '"not" package (synthetic data)')
abline(v = fo2$cpt,  lty = 3, xlab = 't', ylab = 'x')





