

allPlates <- vector("list", length = 3)
for (plateNumber in 1:5){
l1 <- list.files(pattern = paste0("splitfile_plate",plateNumber))
l2 <- list.files(pattern = "Rep1")
l <- intersect(l1,l2)

df1 <- 0
for(myFile  in l){
	df1 <- cbind(df1, read.table(myFile, skip = 2, row.names = 1, header=T, sep = ";"))
}
df1 <- df1[,-1]

l1 <- list.files(pattern = paste0("splitfile_plate",plateNumber))
l2 <- list.files(pattern = "Rep2")
l <- intersect(l1,l2)

df2 <- 0
for(myFile  in l){
	df2 <- cbind(df2, read.table(myFile, skip = 2, row.names = 1, header=T, sep = ";"))
}
df2 <- df2[,-1]


l1 <- list.files(pattern = paste0("splitfile_plate",plateNumber))
l2 <- list.files(pattern = "Rep3")
l <- intersect(l1,l2)

df3 <- 0
for(myFile  in l){
	df3 <- cbind(df3, read.table(myFile, skip = 2, row.names = 1, header=T, sep = ";"))
}
df3 <- df3[,-1]
	if(plateNumber == 1){
	myDataList1 <- vector("list", length = 3)
	myDataList1[[1]] <- df1
	myDataList1[[2]] <- df2
	myDataList1[[3]] <- df3
	}
	if(plateNumber == 2){
	myDataList2 <- vector("list", length = 3)
	myDataList2[[1]] <- df1
	myDataList2[[2]] <- df2
	myDataList2[[3]] <- df3
	}
	if(plateNumber == 3){
	myDataList3 <- vector("list", length = 3)
	myDataList3[[1]] <- df1
	myDataList3[[2]] <- df2
	myDataList3[[3]] <- df3
	}
	if(plateNumber == 4){
	myDataList4 <- vector("list", length = 3)
	myDataList4[[1]] <- df1
	myDataList4[[2]] <- df2
	myDataList4[[3]] <- df3
	}
	if(plateNumber == 5){
	myDataList5 <- vector("list", length = 3)
	myDataList5[[1]] <- df1
	myDataList5[[2]] <- df2
	myDataList5[[3]] <- df3
	}
	colnames(df1) <- colnames(df2) <-colnames(df3) <- paste0("plate_",plateNumber,"_",colnames(df1))
	if(plateNumber == 1){
		allPlates[[1]] <- df1
		allPlates[[2]] <- df2
		allPlates[[3]] <- df3	
	}else{
		allPlates[[1]] <- cbind(allPlates[[1]], df1)
		allPlates[[2]] <- cbind(allPlates[[2]],df2)
		allPlates[[3]] <- cbind(allPlates[[3]],df3)
	}

}




source('~/Dropbox/all_time_course_data/growthPhaseMCMC/growthPhaseMCMC/growthPhaseMCMC.R')
growthCurveAnalysis <- growthPhaseMCMC(myDataList = allPlates, savePlots = "growthPlots", startPoint = 48)

source('~/Dropbox/all_time_course_data/growthPhaseMCMC/growthPhaseMCMC/growthPhaseMCMC_cutoff.R')
growthCurveAnalysis <- growthPhaseMCMC(myDataList = allPlates, savePlots = "growthPlots_cutoff", startPoint = 48, yCutoff = 0.3)


library(mvoutlier)
mahalanobisDistance <- numeric( dim(growthCurveAnalysis[[1]])[1] )
names(mahalanobisDistance) <- rownames(growthCurveAnalysis[[1]])
NAindex <- which(is.na(growthCurveAnalysis[[1]][,1])==TRUE)
myDF <- growthCurveAnalysis[[1]][-NAindex,]
#myDF <- cbind(myDF, growthCurveAnalysis[[3]][-NAindex,])
#myDF <- scale(myDF, scale = TRUE, center = TRUE)
dd <- dd.plot(myDF, quan=0.5)
mvOut <- aq.plot(myDF, alpha=0.01)
mahalanobisDistance[names(dd$md.rob)] <- dd$md.rob
mahalanobisDistance[names(dd$md.cla)] <- dd$md.cla
sortedList <- sort(mahalanobisDistance, decreasing=TRUE)

kreator <- array(data = NA, dim = dim(allPlates[[1]])[c(2,1)])
for(i in 1:length(kreator[,1])){
        kreator[i, ] <- rowMeans(cbind(allPlates[[1]][,i ], allPlates[[2]][,i ], allPlates[[3]][,i ]))
}
kreator <- kreator[-NAindex,]
myCol <- c("green","red")[1 + as.numeric(mvOut$outliers)]
matplot(t(kreator), col = myCol, type = "l", lty = 1, ylab = "average of 3 reps", xaxt = "n", main = paste0("all plates"), xlab = "hours")
axis(1, at = myPositions, labels = myLabels)






myPositions <- c(0,50,100,150,200,250,300)
myLabels <- round(c(0,50,100,150,200,250,300)/6,1)
maxColorValue <- max(mahalanobisDistance)
palette <- colorRampPalette(c("green","darkorange","red"))(maxColorValue)
myPerm <- order(mahalanobisDistance)
myCol <- palette[cut(mahalanobisDistance, maxColorValue)]
kreator <- array(data = NA, dim = dim(allPlates[[1]])[c(2,1)])
for(i in 1:length(kreator[,1])){
        kreator[i, ] <- rowMeans(cbind(allPlates[[1]][,i ], allPlates[[2]][,i ], allPlates[[3]][,i ]))
}
kreator <- kreator[myPerm, ]
blanks <- which(apply(kreator,1,max) < 0.03)
kreator <- kreator[-blanks,]
myCol <- myCol[myPerm]
myCol <- myCol[-blanks]
matplot(t(kreator), col = myCol, type = "l", lty = 1, ylab = "average of 3 reps", xaxt = "n", main = paste0("all plates"), xlab = "hours")
axis(1, at = myPositions, labels = myLabels)
        




myData <- cbind(nrmD[[1]][,137], nrmD[[2]][,137], nrmD[[3]][,137])
sdV <- getVariance(nrmD)
mr <- mhSampler(myData = myData, nIter = 20000, finalIterationPdf = F, mhPropRange = 1, mhSinglePropRange = 50, startPoint = 48, sdValues = sdV, timeScale = 1/6, showProgress = TRUE)


