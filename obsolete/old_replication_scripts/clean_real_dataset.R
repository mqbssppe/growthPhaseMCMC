cat("Reading input data", "\n")

# you should replace it with fungal_time_series_dataset_with_blanks
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

source('fullMCMC_complexityPrior.R')
myDataList <- normalizeTime0(myDataList = myDataList)

blankThreshold <- 0.023
n <- dim(myDataList[[1]])[2]
nTime <- dim(myDataList[[1]])[1]
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
	}
}

exclude36 <- which(alive == 36)
exclude75 <- which(alive == 75)
alive <- alive[-c(exclude36, exclude75)]

myDataList <- lapply(myDataList, function(y)y[,alive])
tmp <- read.table("~/Dropbox/Census/Growth_curves_full_media/Panos_Analysis/number_of_changePoints.txt", fill = TRUE)
AFUB_names <- as.character(tmp[alive,3])
AFUB_names <- cbind(alive, colnames(myDataList[[1]]),AFUB_names)
colnames(AFUB_names) <- c("ROW_NUMBER_ID", "PLATE_WELL_ID", "AFUB_GENE_ID")
write.table(file = paste0(data_folder,"/AFUB_GENE_ID.txt"),  AFUB_names)
write.table(file = paste0(data_folder,"/fungal_time_series_replicate_1.txt"), myDataList[[1]])
write.table(file = paste0(data_folder,"/fungal_time_series_replicate_2.txt"), myDataList[[2]])
write.table(file = paste0(data_folder,"/fungal_time_series_replicate_3.txt"), myDataList[[3]])

system("tar -zcvf fungal_time_series_dataset.tar.gz fungal_time_series_dataset")



