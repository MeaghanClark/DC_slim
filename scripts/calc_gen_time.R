#### calculate average parental age from SLiM output

# M.I. Clark, written 5/25/2022
# last updated 8/18/2022

#####----------------------------------------------------------------------------------------------------------
# Define custom functions
#####----------------------------------------------------------------------------------------------------------
load_genData <- function(filePath, searchTerm, age){
  # testing vars
  # filePath <- "../genTime_output_08152022"
  # searchTerm <- "genTime"
  # age = 2
  
  files <- list.files(filePath, full.names = TRUE) 
  files <- files[grepl(paste0(searchTerm, "_", age, "_", "[:punct:]*"), files)] 
  data <- vector(mode = "list", length = length(filePath))
  
  # load in metadata
  for (i in 1:length(files)){
    data[[i]] <- read.csv(files[i], sep = "\t", header = TRUE)
  }
  return(data)
}

calc_genTime <- function(data){
  # calculate average per rep
  avg <- as.data.frame(matrix(data=NA, ncol= 1, nrow=length(files)))
  colnames(avg) <- c("average_parent_age")
  for (i in 1:length(files)){
    avg[i,1] <- mean(c(data[[i]]$parent1Age, data[[i]]$parent2Age))
  }
  return(avg)
}

#####----------------------------------------------------------------------------------------------------------

# define file paths
filePath <- "../genTime_output_08152022"
# searchTerm <- "genTime"
ages <- c(2, 5, 10, 20)

data_2 <- load_genData("../genTime_output_08152022", "genTime", ages[1])
data_5 <- load_genData("../genTime_output_08152022", "genTime", ages[2])
data_10 <- load_genData("../genTime_output_08152022", "genTime", ages[3])
data_20 <- load_genData("../genTime_output_08152022", "genTime", ages[4])

genTime_2 <- calc_genTime(data_2)
genTime_5 <- calc_genTime(data_5)
genTime_10 <- calc_genTime(data_10)
genTime_20 <- calc_genTime(data_20)

mean(genTime_2[,1])
mean(genTime_5[,1])
mean(genTime_10[,1])
mean(genTime_20[,1])

# calculate global average
write.table(mean(genTime_2[,1]), paste0(filePath, "/gen_time", ages[1], ".txt"), sep = "\t", row.names=FALSE, col.names=FALSE)
write.table(mean(genTime_5[,1]), paste0(filePath, "/gen_time", ages[2], ".txt"), sep = "\t", row.names=FALSE, col.names=FALSE)
write.table(mean(genTime_10[,1]), paste0(filePath, "/gen_time", ages[3], ".txt"), sep = "\t", row.names=FALSE, col.names=FALSE)
write.table(mean(genTime_20[,1]), paste0(filePath, "/gen_time", ages[4], ".txt"), sep = "\t", row.names=FALSE, col.names=FALSE)

# plot distribution
pdf(paste0(filePath, "/gen_time_nWF_hist.pdf"))
hist(genTime_20[,1])
dev.off()

### end script ### 
