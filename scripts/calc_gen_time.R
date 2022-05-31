#### calculate average parental age from SLiM output

# M.I. Clark, 5/25/2022

# define file paths
filePath <- "../gen_time_output_05252022"
searchTerm <- "metaAll"

files <- list.files(filePath, full.names = TRUE) 
files <- files[grepl(paste0("*", searchTerm, "*"), files)] 
data <- vector(mode = "list", length = length(filePath))

# load in metadata
for (i in 1:length(files)){
	data[[i]] <- read.csv(files[i], sep = "\t", header = TRUE)
}

# calculate average per rep
avg <- as.data.frame(matrix(data=NA, ncol= 1, nrow=length(files)))
colnames(avg) <- c("average_parent_age")
for (i in 1:length(files)){
	avg[i,1] <- mean(c(data[[i]]$parent1Age, data[[i]]$parent2Age))
}

# calculate global average
write.table(mean(avg[,1]), paste0(filePath, "/gen_time.txt"), sep = "\t", row.names=FALSE, col.names=FALSE)

# plot distribution
pdf(paste0(filePath, "/gen_time_nWF_hist.pdf"))
hist(avg[,1])
dev.off()

### end script ### 
