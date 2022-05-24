library(MetBrewer)

# define custom functions -------------------------
load_het_data <- function(filePath, searchTerm, fileLength, reps){
  files <- list.files(filePath, full.names = TRUE)
  desired.files <- files[grepl(paste0("*", searchTerm, "*"), files)]
  data.frame <- as.data.frame(matrix(data=NA, ncol= length(desired.files), nrow=fileLength))
  colnames(data.frame) <- 1:length(desired.files)
  for(i in 1:length(desired.files)){
    data.frame[,i] <- read.csv(desired.files[i], sep = "\t")$het
    rownames(data.frame) <- read.csv(desired.files[i], sep = "\t")$generation
  }
  if(ncol(data.frame)==reps){
    return(data.frame)
  }else{print("Number of columns does not match expected reps... make sure your search term is specific enough!")}
  
}
load_Nc_data <- function(filePath, searchTerm, fileLength, reps){
  files <- list.files(filePath, full.names = TRUE)
  desired.files <- files[grepl(paste0("*", searchTerm, "*"), files)]
  data.frame <- as.data.frame(matrix(data=NA, ncol= length(desired.files), nrow=fileLength))
  colnames(data.frame) <- 1:length(desired.files)
  for(i in 1:length(desired.files)){
    data.frame[,i] <- read.csv(desired.files[i], sep = "\t")$Nc
  }
  if(ncol(data.frame)==reps){
    return(data.frame)
  }else{print("Number of columns does not match expected reps... make sure your search term is specific enough!")}
  
}
calc_avg_stdev <- function(data.frame){
  avg <- rowMeans(data.frame)
  dev <- apply(data.frame, 1, sd)
  upper <- avg + dev
  lower <-  avg - dev
  return(list(avg, dev, upper, lower))
}

# establish baseline het---------------------------
filePath <- "../het_test_output"
n.bot <- load_het_data(filePath, "no_bot_200k_[0-9].", 2000, 100)
stats.n.bot <- calc_avg_stdev(n.bot)

pdf("../viz_output/baseline_het.pdf", width = 8, height = 8)
plot(x = NULL, y = NULL, 
     xlab = "time (in generation x 100)", 
     ylab = "heterozygosity", 
     xlim = c(1, 2000), 
     ylim = c(0, 2e-4))
lines(x=1:2000, y=stats.n.bot[[1]], col = met.brewer("Archambault", 1))
polygon(x = c(1:2000, rev(1:2000)), y = c(stats.n.bot[[4]], rev(stats.n.bot[[3]])), col = adjustcolor(met.brewer("Archambault", 1), alpha.f = 0.2), border = NA)
abline(h=mean(stats.n.bot[[1]][1900:2000]), col = "blue", lty = 2)
abline(v=1800, col = "green", lty=2)
legend("topleft", c("nWF without bottleneck", "mean het over last 1000 generations", "generation 180,000"), lty = c(1, 2, 2), col = c(met.brewer("Archambault", 1), "blue", "green"))

#***probably should run the no bottleneck simulation longer
dev.off()
# find desired Nc value for pWF simulation
mu = 1e-8
het.nWF <- mean(stats.n.bot[[1]][1900:2000])
desired.Nc.pWF <- het.nWF / (4*mu)

het.pWF <- mean(stats.n.bot.pWF[[1]][1900:2000])
Ne.pwF <- het.pWF / (4*mu)

# ran with N = 4232 in pWF model
# pWF reps 55 and 56 did not finish running
n <- 4
palette <-"Egypt" #"Derain" #"Archambault" 
n.bot.pWF <- load_het_data(filePath, "no_bot_200K_pWF", 2000, 100)
stats.n.bot.pWF <- calc_avg_stdev(n.bot.pWF)
n.bot.pWF.batch <- load_het_data(filePath, "no_bot_200k_batch", 2000, 100)
stats.n.bot.pWF.batch <- calc_avg_stdev(n.bot.pWF.batch)
n.bot.pWF.pois <- load_het_data(filePath, "no_bot_200k_pois", 2000, 100)
stats.n.bot.pWF.pois <- calc_avg_stdev(n.bot.pWF.pois)

pdf("../viz_output/baseline_het_comparison.pdf", width = 8, height = 8)
plot(x = NULL, y = NULL, 
     xlab = "time (in generation x 100)", 
     ylab = "heterozygosity", 
     xlim = c(1, 2000), 
     ylim = c(0, 3e-4))
lines(x=1:2000, y=stats.n.bot[[1]], col = met.brewer(palette, n)[1])
polygon(x = c(1:2000, rev(1:2000)), y = c(stats.n.bot[[4]], rev(stats.n.bot[[3]])), col = adjustcolor(met.brewer(palette, n)[1], alpha.f = 0.2), border = NA)
lines(x=1:2000, y=stats.n.bot.pWF[[1]], col = met.brewer(palette, n)[2])
polygon(x = c(1:2000, rev(1:2000)), y = c(stats.n.bot.pWF[[4]], rev(stats.n.bot.pWF[[3]])), col = adjustcolor(met.brewer(palette, n)[2], alpha.f = 0.2), border = NA)
lines(x=1:2000, y=stats.n.bot.pWF.batch[[1]], col = met.brewer(palette, n)[3])
polygon(x = c(1:2000, rev(1:2000)), y = c(stats.n.bot.pWF.batch[[4]], rev(stats.n.bot.pWF.batch[[3]])), col = adjustcolor(met.brewer(palette, n)[3], alpha.f = 0.2), border = NA)
lines(x=1:2000, y=stats.n.bot.pWF.pois[[1]], col = met.brewer(palette, n)[4])
polygon(x = c(1:2000, rev(1:2000)), y = c(stats.n.bot.pWF.pois[[4]], rev(stats.n.bot.pWF.pois[[3]])), col = adjustcolor(met.brewer(palette, n)[4], alpha.f = 0.2), border = NA)

#abline(h=mean(stats.n.bot[[1]][1900:2000]), col = "blue", lty = 2)
#abline(v=1800, col = "green", lty=2)

legend("topleft", c("nWF without bottleneck, Nc ~ 7500", "pWF binomial reproduction, Nc= 4232", "pWF, batch reproduction, Nc= 4232", "pWF, poisson reproduction, Nc= 4232"), lty = c(1, 1, 1, 1), col = c(met.brewer(palette,  n)))
dev.off()


# pWF and nWF comparisons---------------------
n <- 3 # number of runs to compare
palette <-"Egypt" #"Derain" #"Archambault" 

filePath <- "../het_test_output"
nWF <- load_het_data(filePath, "nWF_decline", 11899, 100)
pWF <- load_het_data(filePath, "pWF_decline", 11899, 98)

stats.nWF <- calc_avg_stdev(nWF) # avg, dev, upper, lower
stats.pWF <- calc_avg_stdev(pWF)

decline <- cbind.data.frame(1:23, nWF.Nc.stats[[1]][34000:34022], round((nWF.Nc.stats[[1]][34000:34022]/7500)*100, digits = 1))
colnames(decline) <- c("count", "Nc", "percent of N")

pdf("../viz_output/het_testing_run_comparisons_decline_zoom_fit.pdf", width = 8, height = 8)######
plot(x = NULL, y = NULL, 
     xlab = "generation", 
     ylab = "heterozygosity", 
     xlim = c(178000, 185000), 
     ylim = c(0, 2e-4))

lines(x=as.numeric(rownames(nWF)), y=stats.nWF[[1]], col = met.brewer(palette, n)[1])
polygon(x = c(as.numeric(rownames(nWF)), rev(as.numeric(rownames(nWF)))), y = c(stats.nWF[[4]], rev(stats.nWF[[3]])), col = adjustcolor(met.brewer(palette, n)[1], alpha.f = 0.2), border = NA)
lines(x=as.numeric(rownames(pWF)), y=stats.pWF[[1]], col =met.brewer(palette, n)[2])
polygon(x = c(as.numeric(rownames(pWF)), rev(as.numeric(rownames(pWF)))), y = c(stats.pWF[[4]], rev(stats.pWF[[3]])), col = adjustcolor(met.brewer(palette, n)[2], alpha.f = 0.2), border = NA)

abline(v=180000, col = met.brewer(palette, n)[3], lty=2)

legend("topleft", c("nWF, N = 7500", "pWF, N = 4232, Poisson reproduction", "beginning of bottleneck"), lty = c(1, 1, 2) , col = met.brewer(palette, n))
dev.off()

# Nc --- --------------
#what does the decline look like
nWF.Nc <- load_Nc_data(filePath, "nWF_decline", 11899, 100)
pWF.Nc <- load_Nc_data(filePath, "pWF_decline", 11899, 98)
stats.nWF.Nc <- calc_avg_stdev(nWF.Nc) # avg, dev, upper, lower
stats.pWF.Nc <- calc_avg_stdev(pWF.Nc)

pdf("../viz_output/het_testing_run_comparisons_decline_zoom_fit.pdf", width = 8, height = 8)######
plot(x = NULL, y = NULL, 
     xlab = "generation", 
     ylab = "Nc", 
     xlim = c(179980, 180040), 
     ylim = c(0, 8000))
lines(x=(as.numeric(rownames(nWF))+1), y=stats.nWF.Nc[[1]], col = met.brewer(palette, n)[1])
lines(x=as.numeric(rownames(pWF)), y=stats.pWF.Nc[[1]], col = met.brewer(palette, n)[2])

abline(v=180000, col = met.brewer(palette, n)[3], lty=2)

lines(x=as.numeric(rownames(nWF)[1899:1921]), y = 7500*exp(1)^(-0.223006*0:22), lwd=10, col = adjustcolor(met.brewer(palette, n)[3], alpha.f = 0.2))
lines(x=(as.numeric(rownames(pWF)[1899:1921])), y = 4232*exp(1)^(-0.227636*0:22), lwd = 10, col = adjustcolor(met.brewer(palette, n)[3], alpha.f = 0.2))

legend("topright", c("nWF, N = 7500", "pWF, N = 4232, Poisson reproduction", "beginning of bottleneck", "exponential decline"), 
       lty = c(1, 1, 2, 2) , lwd = c(1, 1, 1, 10), col = c(met.brewer(palette, n), met.brewer(palette, n)[3]))

text(x= 180030, y= 2000, labels = "nWF: y=7500e^-0.223006x")
text(x= 180030, y= 1000, labels = "pWF: y=4232e^-0.223006x")

dev.off()

# how long is the decline in nWF? 

filePath <- "../het_test_output"
n=1
nWF.Nc <- load_Nc_data(filePath, "het_test_nWF", 35000, 100)
nWF.Nc.stats <- calc_avg_stdev(nWF.Nc)

plot(x = NULL, y = NULL, 
     xlab = "generation", 
     ylab = "Nc", 
     xlim = c(33990, 34025), 
     ylim = c(0, 8000))
lines(x=1:35000, y=nWF.Nc.stats[[1]], col = met.brewer("Archambault", n)[1])
polygon(x = c(1:35000, rev(1:35000)), y = c(nWF.Nc.stats[[4]], rev(nWF.Nc.stats[[3]])), col = adjustcolor(met.brewer("Archambault", n)[1], alpha.f = 0.2), border = NA)
abline(h=mean(nWF.Nc.stats[[1]][34025:35000]), col = "red", lty = 2)
abline(v=34022, col = "red", lty = 2)

nWF.Nc.stats[[1]][34015:34020]

# 34022 is when Nc stablizes to B --> 22 generations
decline <- cbind.data.frame(1:23, nWF.Nc.stats[[1]][34000:34022], round((nWF.Nc.stats[[1]][34000:34022]/7500)*100, digits = 1))
colnames(decline) <- c("count", "Nc", "percent of N")

# y = ab ^ x
plot(x = NULL, y = NULL, 
     xlab = "generation", 
     ylab = "Nc", 
     xlim = c(33990, 34025), 
     ylim = c(0, 4000))
points(x=1:35000, y=nWF.Nc.stats[[1]], col = met.brewer("Archambault", n)[1], pch = 19)

lines(x=34000:34022, y = 7500*exp(1)^(-0.223006*decline$count), col = "blue")
lines(x=34000:34022, y = 4232*exp(1)^(-0.227636*decline$count), col = "blue")

abline(h = 100, col = "red")
abline(h = 56, col = "blue")

cbind.data.frame(decline, 7495.15*exp(1)^(-0.227006*decline$count))



