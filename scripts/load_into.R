#!/usr/bin/env Rscript

# Command line set up and define arguements------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

Sys.getenv("R_LIBS_USER")

date <- format(Sys.Date(), "%m%d%Y")

# capture command line variables
args <- commandArgs(trailingOnly = TRUE)
print(c("My args are ", args))

# args[1] is filepath, "../theta_output_dist_11082022"
# args[2] is outdir for Robj

indir <- args[2]
outdir <- args[3]
reps <- 100
ages <- c(2, 5, 10, 20)
rVals <- c(2, 10, 100)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load custom functions ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

load_ind_data <- function(filePath, model, reps, rVals, age_vec, keyword = "*het.txt", separator = " ", index = TRUE){
  # testing
  # filePath <- "../theta_output_dist_11082022"
  # model <- "nWF"
  # reps <- 100
  # rVals <- c(2, 10, 100)
  # age_vec <- 5
  # keyword = "*bins.txt"
  # separator = ","
  # index = FALSE
  #-----------------------------------------------
  if(missing(age_vec)){
    age_vec <- 1
  }
  files <- list.files(filePath, full.names = TRUE)
  ind_files <- files[grepl(paste0(keyword), files)]
  ind_files <- ind_files[grepl(paste0("*", model, "*"), ind_files)]
  
  #  if(model == "pWF"){
  #    ind_data <- vector(mode = "list", length = length(rVals)) # list with one element per rVal (3)
  #    for (r in rVals){
  #      desired.files <- ind_files[grepl(paste0("*", model, "_", r, "_[0-9]+_*"), ind_files)]
  #      data <- vector(mode = "list", length = length(reps)) # list with one element per rep (100)
  #      for (i in 1:length(desired.files)){
  #        df <- read.csv(desired.files[i], sep = " ", header = TRUE)
  #        data[[i]] <- df[-1]
  #      }
  #      ind_data[[which(rVals == r)]] <- data
  #    }
  #    names(ind_data) <- sapply(rVals, FUN = function(x){paste0("rVal_", x)})
  #    return(ind_data)
  #  }
  ind_data <- vector(mode = "list", length = length(age_vec)) # list with one element per age (4)
  for (a in age_vec){
    for (r in rVals){
      desired.files <- ind_files[grepl(paste0("*", model, "_", a, "_", r, "_[0-9]+_*"), ind_files)]
      data <- vector(mode = "list", length = length(reps)) # list with one element per rep (100)
      for (i in 1:length(desired.files)){
        df <- read.csv(desired.files[i], sep = separator, header = TRUE)
        if(index == TRUE){data[[i]] <- df[-1]}else{data[[i]] <- df}
      }
      ind_data[[which(age_vec == a)]][[which(rVals == r)]] <- data
    }
    names(ind_data[[which(age_vec == a)]]) <- sapply(rVals, FUN = function(x){paste0("rVal_", x)})
  }
  names(ind_data) <- sapply(age_vec, FUN = function(x){paste0("avg_age_", x)})
  return(ind_data)
  
  # This function returns a hierarchy of lists:
  # If nWF = FALSE, ind_data has 3 entries representing different 
  # R-values. Each of those entries is a list that contains 100 entries representing 
  # different replicates. Each of those entries is a dataframe with the following 
  # columns: "pedigree_id", "het", "gen", and "age" 
  # If nWF = TRUE, ind_data has 4 entries representing average ages 2, 5, 10, and 20. 
  # Each of those entries has 3 entries representing different 
  # R-values. Each of those entries is a list that contains 100 entries representing 
  # different replicates. Each of those entries is a dataframe with the following 
  # columns: "pedigree_id", "het", "gen", and "age" 
} # end function

extract_dist <- function(bins, rVals, age_vec, reps, column, niter){
  # testing vars
  # bins = pi_bins
  # age_vec <- ages[2:4]
  # rVals = rVals
  # reps = 100
  
  bootstraps <- vector(mode = "list", length = length(age_vec)) # list with one element per age (4)
  for (a in 1:length(age_vec)){
    hold_straps <- vector(mode = "list", length = length(rVals)) 
    for (r in 1:length(rVals)){
      for (i in 1:reps){
        data <- bins[[a]][[r]][[i]][,column]
        data <- gsub(",", "", data)
        data <- gsub("\\[", "", data)
        data <- gsub("\\]", "", data)
        straps <- matrix(nrow = length(data), ncol = niter)
        for(c in 1:length(data)){
          straps[c,] <- as.numeric(unlist(strsplit(data[c], split = " ")))
        }
        hold_straps[[r]][[i]] <- straps
      }
    }
    bootstraps[[a]] <- hold_straps
  }# process bootsraps
  return(bootstraps)
}

# Load datafiles into R------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

theta_bins <- load_ind_data(filePath = indir, model = "nWF", reps = reps, rVals = rVals, age_vec = ages, 
              keyword = "*bins.txt", separator = ",", index = FALSE)

# Process data into correct format------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

# for original run
# overall_theta_dist <- extract_dist(theta_bins, rVals = rVals, age_vec = ages, reps = 100, column = "overall_theta_dist", niter = 1e4)
# lower_theta_dist <- extract_dist(theta_bins, rVals = rVals, age_vec = ages, reps = 100, column = "lower_theta_dist", niter = 1e4)
# upper_theta_dist <- extract_dist(theta_bins, rVals = rVals, age_vec = ages, reps = 100, column = "upper_theta_dist", niter = 1e4)

# for theta bin test runs 
upper_ten <- extract_dist(theta_bins, rVals = rVals, age_vec = ages, reps = 100, column = "upper_ten", niter = 1e4)
lower_ten <- extract_dist(theta_bins, rVals = rVals, age_vec = ages, reps = 100, column = "lower_ten", niter = 1e4)
upper_five <- extract_dist(theta_bins, rVals = rVals, age_vec = ages, reps = 100, column = "upper_five", niter = 1e4)
lower_five <- extract_dist(theta_bins, rVals = rVals, age_vec = ages, reps = 100, column = "lower_five", niter = 1e4)
max_age <- extract_dist(theta_bins, rVals = rVals, age_vec = ages, reps = 100, column = "max_age", niter = 1e4)
min_age <- extract_dist(theta_bins, rVals = rVals, age_vec = ages, reps = 100, column = "min_age", niter = 1e4)


theta_sum <-  vector(mode = "list", length = length(ages))
# save summary stats for each rep 
for(a in 1:length(ages)){
  for(r in 1:length(rVals)){
    hold <- vector(mode = "list", length = 100)
    for(i in 1:100){
      hold[[i]] <- theta_bins[[a]][[r]][[i]][,c("timepoint", "upper_ten", "lower_ten", "upper_five", "lower_five", "max_age", "min_age")]
#       hold[[i]] <- theta_bins[[a]][[r]][[i]][,c("timepoint", "overall_pi", "lower_pi", "upper_pi", "overall_theta", "lower_theta", "upper_theta")]
      
    }
    theta_sum[[a]][[r]] <- hold
  }
}


# Save as R object------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
theta_results <- list(overall_theta_dist, lower_theta_dist, upper_theta_dist, theta_sum, date)

save(theta_results, file = paste0(outdir, "/bin_test_results_", date, ".Robj"))

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
