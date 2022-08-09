library(dplyr)


to_comp_test = read.csv(file = "C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_int/to_CompPASS.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


compute_wds = function(data,method="spc") {
  # Data should have Biat, Prey, Replicate, and Spectral Counts columns
  
  # First define K: the number of experimental conditions
  k = length(unique(data[,c("Bait","Replicate")])[[1]])
  
  # Next, do bait-level calculations...
  
}

compute_wds(to_comp_test)