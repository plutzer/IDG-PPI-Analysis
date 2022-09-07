library(dplyr)
library(cRomppass)
library(magrittr)


normalize.wd <- function(xs, norm.factor) {
  # Normalizes a vector of WD scores
  
  xs.f <- Filter(function(x) { ! (is.nan(x) || is.na(x)) }, xs)
  return(xs / quantile(xs.f, norm.factor)[1])
}



resample_AvePSM = function(to_comp,comp_out,n.experiments,norm.factor = 0.98,suffix = "",test = TRUE) {
  # Resamples the data once and computes new WD scores, adding them to the comp_out dataframe
  
  # Permute the intensities within each Bait
  perm = to_comp[order(to_comp$Bait),]
  perm$new_spec = perm %>% group_by(Bait) %>% group_map(~ sample(.x$Spectral.Count,replace = TRUE)) %>% unlist()
  
  # If this is not a test, change the spectral counts to be the permutation
  if (!test) {
    perm$Spectral.Count = perm$new_spec
  }
  
  # Calculate summary statistics for the permuted data
  stats = perm %>%
    group_by(Experiment.ID, Prey, Replicate) %>%
    summarize(Bait = unique(Bait),
              Spectral.Count = max(Spectral.Count),.groups="drop") %>%
    ungroup() %>%
    group_by(Experiment.ID, Prey) %>%
    summarize(Bait = unique(Bait),
              AvePSM = mean(Spectral.Count),
              N.Saw = length(which(Spectral.Count > 0)),.groups="drop") %>%
    ungroup() %>%
    as.data.frame()
  
  # Add Little N to stats df
  little_n_entries = comp_out %>% group_by(Prey) %>% summarize(Little.N = unique(Little.N),.groups="drop")
  stats = stats %>% left_join(little_n_entries,by="Prey")
  
  # Calculate prey statistics needed for WD scores
  preystats = stats %>% filter(Bait != Prey) %>% group_by(Prey) %>%  summarize(
    prey.mean = sum(AvePSM)/n.experiments,
    n.experiments.without.prey = n.experiments - unique(Little.N),
    prey.sum.squared.err = sum(raise_to_power(AvePSM - prey.mean, 2)) +
      (raise_to_power(prey.mean, 2) * n.experiments.without.prey),.groups="drop") %>% ungroup() %>%
    mutate(prey.sd = sqrt(prey.sum.squared.err / (n.experiments - 1)))
  
  # Edge case for preys that are only seen with themselves as bait
  missing_preys = setdiff(unique(stats$Prey),unique(preystats$Prey))
  missing_preystats = stats %>% filter(Prey %in% missing_preys) %>% group_by(Prey) %>% summarize(
    prey.mean = sum(AvePSM)/n.experiments,
    n.experiments.without.prey = n.experiments,
    prey.sum.squared.err = sum(raise_to_power(AvePSM-prey.mean, 2)) + 
      (raise_to_power(prey.mean, 2) * n.experiments.without.prey),.groups="drop") %>% ungroup() %>%
    mutate(prey.sd = sqrt(prey.sum.squared.err / (n.experiments - 1)))
  preystats = rbind(preystats,missing_preystats)
  
  # Merge prey stats with the original stats
  stats = left_join(stats,preystats,by="Prey")
  
  # Make the WD_score calculation
  stats$N.Saw = as.numeric(stats$N.Saw)
  stats$inner_term = (n.experiments/stats$Little.N)*(stats$prey.sd/stats$prey.mean)
  stats$WD = sqrt(stats$AvePSM*raise_to_power(stats$inner_term, stats$N.Saw))
  
  # Normalize the WD scores
  stats$WD = normalize.wd(stats$WD, norm.factor)
  
  # Rename columns and subset for merge
  stats = stats %>% dplyr::select(Experiment.ID,Prey,AvePSM,WD) %>% rename_with(test, .fn = ~paste(.,suffix,sep="_"), .cols = c(AvePSM,WD))
  
  # Merge with the comp_out input
  comp_out = merge(comp_out,stats,by=c("Experiment.ID","Prey"))
  
  # return the comp_out input with new columns
  comp_out
}


run_comppass = function (to_comp_filename,n_iter = 100) {
  
  # Read in the input file
  to_comp = read.csv(file = to_comp_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Produce the starting output
  comp_out = comppass(to_comp, stats = NULL, norm.factor = 0.98)
  
  ###### Permutation FDR by resampling #############
  print("Starting CompPASS resampling...")
  pb = txtProgressBar(min = 0, max = n_iter, initial = 0) 
  
  # Initialize the output with a test run. Allows you to double check that the WD calculation is correct for this data
  comp_resample = resample_AvePSM(to_comp,comp_out,n.experiments = length((to_comp %>% count(Experiment.ID))[[1]]),suffix = "test",test=T)
  
  # Resample
  for (i in 1:n_iter) {
    comp_resample = resample_AvePSM(to_comp,comp_resample,n.experiments = length((to_comp %>% count(Experiment.ID))[[1]]),suffix = as.character(i),test=F)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  # Write resampling results for troubleshooting
  write.table(comp_resample,file=paste(output_dir,'/compPASS_resampled.csv',sep=''),sep="\t")
  
  wds = comp_resample %>% select(contains("WD_") & !contains("test"))
  
  # Compute the false discovery rate
  comp_out$perm_fdrs = apply(as.array(1:length(comp_resample$WD)), MARGIN = 1, FUN = function(i) {sum(wds[i,]>comp_resample$WD[[i]])/resampling_iterations})
  
  # Save the final output
  write.table(comp_out, file=paste(output_dir,'/compPASS.csv',sep=''), sep = "\t")
  
  # Return the output
  comp_out
}

