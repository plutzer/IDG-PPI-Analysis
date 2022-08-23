library(dplyr)
library(cRomppass)
library(magrittr)


to_comp_test = read.csv(file = "C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_int/to_CompPASS.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

comp_out = comppass(to_comp_test,stats = NULL,norm.factor = 0.98)

normalize.wd <- function(xs, norm.factor) {
  xs.f <- Filter(function(x) { ! (is.nan(x) || is.na(x)) }, xs)
  return(xs / quantile(xs.f, norm.factor)[1])
}

resample_AvePSM = function(to_comp,comp_out,n.experiments,norm.factor = 0.98,suffix = "",test = TRUE) {
  # Runs in about 7 seconds
  
  perm = to_comp[order(to_comp_test$Bait),]
  perm$new_spec = perm %>% group_by(Bait) %>% group_map(~ sample(.x$Spectral.Count,replace = TRUE)) %>% unlist()
  
  # If this is not a test, change the spectral counts to be the permutation
  if (!test) {
    perm$Spectral.Count = perm$new_spec
  }
  
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
  
  # Need to add Little N to stats df
  little_n_entries = comp_out %>% group_by(Prey) %>% summarize(Little.N = unique(Little.N),.groups="drop")
  
  stats = stats %>% left_join(little_n_entries,by="Prey")
  
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
  
  # Merge with the original stats
  stats = left_join(stats,preystats,by="Prey")
  
  stats$N.Saw = as.numeric(stats$N.Saw)
  
  stats$inner_term = (n.experiments/stats$Little.N)*(stats$prey.sd/stats$prey.mean)
  stats$WD = sqrt(stats$AvePSM*raise_to_power(stats$inner_term, stats$N.Saw))
  
  stats$WD = normalize.wd(stats$WD, norm.factor) 
  stats = stats %>% select(Experiment.ID,Prey,prey.mean,prey.sd,AvePSM,WD,N.Saw,Little.N) %>% rename_with(test, .fn = ~paste(.,suffix,sep="_"), .cols = c(AvePSM,WD,prey.mean,prey.sd,N.Saw,Little.N))
  
  comp_out = merge(comp_out,stats,by=c("Experiment.ID","Prey"))
  
  comp_out
}

# Test and time...
start = Sys.time()
res = resample_AvePSM(to_comp_test,comp_out,n.experiments = length((to_comp_test %>% count(Experiment.ID))[[1]]),suffix = 1)
end = Sys.time()
start-end


