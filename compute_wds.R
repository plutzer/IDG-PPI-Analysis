library(dplyr)
library(cRomppass)


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
  
  # Merge with the original stats
  stats = left_join(stats,preystats,by="Prey") %>% mutate(
    WD = sqrt(AvePSM*raise_to_power((n.experiments*prey.sd/(Little.N*prey.mean)),N.Saw))
  )
  
  stats$WD = normalize.wd(stats$WD, norm.factor) 
  stats = stats %>% select(Experiment.ID,Prey,AvePSM,WD) %>% rename_with(test, .fn = ~paste(.,suffix,sep="_"), .cols = c(AvePSM,WD))
  
  
  comp_out = merge(comp_out,stats,by=c("Experiment.ID","Prey"))
  
  comp_out
}

# Test and time...
start = Sys.time()
res = resample_AvePSM(to_comp_test,comp_out,n.experiments = length((to_comp_test %>% count(Experiment.ID))[[1]]),suffix = 1)
end = Sys.time()
start-end

# Check to make sure the result is the same...
merged = left_join(comp_out,res,by=c("Experiment.ID","Prey"))























n.experiments = length((to_comp_test %>% count(Experiment.ID))[[1]])

teststats = res %>% filter(Bait != Prey) %>% group_by(Prey) %>%  summarize(
  prey.mean = sum(AvePSM)/n.experiments,
  n.experiments.without.prey = n.experiments - unique(Little.N),
  prey.sum.squared.err <- sum(raise_to_power(AvePSM - prey.mean, 2)) +
    (raise_to_power(prey.mean, 2) * n.experiments.without.prey),.groups = NULL) %>%
  ungroup()










data = to_comp_test

quant_col="Spectral.Count"

k = length((data %>% count(Experiment.ID))[[1]]) # Total # of different experimental conditions

preystats = data %>% group_by(Prey) %>% filter(Bait != Prey) %>% summarize(stdev = sd(.data[[quant_col]]),
                                                  Xbar = sum(.data[[quant_col]])/k)

combostats = data %>% group_by(Bait,Prey,Experiment.ID) %>% summarize(AvgSpec = mean(.data[[quant_col]]),
                                                                      p = n(),
                                                                      .groups="drop")

combostats$f = as.integer(combostats$AvgSpec >= 1)

tallys = combostats %>% group_by(Prey) %>% summarize(total_ints = sum(f),
                                                     avg_stdev = sd(AvgSpec),
                                                     avg_Xbar = mean(AvgSpec),
                                                     .groups="drop")

stats = combostats %>% left_join(preystats,by=c("Prey")) %>% left_join(tallys,by=c("Prey"))

stats = stats %>% mutate(WD = sqrt(`^`((k/total_ints)*(stdev/Xbar),p)*AvgSpec))

# Test merging to see if they agree...
comp_out = comppass(to_comp_test,stats = NULL,norm.factor = 0.98)

merge = left_join(comp_out,stats,by=c("Bait","Prey","Experiment.ID"))








################################################################################
################################################################################
################################################################################
################################################################################

# Data should have Biat, Prey, Replicate, and Spectral Counts columns

# First define K: the number of experimental conditions
k = length((data %>% count(Experiment.ID))[[1]]) # Correct
#k = length(unique(data[,c("Bait")])[[1]])*num_reps

p = data %>% count(Bait,Prey,Experiment.ID) # Correct

data = data %>% left_join(p,by=c("Bait","Prey","Experiment.ID"))

# Need to pre-process the quantification to be the average of the reps (since some reps will be missing)
quantdata = data %>% group_by(Bait,Prey,Experiment.ID) %>% summarize(Quantification = mean(.data[[quant_col]]),.groups="drop")

data = data %>% left_join(quantdata,by=c("Bait","Prey","Experiment.ID"))

# compute f
data$f = as.integer(data$Quantification >= 1) #Correct

data = data %>% left_join(reps_avg_data,by=c("Bait","Prey","Experiment.ID"))

prey_stats = data %>% group_by(Prey) %>% summarize(Gamma = (k/sum(f)),
                                                   Xbar = sum(.data[[quant_col]])/k,
                                                   stdev = sd(.data[[quant_col]])) %>% ungroup() %>%
  mutate(Omega = stdev/Xbar)

# WD score calculation  
data = data %>% left_join(prey_stats,by="Prey") %>% left_join(p,by=c("Bait","Prey"))

data$WD = sqrt(`^`((data$Gamma*data$Omega),data$n)*(data$Quantification))

#mutate(WD = sqrt((`^`((Gamma*Omega),p))*(Quantification)))



###############################
# TESTING

comp_out = comppass(to_comp_test,stats = NULL,norm.factor = 0.98)

my_out = data



# Testing to see if the outputs agree
merged = left_join(comp_out,(my_out),by=c("Bait","Prey"))


