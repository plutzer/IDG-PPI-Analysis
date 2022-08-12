library(dplyr)


to_comp_test = read.csv(file = "C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_int/to_CompPASS.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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


