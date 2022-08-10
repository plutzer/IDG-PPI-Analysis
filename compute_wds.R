library(dplyr)


to_comp_test = read.csv(file = "C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_int/to_CompPASS.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

power <- function(x, y) sign(x) * abs(x)^y


compute_wds = function(data,quant_col="Spectral.Count",num_reps=2) {
  # Data should have Biat, Prey, Replicate, and Spectral Counts columns
  
  # First define K: the number of experimental conditions
  k = length(unique(data[,c("Bait")])[[1]])*num_reps

  # Need to pre-process the quantification to be the average of the reps (since some reps will be missing)
  data = data %>% group_by(Bait,Prey) %>% summarize(Quantification = sum(.data[[quant_col]])/num_reps,.groups="drop")
  
  # compute f
  data$f = as.integer(data$Quantification >= 1)
  
  # Prey-level calculations
  p = data %>% count(Bait,Prey)
  
  prey_stats = data %>% group_by(Prey) %>% summarise(Gamma = (k/sum(f)),
                                                     Xbar = sum(Quantification)/k,
                                                     stdev = sd(Quantification)) %>% ungroup() %>%
                        mutate(Omega = stdev/Xbar)

  # prey_stats$Xbar = (data %>% group_by(Prey) %>% summarise(xbar = sum(.data[[quant_col]])/k))$xbar
  # prey_stats$Omega = ((data %>% group_by(Prey) %>% summarise(stdev = sd(.data[[quant_col]])))$stdev)/prey_stats$Xbar
  
  # WD score calculation  
  data = data %>% left_join(prey_stats,by="Prey") %>% left_join(p,by=c("Bait","Prey"))
  
  data$WD = sqrt(`^`((data$Gamma*data$Omega),data$n)*(data$Quantification))
  
  #mutate(WD = sqrt((`^`((Gamma*Omega),p))*(Quantification)))
  data
}


###############################
# TESTING

comp_out = comppass(to_comp_test,stats = NULL,norm.factor = 0.98)

my_out = compute_wds(to_comp_test)



# Testing to see if the outputs agree
merged = full_join(comp_out,my_out,by=c("Bait","Prey"))
# There's some seriously funky shit goin on here
# WD scores don't agree

# The length of the output is different. The comppass command seems to retain more rows than unique Bait-Prey combos
print(paste("Length of my output: ",length(my_out[[1]])))
print(paste("Length of comppass output: ",length(comp_out[[1]])))
print(paste("Length of unique Bait-Prey combos: "),length(merged %>% count(Bait,Prey)[[1]]))

# What bait-prey combos are in the default that are not in my implementation


raw = to_comp_test %>% subset(select = -Prey.Name)

write.table(raw,file="C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_int/raw_to_comp.tsv",sep="\t")

