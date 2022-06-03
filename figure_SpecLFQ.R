# Title     : figure_SpecLFQ
# Objective : Compare Final output when using spectral counts vs LFQ intensities
# Created by: Smaranda Solomon
# Created on: 5/23/2022

library("venneuler")
library(ggplot2)
library(tidyverse)

spec <- read.csv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Results/20220523_103_Baits/Spectral/Annotated_Merge_All_filtered_DB1N_spec.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)
lfq <- read.csv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Results/20220523_103_Baits/LFQ/Annotated_Merge_All_filtered_DB1N_lfq.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

# plot(venneuler(c(A = 0,          # Draw pairwise venn diagram
#                  B = 0,
#                  "A&B" = 10)))

#Plot comparing all WD scores against each other, no regards to Bait
# spec.wd.full <- spec$WD
# lfq.wd.full <- lfq$WD
# 
# plot(spec.wd.full, lfq.wd.full, xlab = "Spectral", ylab = "LFQ", pch = 16)
# abline(lm(spec.wd.full ~ lfq.wd.full))


#Plot comparing all preys at once
spec2 <- spec %>%
  dplyr::select(Bait.Gene.Name, Prey.Gene.Name, WD) %>%
  rename(WD.Spectral = WD)

lfq2 <- lfq %>%
  dplyr::select(Bait.Gene.Name, Prey.Gene.Name, WD) %>%
  rename(WD.LFQ = WD)

#Combine the two
comb <- full_join(spec2, lfq2, by = c("Bait.Gene.Name", "Prey.Gene.Name"))
  

comb <- comb[order(comb$Bait.Gene.Name), ]



#Just for me to know the number of unique baits
# preyNum <- comb %>% 
#   select(Bait.Gene.Name) %>%
#   unique()

write.csv(comb, 'C:/Users/smaranda/Documents/SmarandaSolomon/Results/20220523_103_Baits/combined_diff_prey.csv')

plots <- ggplot(comb, aes(log2(WD.Spectral), log2(WD.LFQ))) + geom_point()
plots + facet_wrap(~Bait.Gene.Name, ncol = 12)
