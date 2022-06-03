# Title     : Merge_CompPASS_SAINT
# Objective : Merge the results from CompPASS and SAINT into one table
# Created by: Smaranda
# Created on: 8/30/2020

library("tidyverse")

comp <- read.csv(file='output/compPASS.csv', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
saint <- read.csv(file='output/list.txt', sep="\t")

merg <- left_join(comp, saint, by = c("Experiment.ID" = "Bait", "Prey")) 

merg <- merg %>%
  drop_na(SaintScore)

write_csv(merg, 'output/Merge_CompPASS_SAINT.csv')

