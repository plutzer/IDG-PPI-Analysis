# Title     : Merge_CompPASS_SAINT
# Objective : Merge the results from CompPASS and SAINT into one table
# Created by: Smaranda
# Created on: 8/30/2020

library("tidyverse")

comp <- read.csv(file='compPASS.csv', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
saint <- read.csv(file='list.txt', sep="\t")

merg = left_join(comp, saint, by = c("Experiment.ID" = "Bait", "Prey"))

write_csv(merg, 'Merge_CompPASS_SAINT.csv')

