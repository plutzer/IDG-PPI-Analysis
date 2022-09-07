
library(tidyverse)
library(dplyr)


filter_interaction = function(pathname) {
  # Pathname = path of interaction file to filter.
  # Filters interaction file to remove preys with 0 intensity in all replicates of a given bait
  
  interaction = read.delim(file = pathname,sep = '\t',header = F)
  
  interaction_filtered = interaction %>% rename(
    rep = V1,
    bait = V2,
    prey = V3,
    intensity = V4
  ) 
  
  interaction_sums = interaction_filtered %>% group_by(bait,prey) %>%
    summarise(sum_int = sum(intensity),.groups="drop")
  
  interaction_filtered = left_join(interaction_filtered,interaction_sums,by=c("bait","prey")) %>% 
    filter(sum_int != 0) %>% 
    subset(select = -c(sum_int))
  
  write.table(interaction_filtered,file = paste0(str_split(pathname)[[1]][1],"_filtered.txt"),quote = F,sep="\t",row.names = F,col.names = F)
  paste0(str_split(pathname)[[1]][1],"_filtered.txt")
}

