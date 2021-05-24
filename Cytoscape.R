# Title     : Cytoscape
# Objective : create a Cytoscape File
# Created by: Smaranda SOlomon
# Created on: 2/26/2021

#Filtered Cytoscape

to.Cytoscape <- function(all.data.filter, data){
  
all.data.cytoscape <- all.data.filter %>%
  mutate(compPASS = TRUE) %>%
  select(Bait, Prey, compPASS) %>%
  right_join(data, by = c("Bait", "Prey")) %>%
  replace_na(list(compPASS = FALSE))%>%
  select(Bait, Prey, compPASS)


cytoscape <- data %>%
  mutate(saint = ifelse(BFDR <= 0.05 & AvgP >= 0.7, TRUE, FALSE)) %>%
  right_join(all.data.cytoscape, by = c("Bait", "Prey")) %>%
  select(Bait.Gene.Name, Prey.Gene.Name, saint, compPASS, in_BioGRID, is_Bait, Nucleus, Cytoplasm, Cytoskeleton, Endosome, ER, 
         Extracellular, Golgi, Lysosome, Mitochondria, Peroxisome, Plasma_Membrane, Vesicles, Cell_Junction, class_2019)

#cytoscape <- cytoscape %>%
#  dplyr::filter(grepl('CLK3|CLK4|DYRK2|LMTK2|PIP5K1A|PRPF4B|RIOK1|STK19|STK17A', Bait.Gene.Name))

}


################################################################################
#Create file for animated Cytoscape interaction network

#baits <- unique((data %>% filter(!is.na(BaitGene), is_Bait == TRUE, BaitGene == PreyGene))$Bait)
#all.data1 <- data[0,]

#for (mybait in baits) {
#  bait.data <- data %>% filter(Bait == mybait)
#  num.interactors <- min(max(10, nrow(bait.data)*0.05), nrow(bait.data))
#  all.data1 <- all.data1 %>% add_row(bait.data[1:num.interactors, ])
#}

#all.data1 <- all.data1 %>%
#  mutate(compPASS = TRUE) %>%
#  select(Bait, Prey, compPASS) %>%
#  right_join(data, by = c("Bait", "Prey")) %>%
#  replace_na(list(compPASS = FALSE))%>%
#  select(Bait, Prey, compPASS)

#cytoscape <- data %>%
#  mutate(saint = ifelse(BFDR <= 0.05 & AvgP >= 0.7, TRUE, FALSE)) %>%
#  right_join(all.data1, by = c("Bait", "Prey")) %>%
#  select(Bait.Gene.Name, Prey.Gene.Name, saint, compPASS, in_BioGRID, is_Bait)

#cytoscape <- cytoscape %>%
#  dplyr::filter(grepl('CLK3|CLK4|DYRK2|LMTK2|PIP5K1A|PRPF4B|RIOK1|STK19|STK17A', Bait.Gene.Name))


#Write file with filtered experiments, only Dark kinases
#write_csv(cytoscape, 'output/Animated_cytoscape_poster.csv')