# Title     : Cytoscape_REVISED
# Objective : create a Cytoscape File, REVISED 4/10/2022
# Created by: Smaranda Solomon
# Created on: 4/10/2022

#to.Cyto <- function(st, fn){
  
data <- read.csv(file='output/Annotated_Merge_All_filtered.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Add all columns needed for Node file
node <- data %>%
  dplyr::select(Bait.Gene.Name, First.Bait.GeneID, Prey.Gene.Name, First.Prey.GeneID, FoldChange, AvgSpec, 
         Nucleus, Cytoplasm, Cytoskeleton, Endosome, ER, 
         Extracellular, Golgi, Lysosome, Mitochondria, Peroxisome, Plasma_Membrane, Vesicles, Cell_Junction,
         class_2019, is_Bait) %>%
  drop_na(Prey.Gene.Name, Bait.Gene.Name) #Remove NA's in Bait and Prey)

node <- node[!duplicated(node[,c("Bait.Gene.Name","Prey.Gene.Name")]),] #Remove duplicate Bait/Prey pairs

#Add all columns needed for Edge file
edge <- data %>%
  dplyr::select(Bait.Gene.Name, Prey.Gene.Name, Prey.Gene.Synonym, in_BioGRID, SaintScore, BFDR,
         WD, Little.N) %>%
  drop_na(Prey.Gene.Name, Bait.Gene.Name) #Remove NA's in Bait and Prey)
  
#}
