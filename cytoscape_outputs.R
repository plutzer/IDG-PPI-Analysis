
library(reticulate)
library(cRomppass)
library(tidyverse)
library(dplyr)


create_cytoscape_outputs = function (annotated,annotated_filtered,mv_biogrid_path,output_dir) {
  #Load in the biogrid file
  mv_biogrid = read.csv(file = mv_biogrid_path, header = TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606)
  
  # Get a list of baits
  baits <- unique((annotated %>% filter(!is.na(BaitGene), is_Bait == TRUE, BaitGene == PreyGene))$Bait)
  
  # Get a list of genes from data
  geneList = as.vector(annotated$`First.Prey.GeneID`[!is.na(annotated$`First.Prey.GeneID`)])
  
  interactions <- mv_biogrid %>%
    # select columns needed from bioGrid
    select(Entrez.Gene.Interactor.A, Entrez.Gene.Interactor.B) %>%
    # only take proteins with geneIDs
    filter(Entrez.Gene.Interactor.A != "-" & Entrez.Gene.Interactor.B != "-") %>% 
    # no self interactions
    filter(Entrez.Gene.Interactor.A != Entrez.Gene.Interactor.B) %>%
    # always putting the smaller geneID first
    mutate(interactor_min = pmin(as.numeric(Entrez.Gene.Interactor.A), as.numeric(Entrez.Gene.Interactor.B)), 
           interactor_max = pmax(as.numeric(Entrez.Gene.Interactor.A), as.numeric(Entrez.Gene.Interactor.B))) %>%
    select(interactor_min, interactor_max) %>%
    unique()
  
  # Create base nodes and edges files.
  nodes_combined = select(annotated,Prey.Gene.Name,class,is_Bait,WD,perm_fdrs,FoldChange,GO.Slim,Prey.Gene.Synonym,First.Bait.GeneID,First.Prey.GeneID,Bait.Gene.Name)
  edges_combined = select(annotated,Experiment.ID,Prey.Gene.Name,in.BioGRID,Multivalidated,Evidence.Weight,Experimental.Systems,Authors,Publications,WD,perm_fdrs,Bait.Gene.Name)
  nodes_combined_filtered = select(annotated_filtered,Prey.Gene.Name,class,is_Bait,WD,perm_fdrs,FoldChange,GO.Slim,Prey.Gene.Synonym,First.Bait.GeneID,First.Prey.GeneID,Bait.Gene.Name)
  edges_combined_filtered = select(annotated_filtered,Experiment.ID,Prey.Gene.Name,in.BioGRID,Multivalidated,Evidence.Weight,Experimental.Systems,Authors,Publications,WD,perm_fdrs,Bait.Gene.Name)
  
  # Add the interaction column
  edges_combined$Interaction = "Strep APMS"
  edges_combined_filtered$Interaction = "Strep APMS"
  
  # Create cytoscape directory
  dir.create(paste(output_dir,"/Cytoscape_Outputs",sep=''))
  
  #Write the combined nodes files... combined edges file is too big.
  write_csv(nodes_combined, paste(output_dir,"/Cytoscape_Outputs/Combined_nodes.csv",sep=''))
  write_csv(nodes_combined_filtered, paste(output_dir,"/Cytoscape_Outputs/Combined_nodes_filtered.csv",sep=''))
  
  for (mybait in baits) {
    
    bait.data.filter <- annotated %>% filter(Bait == mybait)
    num.interactors <- min(max(10, nrow(bait.data.filter)*0.05), nrow(annotated_filtered))
    bait.data.filter.comp <- bait.data.filter[1:num.interactors, ]
    # all.data.filter <- all.data.filter %>% add_row(bait.data.filter.comp)
    
    #Filtering for interactions
    prey.prey.inter <- filter(interactions, (`interactor_min` %in%  bait.data.filter.comp$First.Prey.GeneID), (`interactor_max` %in%  bait.data.filter.comp$First.Prey.GeneID)) %>%
      rename(c(interactor_min = "Prey.1.Entrez.ID", interactor_max = "Prey.2.Entrez.ID")) %>% #Changing column names
      mutate(Prey.1.Entrez.ID = as.character(Prey.1.Entrez.ID), 
             Prey.2.Entrez.ID = as.character(Prey.2.Entrez.ID)) #Changing data type from double to character to be left_joined with all.data.filter
    
    #Left_joining table with prey 1 and adding Uniprot and Nice Prey name columns
    prey.prey.join <- left_join(prey.prey.inter, bait.data.filter.comp, by=c("Prey.1.Entrez.ID" = "Prey.GeneID"))%>%
      select(Prey.1.Entrez.ID, Prey.2.Entrez.ID, Canonical.First.Prey.Uniprot, Prey.Gene.Name, WD) %>%
      rename(c(Canonical.First.Prey.Uniprot = "Prey.1.Uniprot", Prey.Gene.Name = "Prey.1.Gene.Name", WD = "WD.1"))
    
    #Left_joining table with prey 2 and adding Uniprot and Nice Prey name columns
    prey.prey.final <- left_join(prey.prey.join, bait.data.filter.comp, by=c("Prey.2.Entrez.ID" = "Prey.GeneID")) %>%
      select(Prey.1.Gene.Name, Prey.1.Uniprot, Prey.1.Entrez.ID, WD.1, Prey.Gene.Name, Canonical.First.Prey.Uniprot, Prey.2.Entrez.ID, WD,perm_fdrs) %>%
      rename(c(Canonical.First.Prey.Uniprot = "Prey.2.Uniprot", Prey.Gene.Name = "Prey.2.Gene.Name", WD = "WD.2")) %>%
      unique()
    
    ###############################################################################
    ##############################################################################
    ##############################################################################
    ##############################################################################
    if(nrow(prey.prey.final) != 0){
      prey.prey.final$Source <- 'Biogrid'
    }
    else {
      print("No prey-prey!")
      prey.prey.final$Source <- character(0)
      
    }
    
    # Check to see if its a contol bait
    if (mybait %in% test_baits$Bait.ID) {
      
      bait = filter(bait.data.filter,Bait == mybait)$Bait.Gene.Name[[1]]
      
      print(bait)######
      
      dir.create(paste(output_dir,"/Cytoscape_Outputs/",bait,sep=''))
      nodes_combined %>%
        filter(Bait.Gene.Name == bait) %>%
        write.csv(paste(output_dir,"/Cytoscape_Outputs/",bait,"/",bait,"_nodes.csv",sep=''),na="",row.names=F)
      nodes_combined_filtered %>%
        filter(Bait.Gene.Name == bait) %>%
        write.csv(paste(output_dir,"/Cytoscape_Outputs/",bait,"/",bait,"_nodes_filtered.csv",sep=''),na="",row.names=F)
      
      edges_bait = edges_combined %>%
        filter(Bait.Gene.Name == bait)
      edges_bait_filtered = edges_combined_filtered %>%
        filter(Bait.Gene.Name == bait)
      
      # Merge the prey-prey information with the nodes
      prey.prey.final$Prey.1.Entrez.ID = as.character(prey.prey.final$Prey.1.Entrez.ID)
      prey.prey.final$Prey.2.Entrez.ID = as.character(prey.prey.final$Prey.2.Entrez.ID)
      
      
      prey.prey.final.biogrid = left_join(prey.prey.final,summ_biogrid, by=c('Prey.1.Entrez.ID'='Entrez.Gene.Interactor.A','Prey.2.Entrez.ID'='Entrez.Gene.Interactor.B'))
      
      # Add the relevant columns to the edges dataframe
      prey.prey.select = select(prey.prey.final.biogrid,Prey.1.Gene.Name,in.BioGRID,Multivalidated,Evidence.Weight,Experimental.Systems,Authors,Publications,WD.1,Prey.2.Gene.Name,Source,WD.2,perm_fdrs)
      
      edges_bait$X = NULL
      
      # What if prey.prey.select has 0 rows?
      if (nrow(prey.prey.select) >= 1){
        if (length(edges_bait)>=1){
          prey.prey.select$Experiment.ID = edges_bait$Experiment.ID[[1]] # These shouldn't change so just first is fine
        }
        else {
          prey.prey.select$Experiment.ID = '-'
        }
      } else {
        prey.prey.select = cbind(prey.prey.select,data.frame(Experiment.ID = character(),stringsAsFactors = F))
      }
      print(colnames(prey.prey.select))
      # Put Experiment ID at the start of the dataframe
      prey.prey.select = prey.prey.select[, c(13,1:10,12,11)]
      edges_bait = edges_bait[, c(1:9,11,12,10)]
      edges_bait_filtered = edges_bait_filtered[, c(1:9,11,12,10)]
      
      # Rename the columns for compatibility...
      colnames(prey.prey.select) = colnames(edges_bait)
      
      # Combine and add the WD2 column
      total_edges = bind_rows(edges_bait,prey.prey.select)
      total_edges_filtered = bind_rows(edges_bait_filtered,prey.prey.select)
      colnames(total_edges)[[13]] = "WD2"
      colnames(total_edges_filtered)[[13]] = "WD2"
      
      
      
      write.csv(total_edges, paste(output_dir,"/Cytoscape_Outputs/",bait,"/",bait,"_edges.csv",sep=''),na="",row.names=F)
      write.csv(total_edges_filtered, paste(output_dir,"/Cytoscape_Outputs/",bait,"/",bait,"_edges_filtered.csv",sep=''),na="",row.names=F)
    }
  }
}