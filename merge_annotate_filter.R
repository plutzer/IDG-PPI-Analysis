
library(org.Hs.eg.db)
library(reticulate)
library(cRomppass)
library(tidyverse)
library(DarkKinaseTools)
library(dplyr)

merge_scores = function(SAINT_list_path,comp_out,output_dir) {
  # Read in output list from SAINT
  saint <- read.csv(file=SAINT_list_path, sep="\t")

  # Merge with comppass output
  merge = left_join(as.data.frame.matrix(comp_out), saint, by = c("Experiment.ID" = "Bait", "Prey"))
  
  # Write output and return the dataframe
  write_csv(merge, paste(output_dir,'/Merge_CompPASS_SAINT.csv',sep=''))
  return(merge)
}

annotate_uniprot_go = function(merge,uniprot_map_path) {
  # Function adapted from Smaranda's code
  
  uniprot.mapping <- read_tsv(uniprot_map_path)
  
  # List of unique baits
  unique.Baits <- unique(merge$Bait)
  
  #Prey
  # extract Prey Uniprot Identifier
  data <- separate(data = merge, col = `Prey`, into = c("First.Prey.Uniprot"), sep=";", remove=F, extra="drop")
  
  # Remove isoform from the Prey Uniprot Identifier
  data <- separate(data, "First.Prey.Uniprot", c("Canonical.First.Prey.Uniprot"), sep="-", remove=F, extra="drop")
  
  # fill Prey Entrez GeneID
  data <- left_join(data, uniprot.mapping, by=c("Canonical.First.Prey.Uniprot" = "UniProt"))
  
  #Bait
  ## extract Bait Uniprot Identifier
  data <- separate(data = data, col = `Bait`, into = c("First.Bait.Uniprot"), sep=";", remove=F, extra="drop")
  
  # Remove isoform from the Bait Uniprot Identifier
  data <- separate(data, "First.Bait.Uniprot", c("Canonical.First.Bait.Uniprot"), sep="-", remove=F, extra="drop")
  
  # fill Bait Entrez GeneID
  data <- left_join(data, uniprot.mapping, by=c("Canonical.First.Bait.Uniprot" = "UniProt"))
  
  #Move and rename columns
  data <- data[c(1:7, 32:39, 8:31)]
  data <- dplyr::rename(data, c("Prey.Gene.Name" = "Gene_Name.x", "Prey.Gene.Synonym" = "Gene_Synonym.x", "Prey.GeneID" = "GeneID.x", "First.Prey.GeneID" = "First_GeneID.x" ,
                                "Bait.Gene.Name" = "Gene_Name.y", "Bait.Gene.Synonym" = "Gene_Synonym.y", "Bait.GeneID" = "GeneID.y", "First.Bait.GeneID" = "First_GeneID.y"))
  data <- data[,c(1,12,13,2:4,14,15,8,9,5:7,10,11,16:39)]
  
  ################################################################################
  
  #Annotations
  xx <- as.list(org.Hs.egGO2ALLEGS)
  x <-c("GO:0005634","GO:0005737","GO:0005856","GO:0005768","GO:0005783","GO:0005576","GO:0005794",
        "GO:0005764","GO:0005739","GO:0005777","GO:0005886","GO:0031982","GO:0005911")
  go <- c("Nucleus","Cytoplasm","Cytoskeleton","Endosome","ER","Extracellular","Golgi","Lysosome",
          "Mitochondria","Peroxisome","Plasma_Membrane","Vesicles","Cell_Junction")
  
  # 
  for(i in 1:length(go)) {
    goids <- xx[x[i]]
    goslim <- tibble(Prey.GeneID = goids[[1]], Evidence = names(goids[[1]]))
    #goslim <- filter(goslim, "IEA" != Evidence)
    goslim <- (goslim %>% 
                 group_by(Prey.GeneID)%>% 
                 summarise(n=n()))
    data <- left_join(data, goslim, by=c("Prey.GeneID" = "Prey.GeneID")) %>% 
      mutate(
        goNames = !is.na(n) 
      ) %>% 
      select(-`n`) %>% 
      rename(goNames = go[i])
  }
  
  #Merge GO Slim in Single Column
  data$GO.Slim <- apply(data, 1, 
                        function(x) {
                          str_c(go[as.logical(c(x[["Nucleus"]],x[["Cytoplasm"]],x[["Cytoskeleton"]], x[["Endosome"]], x[["ER"]],
                                                x[["Extracellular"]], x[["Golgi"]], x[["Lysosome"]], x[["Mitochondria"]],
                                                x[["Peroxisome"]], x[["Plasma_Membrane"]], x[["Vesicles"]], x[["Cell_Junction"]]))], 
                                collapse = ";")
                        }
  )
  
  #Annotate Dark Kinases
  data <- left_join(data, all_kinases, by=c("First.Prey.GeneID" = "entrez_id"))
  
  #Is it a Bait Column
  data$is_Bait <- apply(data, 1, function(x) {
    any(str_detect(unique.Baits, x[["Canonical.First.Prey.Uniprot"]]))
  }
  )
  
  #Nice bait name column addition by grabbing from prey name if it pairs with is_bait column
  baitTable <- data %>%
    filter(is_Bait == TRUE) %>%
    filter(str_detect(Bait, Prey)) %>%
    select(Bait, BaitGene = PreyGene) %>%
    distinct() 
  
  data <- left_join(data, baitTable, by="Bait")
  
  return(data)
}


preprocess_biogrid = function(biogrid_all_path, biogrid_mv_path,output_dir) {
  all_biogrid = read.csv(file = biogrid_all_path, header = TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606) %>%
    filter(Experimental.System.Type == "physical")
  
  mv_biogrid = read.csv(file = biogrid_mv_path, header = TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606)
  
  ## Need to mirror the biogrid files
  biogrid_mirror = all_biogrid
  biogrid_mirror$Entrez.Gene.Interactor.A = all_biogrid$Entrez.Gene.Interactor.B
  biogrid_mirror$Entrez.Gene.Interactor.B = all_biogrid$Entrez.Gene.Interactor.A
  all_biogrid = rbind(all_biogrid,biogrid_mirror)
  
  mv_mirror = mv_biogrid
  mv_mirror$Entrez.Gene.Interactor.A = mv_biogrid$Entrez.Gene.Interactor.B
  mv_mirror$Entrez.Gene.Interactor.B = mv_biogrid$Entrez.Gene.Interactor.A
  mv_biogrid = rbind(mv_biogrid,mv_mirror)
  # Make this column to merge in
  mv_biogrid$multivalidated = TRUE
  # Simplify the mv_biogrid data before merge
  mv_biogrid = mv_biogrid %>%
    group_by(Entrez.Gene.Interactor.A,Entrez.Gene.Interactor.B) %>%
    summarize(Multivalidated = any(multivalidated),.groups = "drop")
  
  # Merge
  all_biogrid = left_join(all_biogrid,mv_biogrid,by = c("Entrez.Gene.Interactor.A","Entrez.Gene.Interactor.B"))
  
  # Group and summarize
  summ_biogrid = all_biogrid %>%
    group_by(Entrez.Gene.Interactor.A,Entrez.Gene.Interactor.B) %>%
    summarize(in.BioGRID = (length(Entrez.Gene.Interactor.A)>=1),
              Evidence.Weight = length(Entrez.Gene.Interactor.A),
              Experimental.Systems = paste(Experimental.System,collapse = ';'),
              Authors = paste(Author,collapse = ';'),
              Publications = paste(Publication.Source, collapse = ';'),
              Entrez.Gene.Interactor.A = Entrez.Gene.Interactor.A[[1]],
              Entrez.Gene.Interactor.B = Entrez.Gene.Interactor.B[[1]],
              Multivalidated = any(Multivalidated),
              .groups="drop"
    )
  
  # Write the summary biogrid file
  write_csv(summ_biogrid, paste(output_dir,'/biogrid_summary.csv',sep=''))
}














