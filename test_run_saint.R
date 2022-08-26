

# library("devtools")
# devtools::install_github("smarasolo/cRomppass")


library(reticulate)
library(cRomppass)
library(tidyverse)
library(DarkKinaseTools)
library(org.Hs.eg.db)
library(dplyr)
source("C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/compute_wds.R")
#source("C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/Cytoscape.R")

# I have absolutely no idea what this is supposed to do...
select <- get(x="select", pos = "package:dplyr")

# Set paths
SAINT_path = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/build/SAINTexpress-spc.exe'
ED_path = 'C:/Users/plutzer/Work/IDG_pipeline/ED_DB.csv'
PG_path = 'C:/Users/plutzer/Work/IDG_pipeline/proteinGroups.txt'
output_dir = 'C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_blank'
uniprot_map_path = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/uniprot_mapping.tsv.zip'
biogrid_mv_path = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/BIOGRID-MV-Physical-4.4.211.tab3.txt'
biogrid_all_path = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/BIOGRID-ALL-4.4.211.tab3.txt'


resampling_iterations = 100

setwd(output_dir)

# For some unknown reason, this line is necessary even though it does nothing.
py_run_string("print(\"Python is Running.\")")

# Generate SAINT inputs from old scoreAPMS python script
system(paste("python C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/score_APMS_noSAINT.py","--experimentalDesign",ED_path,"--proteinGroups",PG_path,"--outputPath",output_dir))

interaction_path = paste(output_dir,'/interaction.txt',sep = '')
prey_path = paste(output_dir,'/prey.txt',sep = '')
bait_path = paste(output_dir,'/bait.txt',sep = '')

# Run SAINT
system(paste(SAINT_path,interaction_path,prey_path,bait_path))

# Smaranda's CompPASS script
filename = paste(output_dir,"/to_CompPASS.csv",sep='')
to_comp = read.csv(file = filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

comp_out = comppass(to_comp, stats = NULL, norm.factor = 0.98)
write.table(comp_out, file=paste(output_dir,'/compPASS.csv',sep=''), sep = "\t")

###  Comppass FDR goes here
print("Starting CompPASS resampling...")
pb = txtProgressBar(min = 0, max = resampling_iterations, initial = 0) 
comp_resample = resample_AvePSM(to_comp,comp_out,n.experiments = length((to_comp %>% count(Experiment.ID))[[1]]),suffix = "test",test=T)
for (i in 1:resampling_iterations) { # Adjust number of iterations as needed
  comp_resample = resample_AvePSM(to_comp,comp_resample,n.experiments = length((to_comp %>% count(Experiment.ID))[[1]]),suffix = as.character(i),test=F)
  setTxtProgressBar(pb,i)
}
close(pb)
write.table(comp_resample,file=paste(output_dir,'/resample_comp_out.csv',sep=''),sep="\t")

wds = comp_resample %>% select(contains("WD_") & !contains("test"))

comp_out$perm_fdrs = apply(as.array(1:length(comp_resample$WD)), MARGIN = 1, FUN = function(i) {sum(wds[i,]>comp_resample$WD[[i]])/resampling_iterations})


## 

# Smaranda's Merge script
saint <- read.csv(file=paste(output_dir,'/list.txt',sep=''), sep="\t")

merge = left_join(as.data.frame.matrix(comp_out), saint, by = c("Experiment.ID" = "Bait", "Prey"))

write_csv(merge, paste(output_dir,'/Merge_CompPASS_SAINT.csv',sep=''))


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
data <- rename(data, c("Gene_Name.x" = "Prey.Gene.Name", "Gene_Synonym.x" = "Prey.Gene.Synonym", "GeneID.x" = "Prey.GeneID", "First_GeneID.x" = "First.Prey.GeneID", 
                       "Gene_Name.y" = "Bait.Gene.Name", "Gene_Synonym.y" = "Bait.Gene.Synonym", "GeneID.y" = "Bait.GeneID", "First_GeneID.y" = "First.Bait.GeneID"))
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

################################################################################

# Newest biogrid solution
all_biogrid = read.csv(file = biogrid_all_path, header = TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606) %>%
  filter(Experimental.System.Type == "physical")

mv_biogrid = read.csv(file = biogrid_all_path, header = TRUE, sep="\t", stringsAsFactors = FALSE) %>%
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

data = left_join(data,summ_biogrid, by=c('Bait.GeneID'='Entrez.Gene.Interactor.A','Prey.GeneID'='Entrez.Gene.Interactor.B'))

################################################################################
# Filter out control baits using ED file
test_baits = read.csv(ED_path, header = TRUE, stringsAsFactors = FALSE) %>%
  filter(Type == 'T')

#Filter 
data.filter <- filter(data, BFDR <= 0.05, AvgP >= 0.7) %>%
  filter(Experiment.ID %in% test_baits$Bait)

data.testbaits = data %>% filter(Experiment.ID %in% test_baits$Bait)

# Sort by comppass scores
data.filter.comp <- arrange(data.filter, desc(WD))
data.testbaits.comp = arrange(data.testbaits, desc(WD))

# Write these output files
write.csv(data.filter.comp, paste(output_dir,"/Annotated_Merge_Saint_filter.csv",sep='',na="",row.names=F))
write.csv(data.testbaits.comp, paste(output_dir,"/Annotated_Merge_Saint.csv",sep='',na="",row.names=F))

# Prey-Prey
#Create separate data tables for prey-prey interactions

baits <- unique((data.filter.comp %>% filter(!is.na(BaitGene), is_Bait == TRUE, BaitGene == PreyGene))$Bait)

geneList = as.vector(data$`First.Prey.GeneID`[!is.na(data$`First.Prey.GeneID`)])

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

nodes_combined = select(data.testbaits,Prey.Gene.Name,class,is_Bait,WD,perm_fdrs,FoldChange,GO.Slim,Prey.Gene.Synonym,First.Bait.GeneID,First.Prey.GeneID,Bait.Gene.Name)
edges_combined = select(data.testbaits,Experiment.ID,Prey.Gene.Name,in.BioGRID,Multivalidated,Evidence.Weight,Experimental.Systems,Authors,Publications,WD,perm_fdrs,Bait.Gene.Name)
nodes_combined_filtered = select(data.filter,Prey.Gene.Name,class,is_Bait,WD,perm_fdrs,FoldChange,GO.Slim,Prey.Gene.Synonym,First.Bait.GeneID,First.Prey.GeneID,Bait.Gene.Name)
edges_combined_filtered = select(data.filter,Experiment.ID,Prey.Gene.Name,in.BioGRID,Multivalidated,Evidence.Weight,Experimental.Systems,Authors,Publications,WD,perm_fdrs,Bait.Gene.Name)

# Commented out the combined edges... that shit was too big

# Add the interaction column
edges_combined$Interaction = "Strep APMS"
edges_combined_filtered$Interaction = "Strep APMS"

dir.create(paste(output_dir,"/Cytoscape_Outputs",sep=''))

write_csv(nodes_combined, paste(output_dir,"/Cytoscape_Outputs/Combined_nodes.csv",sep=''))
write_csv(nodes_combined_filtered, paste(output_dir,"/Cytoscape_Outputs/Combined_nodes_filtered.csv",sep=''))

for (mybait in baits) {
  
  #mybait <- "P24941"
  
  bait.data.filter <- data.testbaits.comp %>% filter(Bait == mybait)
  num.interactors <- min(max(10, nrow(bait.data.filter)*0.05), nrow(data.filter.comp))
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
    
    dir.create(paste(output_dir,"/Cytoscape_Outputs/",bait,sep=''),na="",row.names=F)
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
    
    if (length(edges_bait)>=1){
      prey.prey.select$Experiment.ID = edges_bait$Experiment.ID[[1]] # These shouldn't change so just first is fine
    }
    else {
      prey.prey.select$Experiment.ID = '-'
    }
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


