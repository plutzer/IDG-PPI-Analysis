

# library("devtools")
# devtools::install_github("smarasolo/cRomppass")


library(reticulate)
library(cRomppass)
library(tidyverse)
library(DarkKinaseTools)
library(org.Hs.eg.db)
source("C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/Cytoscape.R")



# Set paths
SAINT_path = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/build/SAINTexpress-spc.exe'
ED_path = 'C:/Users/plutzer/Work/IDG_pipeline/ED_DB.csv'
PG_path = 'C:/Users/plutzer/Work/IDG_pipeline/proteinGroups.txt'
output_dir = 'C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_blank'
uniprot_map_path = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/uniprot_mapping.tsv.zip'
biogrid_path = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/BIOGRID-MV-Physical-4.4.211.tab3.txt'

setwd(output_dir)

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

# Smaranda's Merge script
saint <- read.csv(file=paste(output_dir,'/list.txt',sep=''), sep="\t")

merge = left_join(as.data.frame.matrix(comp_out), saint, by = c("Experiment.ID" = "Bait", "Prey"))

write_csv(merge, paste(output_dir,'/Merge_CompPASS_SAINT.csv',sep=''))

# Smaranda's Annotate Script

# Create Prey-Prey Directory
if(!dir.exists(paste(output_dir, '/Prey_Prey_Interactions/',sep=''))){
  dir.create(paste(output_dir, '/Prey_Prey_Interactions/',sep=''))
} else {
  FALSE
}

# I have absolutely no idea what this is supposed to do...
select <- get(x="select", pos = "package:dplyr")


# # This needs to be updated because I don't have this file...
uniprot.mapping <- read_tsv(uniprot_map_path)

# List of unique baits
unique.Baits <- unique(merge$Bait)

################################################################################
#Entrez GeneID

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
data <- data[c(1:7, 32:38, 8:31)]
data <- rename(data, c("Gene_Name.x" = "Prey.Gene.Name", "Gene_Synonym.x" = "Prey.Gene.Synonym", "GeneID.x" = "Prey.GeneID", "First_GeneID.x" = "First.Prey.GeneID", 
                       "Gene_Name.y" = "Bait.Gene.Name", "Gene_Synonym.y" = "Bait.Gene.Synonym", "GeneID.y" = "Bait.GeneID", "First_GeneID.y" = "First.Bait.GeneID"))
data <- data[,c(1,12,13,2:4,14,15,8,9,5:7,10,11,16:38)]

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

#Priority GO Column: Either matched with bait or from uniprot
#data$GO.Priority <- apply(data, 1, function(x)){
#}

#bait.data.filter <- data.filter %>% filter(Bait == mybait)
#num.interactors <- min(max(10, nrow(bait.data.filter)*0.05), nrow(bait.data.filter))
#all.data.filter <- all.data.filter %>% add_row(bait.data.filter[1:num.interactors, ])


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
#BioGrid Annotations

#Biogrid read  in
biogrid <- read.csv(file = biogrid_path, header = TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606)

#Is In BioGrid (T/F)
geneList = as.vector(data$`First.Prey.GeneID`[!is.na(data$`First.Prey.GeneID`)])

interactions <- biogrid %>%
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

#CHANGE AS.NUMERIC AND PUT IT ABOVE SO IT GOES FASTER - this doesn't work
# data$in_BioGRID <- apply(data, 1, function(x) {
#   nrow(filter(interactions, (`interactor_min` == as.numeric(x[["Prey.GeneID"]]) & `interactor_max` == as.numeric(x[["Bait.GeneID"]])) |
#                 (`interactor_max` == as.numeric(x[["Prey.GeneID"]]) & `interactor_min` == as.numeric(x[["Bait.GeneID"]]))
#   )) > 0
# })

# This solution works and is much faster
bait_prey_pairs = transpose(as.list(data[, c("Bait.GeneID","Prey.GeneID")]))
data$in_BioGRID = lapply(bait_prey_pairs, FUN=function(x){(any(interactions$interactor_min==as.numeric(x$Bait.GeneID) & interactions$interactor_max==as.numeric(x$Prey.GeneID)) |
                                                                  any(interactions$interactor_max==as.numeric(x$Bait.GeneID) & interactions$interactor_min==as.numeric(x$Prey.GeneID)))})
################################################################################
#Filter SAINT

#Write file with no filter
write_csv(data, paste(output_dir,'/Annotated_Merge_NO_FILTER.csv',sep=''))

#Filter SAINT
data.filter <- filter(data, BFDR <= 0.05, AvgP >= 0.7)

write_csv(data.filter, paste(output_dir,'/Annotated_Merge_Saint_filter.csv',sep=''))

################################################################################

#Filter CompPASS
#Arrange first
data.filter.comp <- arrange(data.filter, desc(WD))

#Filter by top 5% or top min *KEEPING Everything*
baits <- unique((data.filter.comp %>% filter(!is.na(BaitGene), is_Bait == TRUE, BaitGene == PreyGene))$Bait)
all.data.filter <- data.filter.comp[0,]

#Create separate data tables for prey-prey interactions
for (mybait in baits) {
  
  #mybait <- "P24941"
  
  bait.data.filter <- data.filter.comp %>% filter(Bait == mybait)
  num.interactors <- min(max(10, nrow(bait.data.filter)*0.05), nrow(data.filter.comp))
  bait.data.filter.comp <- bait.data.filter[1:num.interactors, ]
  all.data.filter <- all.data.filter %>% add_row(bait.data.filter.comp)
  
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
    select(Prey.1.Gene.Name, Prey.1.Uniprot, Prey.1.Entrez.ID, WD.1, Prey.Gene.Name, Canonical.First.Prey.Uniprot, Prey.2.Entrez.ID, WD) %>%
    rename(c(Canonical.First.Prey.Uniprot = "Prey.2.Uniprot", Prey.Gene.Name = "Prey.2.Gene.Name", WD = "WD.2")) %>%
    unique()
  
  #Add validation column (from Biogrid in this case) if dataframe has data in it
  if(nrow(prey.prey.final) != 0){
    prey.prey.final$Source <- 'Biogrid'
  }
  
  #write individual csv files for each bait
  write_csv(prey.prey.final, str_c(paste(output_dir,'/Prey_Prey_Interactions/',sep=''), paste(unique(bait.data.filter$Bait.Gene.Name),"_",mybait,'.csv', sep = "")))
}

#Write file with filtered experiments, keeping everything
write_csv(all.data.filter, paste(output_dir,'/Annotated_Merge_All_filtered.csv',sep=''))

################################################################################
#Filter by top 5% or top min *For DKK - no controls *
baits.dkk <- unique((data.filter.comp %>% filter(!is.na(BaitGene), is_Bait == TRUE, BaitGene == PreyGene, !is.na(class_2019)))$Bait)
all.data.dkk <- data.filter.comp[0,]

for (mybait in baits.dkk) {
  bait.data.filter <- data.filter.comp %>% filter(Bait == mybait)
  num.interactors <- min(max(10, nrow(bait.data.filter)*0.05), nrow(bait.data.filter))
  all.data.dkk <- all.data.dkk %>% add_row(bait.data.filter[1:num.interactors, ])
}

all.data.dkk <- all.data.dkk %>%
  dplyr::filter(!grepl('CSNK1G1|CSNK1G2|CSNK1G3', Experiment.ID))

################################################################################
#Write file with filtered experiments, only Dark kinases
write_csv(all.data.dkk, paste(output_dir,'/Annotated_Merge_filtered_DKK.csv',sep=''))

cyto <- to.Cytoscape(all.data.filter, data)




