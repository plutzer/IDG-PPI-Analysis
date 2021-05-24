# Title     : annotate_filter
# Objective : Add GO annotation to the merged SAINT and CompPASS file
# Created by: Smaranda
# Created on: 9/3/2020

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")

library(tidyverse)
library(org.Hs.eg.db)
library(DarkKinaseTools)
source("Cytoscape.R")

select <- get(x="select", pos = "package:dplyr")

uniprot.mapping <- read_tsv("annotations/uniprot_mapping.tsv.zip")

st <- read.csv(file='output/Merge_CompPASS_SAINT.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Bait Uniprot
#Grab Unique Baits
unique.Baits <- unique(st$Bait)

################################################################################
#Entrez GeneID

#Prey
# extract Prey Uniprot Identifier
data <- separate(data = st, col = `Prey`, into = c("First.Prey.Uniprot"), sep=";", remove=F, extra="drop")

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
data <- data[,c(1:7, 32:39, 8:31)]
data <- rename(data, c(Gene_Name.x = "Prey.Gene.Name", Gene_Synonym.x = "Prey.Gene.Synonym", GeneID.x = "Prey.GeneID", First_GeneID.x = "First.Prey.GeneID", 
                               Gene_Name.y = "Bait.Gene.Name", Gene_Synonym.y = "Bait.Gene.Synonym", GeneID.y = "Bait.GeneID", First_GeneID.y = "First.Bait.GeneID"))

################################################################################
#Annotations

xx <- as.list(org.Hs.egGO2ALLEGS)
x <-c("GO:0005634","GO:0005737","GO:0005856","GO:0005768","GO:0005783","GO:0005576","GO:0005794",
      "GO:0005764","GO:0005739","GO:0005777","GO:0005886","GO:0031982","GO:0005911")
go <- c("Nucleus","Cytoplasm","Cytoskeleton","Endosome","ER","Extracellular","Golgi","Lysosome",
        "Mitochondria","Peroxisome","Plasma_Membrane","Vesicles","Cell_Junction")


for(i in 1:length(go)) {
  goids <- xx[x[i]]
  goslim <- tibble(Prey.GeneID = goids[[1]], Evidence = names(goids[[1]]))
  goslim <- filter(goslim, "IEA" != Evidence)
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
data$GO.Slim <- apply(data, 1, function(x) {
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
biogrid <- read.csv(file = 'C:/Users/smaranda/Documents/SmarandaSolomon/BIOGRID-MV-Physical-4.2.193.tab3.txt', header = TRUE, sep="\t", stringsAsFactors = FALSE) %>%
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

data$in_BioGRID <- apply(data, 1, function(x) {
  nrow(filter(interactions, (`interactor_min` == as.numeric(x[["Prey.GeneID"]]) & `interactor_max` == as.numeric(x[["Bait.GeneID"]])) |
                (`interactor_max` == as.numeric(x[["Prey.GeneID"]]) & `interactor_min` == as.numeric(x[["Bait.GeneID"]]))
  )) > 0
}) 

################################################################################
#Filter and Write files

#Write file with no filter
write_csv(data, 'output/Annotated_Merge_NO_FILTER.csv')

#Filter SAINT
data.filter <- filter(data, BFDR <= 0.05, AvgP >= 0.7)

#Filter CompPASS
data.filter <- arrange(data.filter, desc(WD))

#Write file with all experiments
write_csv(data.filter, 'output/Annotated_Merge_Saint_filter.csv')

#Filter by top 5% or top min *KEEPING Everything*
baits <- unique((data.filter %>% filter(!is.na(BaitGene), is_Bait == TRUE, BaitGene == PreyGene))$Bait)
all.data.filter <- data.filter[0,]

for (mybait in baits) {
  bait.data.filter <- data.filter %>% filter(Bait == mybait)
  num.interactors <- min(max(10, nrow(bait.data.filter)*0.05), nrow(bait.data.filter))
  all.data.filter <- all.data.filter %>% add_row(bait.data.filter[1:num.interactors, ])
  
  prey.prey.inter <- filter(interactions, (`interactor_min` %in%  bait.data.filter$First.Prey.GeneID), (`interactor_max` %in%  bait.data.filter$First.Prey.GeneID))
  write_csv(prey.prey.inter, str_c('output/Prey_Prey_Interactions/', mybait, '.csv'))
}

#Write file with filtered experiments, keeping everything
write_csv(all.data.filter, 'output/Annotated_Merge_All_filtered.csv')

################################################################################
#Filter by top 5% or top min *For DKK - no controls *
baits.dkk <- unique((data.filter %>% filter(!is.na(BaitGene), is_Bait == TRUE, BaitGene == PreyGene, !is.na(class_2019)))$Bait)
all.data.dkk <- data.filter[0,]

for (mybait in baits.dkk) {
  bait.data.filter <- data.filter %>% filter(Bait == mybait)
  num.interactors <- min(max(10, nrow(bait.data.filter)*0.05), nrow(bait.data.filter))
  all.data.dkk <- all.data.dkk %>% add_row(bait.data.filter[1:num.interactors, ])
}

all.data.dkk <- all.data.dkk %>%
  dplyr::filter(!grepl('CSNK1G1|CSNK1G2|CSNK1G3', Experiment.ID))

################################################################################
#Write file with filtered experiments, only Dark kinases
write_csv(all.data.dkk, 'output/Annotated_Merge_filtered_DKK.csv')

cyto <- to.Cytoscape(all.data.filter, data)

#Write file with filtered experiments, only Dark kinases
write_csv(cyto, 'output/Cytoscape.csv')
