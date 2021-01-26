# Title     : annotate
# Objective : Add GO annotation to the merged SAINT and CompPASS file
# Created by: Smaranda
# Created on: 9/3/2020

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")

library(tidyverse)
library(org.Hs.eg.db)
library(DarkKinaseTools)

select <- get(x="select", pos = "package:dplyr")

uniprot.mapping <- read_tsv("annotations/uniprot_mapping.tsv.zip")

st <- read.csv(file='output/Merge_CompPASS_SAINT.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Bait Uniprot
#Grab Unique Baits
unique.Baits <- unique(st$Bait)

# extract Uniprot Identifier
proteins <- separate(data = st, col = `Prey`, into = c("Uniprot.Identifier"), sep=";", remove=F, extra="drop")

# Remove isoform from the Uniprot Identifier
proteins <- separate(proteins, "Uniprot.Identifier", c("Unique.Uniprot.Identifier"), sep="-", remove=F, extra="drop")

# fill Entrez GeneID
proteins <- left_join(proteins, uniprot.mapping, by=c("Unique.Uniprot.Identifier" = "UniProt"))

# renaming columns
#proteins <- rename(proteins, c("Uniprot.Identifier" = "First.Prey.Uniprot", "Unique.Uniprot.Identifier" = "Canon.First.Prey.Uniprot"))

################################################################################

#GO ID's

xx <- as.list(org.Hs.egGO2ALLEGS)
x <-c("GO:0005634","GO:0005737","GO:0005856","GO:0005768","GO:0005783","GO:0005576","GO:0005794",
      "GO:0005764","GO:0005739","GO:0005777","GO:0005886","GO:0031982","GO:0005911")
go <- c("Nucleus","Cytoplasm","Cytoskeleton","Endosome","ER","Extracellular","Golgi","Lysosome",
        "Mitochondria","Peroxisome","Plasma_Membrane","Vesicles","Cell_Junction")


for(i in 1:length(go)) {
  goids <- xx[x[i]]
  goslim <- tibble(Gene.ID = goids[[1]], Evidence = names(goids[[1]]))
  goslim <- filter(goslim, "IEA" != Evidence)
  goslim <- (goslim %>% 
               group_by(Gene.ID)%>% 
               summarise(n=n()))
  proteins <- left_join(proteins, goslim, by=c("GeneID" = "Gene.ID")) %>% 
    mutate(
      goNames = !is.na(n) 
    ) %>% 
    select(-`n`) %>% 
    rename(goNames = go[i])
}

#Merge GO Slim in Single Column
proteins$GO.Slim <- apply(proteins, 1, function(x) {
  str_c(go[as.logical(c(x[["Nucleus"]],x[["Cytoplasm"]],x[["Cytoskeleton"]], x[["Endosome"]], x[["ER"]],
                        x[["Extracellular"]], x[["Golgi"]], x[["Lysosome"]], x[["Mitochondria"]],
                        x[["Peroxisome"]], x[["Plasma_Membrane"]], x[["Vesicles"]], x[["Cell_Junction"]]))], 
        collapse = ";")
  }
)

#Filter SAINT
proteins <- filter(proteins, BFDR <= 0.05, AvgP >= 0.7)

#Filter CompPASS
proteins <- arrange(proteins, desc(WD))

#Annotate Dark Kinases
proteins <- left_join(proteins, all_kinases, by=c("First_GeneID" = "entrez_id"))

#Is it a Bait Column
proteins$is_Bait <- apply(proteins, 1, function(x) {
  any(str_detect(unique.Baits, x[["Unique.Uniprot.Identifier"]]))
  }
)

#Nice bait name column addition by grabbing from prey name if it pairs with is_bait column
baitTable <- proteins %>%
  filter(is_Bait == TRUE) %>%
  filter(str_detect(Bait, Prey)) %>%
  select(Bait, BaitGene = PreyGene) %>%
  distinct() 

proteins <- left_join(proteins, baitTable, by="Bait")

#Write file with all experiments
write_csv(proteins, 'output/Annotated_Merge_ALL.csv')

#Filter by top 5% or top min *KEEPING Everything*

baits <- unique((proteins %>% filter(!is.na(BaitGene), is_Bait == TRUE, BaitGene == PreyGene))$Bait)
all.data <- proteins[0,]

for (mybait in baits) {
  
  bait.data <- proteins %>% filter(Bait == mybait)
  num.interactors <- min(max(10, nrow(bait.data)*0.05), nrow(bait.data))
  all.data <- all.data %>% add_row(bait.data[1:num.interactors, ])
}

#Write file with filtered experiments, keeping everything
write_csv(all.data, 'output/Annotated_Merge_filtered.csv')

#Filter by top 5% or top min *For DKK - no controls *

baits <- unique((proteins %>% filter(!is.na(BaitGene), is_Bait == TRUE, BaitGene == PreyGene, !is.na(class_2019)))$Bait)
all.data <- proteins[0,]

for (mybait in baits) {
  bait.data <- proteins %>% filter(Bait == mybait)
  num.interactors <- min(max(10, nrow(bait.data)*0.05), nrow(bait.data))
  all.data <- all.data %>% add_row(bait.data[1:num.interactors, ])
}

all.data <- all.data %>%
  dplyr::filter(!grepl('CSNK1G1|CSNK1G2|CSNK1G3', Experiment.ID))

#Write file with filtered experiments, only Dark kinases
write_csv(all.data, 'output/Annotated_Merge_filtered_DKK.csv')

