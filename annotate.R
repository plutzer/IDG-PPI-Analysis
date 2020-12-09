# Title     : annotate
# Objective : Add GO annotation to the merged SAINT and CompPASS file
# Created by: Smaranda
# Created on: 9/3/2020

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
baitName <- tibble(proteins$Unique.Uniprot.Identifier, Nice.Bait.Name = NA, proteins$PreyGene, proteins$is_Bait) %>%
  mutate(Nice.Bait.Name = ifelse(proteins$is_Bait == TRUE, proteins$Unique.Uniprot.Identifier, Nice.Bait.Name)) %>%
  dplyr::rename(Unique.Uniprot.Identifier = "proteins$Unique.Uniprot.Identifier", Bait.Gene = "proteins$PreyGene" , is_Bait = "proteins$is_Bait")

baitName <- filter(baitName, is_Bait == TRUE) %>%
  select(-is_Bait)

proteins <- left_join(proteins, baitName, by=c("Unique.Uniprot.Identifier" = "Unique.Uniprot.Identifier"))

#Write file with all experiments
write_csv(proteins, 'output/Annotated_Merge.csv')

#Filter Top Percent of CompPASS 
top <- 0.05*nrow(proteins)

proteins <- proteins[1:top,]

#Write file with only top percent of CompPass
write_csv(proteins, 'output/Annotated_Merge_top5.csv')
