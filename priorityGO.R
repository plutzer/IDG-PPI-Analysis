# Title     : priorityGO
# Objective : Make a priority Go column based on if matched to bait, if not then the first
# Created by: Smaranda Solomon
# Created on: 4/13/22

library(tidyverse)
library(org.Hs.eg.db)

select <- get(x="select", pos = "package:dplyr")
#rename <- get(x = "rename", pos = "package:dplyr")

#priority.GO <- function(all.data.filter, data){
  
st <- read.csv(file='output/bait_cyto.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

bait.GO <- st %>%
  select(Bait, Bait.Gene.Name, Bait.GeneID)

xx <- as.list(org.Hs.egGO2ALLEGS)
x <-c("GO:0005634","GO:0005737","GO:0005856","GO:0005768","GO:0005783","GO:0005576","GO:0005794",
      "GO:0005764","GO:0005739","GO:0005777","GO:0005886","GO:0031982","GO:0005911")
go <- c("Nucleus","Cytoplasm","Cytoskeleton","Endosome","ER","Extracellular","Golgi","Lysosome",
        "Mitochondria","Peroxisome","Plasma_Membrane","Vesicles","Cell_Junction")


for(i in 1:length(go)) {
  goids <- xx[x[i]]
  goslim <- tibble(Bait.GeneID = goids[[1]], Evidence = names(goids[[1]]))
  goslim <- (goslim %>% 
               group_by(Bait.GeneID)%>% 
               summarise(n=n()))
  bait.GO <- left_join(bait.GO, goslim, by=c("Bait.GeneID" = "Bait.GeneID")) %>% 
    mutate(
      goNames = !is.na(n) 
    ) %>% 
    select(-`n`) %>% 
    rename(goNames = go[i])
}


#Merge GO Slim in Single Column
bait.GO$GO.Slim <- apply(bait.GO, 1, function(x) {
  str_c(go[as.logical(c(x[["Nucleus"]],x[["Cytoplasm"]],x[["Cytoskeleton"]], x[["Endosome"]], x[["ER"]],
                        x[["Extracellular"]], x[["Golgi"]], x[["Lysosome"]], x[["Mitochondria"]],
                        x[["Peroxisome"]], x[["Plasma_Membrane"]], x[["Vesicles"]], x[["Cell_Junction"]]))], 
        collapse = ";")
}
)


#}

