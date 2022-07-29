
library(dplyr)

filename = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/BIOGRID-ALL-4.4.211.tab3.txt'
dataframe = data

# Open the BioGRID file
biogrid_toy = read.csv(file = filename, header = TRUE, sep="\t", stringsAsFactors = FALSE) %>% 
  filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606)

# Is the biogrid dataset symmetrical?
# for (bait in unique(biogrid_toy$Entrez.Gene.Interactor.A)) {
#   if (length(filter(biogrid_toy,biogrid_toy$Entrez.Gene.Interactor.A == bait)[[1]]) == length(filter(biogrid_toy,biogrid_toy$Entrez.Gene.Interactor.B == bait)[[1]])) {
#     foo='bar'
#   }
#   else {
#     print(F)
#   }
# }
# NO IT IS NOT

# I'm just gonna make it symmetrical
biogrid_mirror = biogrid_toy
biogrid_mirror$Entrez.Gene.Interactor.A = biogrid_toy$Entrez.Gene.Interactor.B
biogrid_mirror$Entrez.Gene.Interactor.B = biogrid_toy$Entrez.Gene.Interactor.A
biogrid_test = rbind(biogrid_toy,biogrid_mirror)

# Summarize the biogrid data
test = biogrid_test %>%
  group_by(Entrez.Gene.Interactor.A,Entrez.Gene.Interactor.B) %>%
  summarize(all = (length(Entrez.Gene.Interactor.A)>=1),
            biogrid_evidence_weight = length(Entrez.Gene.Interactor.A),
            exp_systems = paste(Experimental.System,collapse = ';'),
            authors = paste(Author,collapse = ';'),
            publications = paste(Publication.Source, collapse = ';'),
            Entrez.Gene.Interactor.A = Entrez.Gene.Interactor.A[[1]],
            Entrez.Gene.Interactor.B = Entrez.Gene.Interactor.B[[1]],
            .groups="drop"
  )

# Merge the dataframes on the IDs 
dataframe = left_join(dataframe,test, by=c('Bait.GeneID'='Entrez.Gene.Interactor.A','Prey.GeneID'='Entrez.Gene.Interactor.B'))




# Get a list of every unique bait-prey pair
pairs = paste(biogrid$Entrez.Gene.Interactor.A,biogrid$Entrez.Gene.Interactor.B,sep = '-')

return(dataframe)
