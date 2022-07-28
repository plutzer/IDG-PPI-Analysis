
library(dplyr)

filename = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/BIOGRID-ALL-4.4.211.tab3.txt'


function process_bioGRID(filename) {
  # Open the BioGRID file
  biogrid = read.csv(file = filename, header = TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606)

  # Get a list of every unique bait-prey pair
  pairs = paste(biogrid$Entrez.Gene.Interactor.A,biogrid$Entrez.Gene.Interactor.B,sep = '-')
  
  
}
