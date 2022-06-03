# Title     : csnk1g_figures_NEW
# Objective : Output the file for cytoscape and a heatmap for CSNK1g1 paper
# Created by: Smaranda Solomon
# Created on: 10/29/2021

library(devtools)
library(tidyverse)
library(circlize)

organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

select <- get(x="select", pos = "package:dplyr")
rename <- get(x="rename", pos = "package:dplyr")

################################################################################ 

sec <- read.csv(file='output/Annotated_Merge_All_filtered.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::filter(grepl('CSNK1G1|CSNK1G2|CSNK1G3', Bait.Gene.Name)) %>%
  select(Experiment.ID, Bait.Gene.Name, Prey.Gene.Name, BFDR, AvgP, WD, Spec,Nucleus,Cytoplasm,Cytoskeleton,Endosome,ER,Extracellular,Golgi,Lysosome,
         Mitochondria,Peroxisome,Plasma_Membrane,Vesicles,Cell_Junction)

st <- read.csv(file='output/Annotated_Merge_Saint_filter.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

just.saint <- st %>% 
  dplyr::filter(grepl('AAK1|AXIN1|AXIN2|APC|LRP5|LRP6|DVL1|DVL2|DVL3|CTNNB1|FZD1|FZD2|FZD3|FZD4|FZD5|FZD6|FZD7|FZD8|
                      GSK3A|GSK3B|AMER1|CSNK1A|CSNK1D|CSNK1E|CSNK1G1|CSNK1G2|CSNK1G3|ROR1|ROR2|RYK|ANKRD6|CELSR1|CELSR2|
                      FAT|VANGL|NOTCH', Prey.Gene.Name)) %>%
  select(Experiment.ID, Bait.Gene.Name, Prey.Gene.Name, BFDR, AvgP, WD, Spec,Nucleus,Cytoplasm,Cytoskeleton,Endosome,ER,Extracellular,Golgi,Lysosome,
         Mitochondria,Peroxisome,Plasma_Membrane,Vesicles,Cell_Junction, Cell_Line, hits, Neg, Pos.high, Pos.low, Att)

saint.comp <- full_join(sec, just.saint) %>%
  dplyr::filter(grepl('CSNK1G1|CSNK1G2|CSNK1G3', Bait.Gene.Name)) %>%
  unique()

saint.comp <- saint.comp %>%
  unite("hits_Neg_Pos.high_Pos.low_Att", hits:Att) 

################################################################################

screen <- read.csv(file='output/Annotated_Merge_Saint_filter.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  select(Experiment.ID, Bait.Gene.Name, Prey.Gene.Name, Spec, hits, Neg, Pos.high, Pos.low, Att, Biechele.Fold, AGGF.Hit) %>%
  filter(grepl('CSNK1G1|CSNK1G2|CSNK1G3', Bait.Gene.Name)) %>%
  mutate(screen=TRUE)

screen <- screen %>%
  filter(!is.na(hits) | !is.na(Neg) | !is.na(Pos.high) | !is.na(Pos.low) | !is.na(Att) | !is.na(Biechele.Fold) | AGGF.Hit == TRUE) %>%
  filter(!is.na(hits) | Neg == TRUE | Pos.high == TRUE | Pos.low == TRUE | Att == TRUE | Biechele.Fold == TRUE | AGGF.Hit == TRUE) %>%
  unique()

extra.prey <- read.csv(file='output/Annotated_Merge_NO_FILTER.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::filter(grepl('CSNK1G1|CSNK1G2|CSNK1G3', Bait.Gene.Name)) %>%
  dplyr::filter(Prey.Gene.Name %in% unique(saint.comp$Prey.Gene.Name)) %>%
  full_join(., screen)%>%
  select(Bait.Gene.Name, Prey.Gene.Name, Spec, screen) %>%
  unique()

extra.prey <- extra.prey %>%
  arrange(., Prey.Gene.Name, screen) %>%
  group_by(Prey.Gene.Name) %>%
  fill(screen)

write_csv(extra.prey, 'output/extra.prey.csv')

################################################################################

saint.and.comp <- sec %>%
  select(Prey.Gene.Name) %>%
  mutate(source = "Saint.CompPASS") %>%
  unique()

wnt.and.saint <- just.saint %>%
  select(Prey.Gene.Name)%>%
  mutate(source = "Saint.WNT") %>%
  unique()

screens.and.saint <- screen %>%
  select(Prey.Gene.Name) %>%
  mutate(source = "Saint.Screens") %>%
  unique()

first <- full_join(saint.and.comp, wnt.and.saint, by = "Prey.Gene.Name") %>%
  select(Prey.Gene.Name, source.x)

first[is.na(first)] <- "Saint.WNT"

second <- full_join(first, screens.and.saint, by = "Prey.Gene.Name") %>%
  select(Prey.Gene.Name, source.x)

second[is.na(second)] <- "Saint.Screens"

levels <- full_join(extra.prey, second)

levels <- levels[complete.cases(levels[,1]),]

levels$source.x <- ifelse(levels$Prey.Gene.Name=="CSNK1G1", levels$source.x=="Bait", levels$source.x)

levels[levels == FALSE] <- "Bait"

levels$source.x <- ifelse(levels$Prey.Gene.Name=="CSNK1G2", levels$source.x=="Bait", levels$source.x)

levels[levels == FALSE] <- "Bait"

levels$source.x <- ifelse(levels$Prey.Gene.Name=="CSNK1G3", levels$source.x=="Bait", levels$source.x)

levels[levels == FALSE] <- "Bait"

levels <- levels %>%
  select(-Bait.Gene.Name, -Spec, -screen) %>%
  unique()

#################################################
#Heatmap made from LFQ intensities

lfq <- read_tsv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Experimental_Data/20211014_CK1G_All_Baits/proteinGroups.txt', guess_max = 10000) %>%
  select(`Gene names`, `LFQ intensity CSNK1G1-30mB-MA_1`, `LFQ intensity CSNK1G1-30mB-MA_2`, `LFQ intensity CSNK1G1-30mB-MA_3`, 
         `LFQ intensity CSNK1G2-30mB-MA_1`, `LFQ intensity CSNK1G2-30mB-MA_2`, `LFQ intensity CSNK1G2-30mB-MA_3`,
         `LFQ intensity CSNK1G3-30mB-MA_1`, `LFQ intensity CSNK1G3-30mB-MA_2`, `LFQ intensity CSNK1G3-30mB-MA_3`)

heatmap.lfq <- left_join(extra.prey, lfq, by = c("Prey.Gene.Name" = "Gene names")) %>%
  select(-Spec, -Bait.Gene.Name) %>%
  unique()

heatmap.lfq <- heatmap.lfq %>%
  filter(rowSums(across(matches("LFQ"))) > 0)

heatmap.lfq <- heatmap.lfq[,c(1, 3:11, 2)]

heatmap.levels.lfq <- left_join(heatmap.lfq, levels, by = "Prey.Gene.Name")

bait.order <- order(rowSums(heatmap.lfq[,-1]),decreasing=T)

heatmap.lfq <- heatmap.lfq[bait.order,]

heatmap.levels.lfq <- heatmap.levels.lfq[bait.order,]

#################################################
#Tracks

#Screens
track.names <- heatmap.lfq %>%
  select(screen)

track.names[is.na(track.names)] <- FALSE

track.names2 <- as.data.frame(track.names)

rownames(track.names2) <- track.names$Prey.Gene.Name 

track.names2$Prey.Gene.Name <- NULL

#WNT Related
wnt.and.saint$source <- TRUE

wnt.and.saint.track <- heatmap.lfq %>%
  left_join(wnt.and.saint, by = "Prey.Gene.Name") %>%
  select(source)

wnt.and.saint.track[is.na(wnt.and.saint.track)] <- FALSE

wnt.and.saint.track <- as.data.frame(wnt.and.saint.track)

rownames(wnt.and.saint.track) <- wnt.and.saint.track$Prey.Gene.Name 

wnt.and.saint.track$Prey.Gene.Name <- NULL

#################################################
#Matrix and heatmap creation

heatmap.lfq <- heatmap.lfq %>%
  select(-screen)

write_csv(heatmap.levels.lfq, 'output/csnk1g_heatmap_matrix_intensity.csv')

heatmap.mat.lfq <- as.matrix(heatmap.lfq) 

heatmap.mat.lfq <- heatmap.mat.lfq[,-1]

heatmap.mat.lfq <- matrix(as.numeric(unlist(heatmap.mat.lfq)), nrow = nrow(heatmap.mat.lfq))

rownames(heatmap.mat.lfq) <- heatmap.lfq$Prey.Gene.Name

heatmat.lfq.log <-  log2(heatmap.mat.lfq + 1)

colnames(heatmat.lfq.log) <- c("CSNK1G1_1", "CSNK1G1_2", "CSNK1G1_3", "CSNK1G2_1", "CSNK1G2_2", "CSNK1G2_3", "CSNK1G3_1", "CSNK1G3_2", "CSNK1G3_3")

#Color for screen track
colors <- list("screen" = c("TRUE" = "#C7706F", "FALSE" = "white"))

#Screen track command
track <- ComplexHeatmap::HeatmapAnnotation(df = track.names2, which = "row", col = colors, simple_anno_size = unit(.3, "cm"))

#Color for wnt track
colors.wnt <- list("source" = c("TRUE" = "#B0D656", "FALSE" = "white"))

#Wnt track command
wnt.track <- ComplexHeatmap::HeatmapAnnotation(df = wnt.and.saint.track, which = "row", col = colors.wnt, simple_anno_size = unit(.3, "cm"))

ComplexHeatmap::Heatmap(heatmat.lfq.log,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        column_order = order(as.numeric(gsub("CSNK1G", "", colnames(heatmat.lfq.log)))),
                        row_split = heatmap.levels.lfq$source.x,
                        column_split = rep(c("CSNK1G1", "CSNK1G2", "CSNK1G3"), each = 3),
                        width = unit(6, "cm"), heatmap_height = unit(25, "cm"), 
                        col = colorRamp2(c(0, 18, 19.3, 20.6, 21.9, 23.1, 14.4, 15.7, 29, 30), 
                                         c("#DDDDDD", "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b")),
                        right_annotation = c(track, wnt.track))

################################################################################
#Cytoscape Prey-Prey

heatmap.cyto <- heatmap.lfq %>%
  left_join(st, by = "Prey.Gene.Name") %>%
  select(Prey.Gene.Name, Prey.GeneID, Bait.Gene.Name, Bait.GeneID,Nucleus,Cytoplasm,Cytoskeleton,Endosome,
         ER,Extracellular,Golgi,Lysosome,Mitochondria,Peroxisome,Plasma_Membrane,Vesicles,Cell_Junction)%>%
  filter(Bait.Gene.Name %in% c("CSNK1G1", "CSNK1G2", "CSNK1G3")) %>%
  mutate(source = "MS/MS") %>%
  unique()

write_csv(heatmap.cyto, 'output/CSNK1G/cytoscape.csv')


#FIX THIS
heatmap.lfq.cyto <- heatmap.cyto %>%
  select(Prey.Gene.Name, Prey.GeneID, Bait.Gene.Name, Bait.GeneID)%>%
  filter(Bait.Gene.Name %in% c("CSNK1G1", "CSNK1G2", "CSNK1G3")) %>%
  unique()

###################################
#BioGrid Annotations

#Biogrid read  in
biogrid <- read.csv(file = 'C:/Users/smaranda/Documents/SmarandaSolomon/BIOGRID/BIOGRID-MV-Physical-4.4.203.tab3.txt', header = TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606)

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

#CHANGE AS.NUMERIC AND PUT IT ABOVE SO IT GOES FASTER
heatmap.lfq.cyto$in_BioGRID <- apply(heatmap.lfq.cyto, 1, function(x) {
  nrow(filter(interactions, (`interactor_min` == as.numeric(x[["Prey.GeneID"]]) & `interactor_max` == as.numeric(x[["Bait.GeneID"]])) |
                (`interactor_max` == as.numeric(x[["Prey.GeneID"]]) & `interactor_min` == as.numeric(x[["Bait.GeneID"]]))
  )) > 0
}) 

################################################################################
#Create separate data tables for prey-prey interactions
  
  #Filtering for interactions
  prey.prey.inter <- filter(interactions, (`interactor_min` %in%  heatmap.lfq.cyto$Prey.GeneID), (`interactor_max` %in%  heatmap.lfq.cyto$Prey.GeneID)) %>%
    rename(Prey.1.Entrez.ID = interactor_min, Prey.2.Entrez.ID = interactor_max) %>%
  mutate(Prey.1.Entrez.ID = as.character(Prey.1.Entrez.ID), 
         Prey.2.Entrez.ID = as.character(Prey.2.Entrez.ID))

  #Left_joining table with prey 1
  prey.prey.join <- left_join(prey.prey.inter, heatmap.lfq.cyto, by=c("Prey.1.Entrez.ID" = "Prey.GeneID")) %>%
    select(Prey.1.Entrez.ID, Prey.2.Entrez.ID)
  
  #Left_joining table with prey 2 and adding Uniprot and Nice Prey name columns
  prey.prey.final <- left_join(prey.prey.join, heatmap.lfq.cyto, by=c("Prey.2.Entrez.ID" = "Prey.GeneID"))%>%
    select(Prey.1.Entrez.ID, Prey.2.Entrez.ID) %>%
    mutate(source = "BioGRID") %>%
    unique()
  
  #write individual csv files for each bait
  write_csv(prey.prey.final, 'output/CSNK1G/prey_prey.csv')

################################################################################
#Enrichment Analysis
  
  library(clusterProfiler)
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("clusterProfiler")
  
  # reading in input
  univ <- read.csv(file='output/Annotated_Merge_NO_FILTER.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)
  saint <- read.csv(file='output/Annotated_Merge_Saint_filter.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
    dplyr::filter(grepl('CSNK1G1|CSNK1G2|CSNK1G3', Bait.Gene.Name))
  
  # 
  original_gene_list <- univ$Prey.GeneID
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  gene_list <- unique(gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  # Exctract significant results
  sig_genes_df = saint$Prey.GeneID
  
  # omit NA values
  genes <- na.omit(sig_genes_df)
  
  genes <- unique(genes)
  
  
  # xx <- as.list(org.Hs.egGO2ALLEGS)
  # 
  # idx <- tibble()
  # 
  # for(i in 1:length(xx)) {
  #   GO <- xx[[i]]
  #   row <- tibble(GeneID = GO, Source = names(GO), GO.Annotation = names(xx)[i])
  #   filter_row <- filter(row, Source %in% c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP"))
  #   summarise(
  #     #drop evidence columns and take unique rows
  #   )
  #   idx <- rbind(idx, filter_row)
  # }
  
  
  
  go_enrich <- enrichGO(gene = genes,
                        universe = gene_list,
                        OrgDb = org.Hs.eg.db, 
                        keyType = 'ENTREZID',
                        readable = T,
                        ont = "CC",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

  write_csv(go_enrich@result, 'output/go_enrichment.csv')