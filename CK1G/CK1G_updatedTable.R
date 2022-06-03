# Title     : CK1G_updatedTable
# Objective : Add the extra columns to the heatmap
# Created by: Smaranda Solomon
# Created on: 2/25/2022

library(devtools)
library(tidyverse)
library(circlize)

organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

select <- get(x="select", pos = "package:dplyr")
rename <- get(x="rename", pos = "package:dplyr")

################################################################################ 
#Read in Megan Data and filter

sec <- read.csv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Results/CSNK1G/3_Replicates/Annotated_Merge_All_filtered.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::filter(grepl('CSNK1G1|CSNK1G2|CSNK1G3', Bait.Gene.Name)) %>%
  select(Experiment.ID, Bait.Gene.Name, Prey.Gene.Name, BFDR, AvgP, WD, Spec,Nucleus,Cytoplasm,Cytoskeleton,Endosome,ER,Extracellular,Golgi,Lysosome,
         Mitochondria,Peroxisome,Plasma_Membrane,Vesicles,Cell_Junction)

st <- read.csv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Results/CSNK1G/3_Replicates/Annotated_Merge_Saint_filter.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

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
#Read in Dhaval data and filter

db <- read.csv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Results/CSNK1G/DB/Annotated_Merge_Saint_filter.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::filter(grepl('CSNK1G3', Bait.Gene.Name))

db.wnt <- db %>% 
  dplyr::filter(grepl('AAK1|AMER1|ANKRD6|APC|AXIN1|AXIN2|CELSR1|CELSR2|CSNK1A|CSNK1D|CSNK1E|CSNK1G1|CSNK1G2|CSNK1G3|CTNNB1|DVL1|DVL2|DVL3|FAT|FZD1|FZD2|FZD3|FZD4|FZD5|FZD6|FZD7|FZD8|
                      GSK3A|GSK3B|LRP5|LRP6|NOTCH|ROR1|ROR2|RYK|VANGL', Prey.Gene.Name)) %>%
  select(Bait.Gene.Name, Prey.Gene.Name, Spec) 

saint.comp.db <- full_join(saint.comp, db.wnt, by = "Prey.Gene.Name")
saint.comp.db <- mutate_at(saint.comp.db, "Bait.Gene.Name.x", ~replace(., is.na(.), 'CSNK1G3'))
saint.comp.db$Spec.x[is.na(saint.comp.db$Spec.x)] <- saint.comp.db$Spec.y[is.na(saint.comp.db$Spec.x)] 

saint.comp.db <- saint.comp.db %>%
  select(-Bait.Gene.Name.y, -Spec.y) %>%
  rename(Bait.Gene.Name = Bait.Gene.Name.x, Spec = Spec.x)

################################################################################
#Read in Screen Data and filter

screen <- read.csv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Results/CSNK1G/3_Replicates/Annotated_Merge_Saint_filter.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  select(Experiment.ID, Bait.Gene.Name, Prey.Gene.Name, Spec, hits, Neg, Pos.high, Pos.low, Att, Biechele.Fold, AGGF.Hit) %>%
  filter(grepl('CSNK1G1|CSNK1G2|CSNK1G3', Bait.Gene.Name)) %>%
  mutate(screen=TRUE)

screen <- screen %>%
  filter(!is.na(hits) | !is.na(Neg) | !is.na(Pos.high) | !is.na(Pos.low) | !is.na(Att) | !is.na(Biechele.Fold) | AGGF.Hit == TRUE) %>%
  filter(!is.na(hits) | Neg == TRUE | Pos.high == TRUE | Pos.low == TRUE | Att == TRUE | Biechele.Fold == TRUE | AGGF.Hit == TRUE) %>%
  unique()

extra.prey <- read.csv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Results/CSNK1G/3_Replicates/Annotated_Merge_NO_FILTER_3.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::filter(grepl('CSNK1G1|CSNK1G2|CSNK1G3', Bait.Gene.Name)) %>%
  dplyr::filter(Prey.Gene.Name %in% unique(saint.comp.db$Prey.Gene.Name)) %>%
  full_join(., screen)%>%
  select(Bait.Gene.Name, Prey.Gene.Name, Spec, screen) %>%
  unique()

extra.prey <- extra.prey %>%
  arrange(., Prey.Gene.Name, screen) %>%
  group_by(Prey.Gene.Name) %>%
  fill(screen)

write_csv(extra.prey, 'output/extra.prey.csv')

################################################################################
#Filter Screen data for tracks

saint.and.comp <- sec %>%
  select(Prey.Gene.Name) %>%
  mutate(source = "Saint.CompPASS") %>%
  unique()

wnt.and.saint <- just.saint %>%
  select(Prey.Gene.Name)
  
wnt.and.saint <- full_join(wnt.and.saint,  db.wnt, by = "Prey.Gene.Name") %>%
  select(Prey.Gene.Name) %>%
  mutate(source = "Saint.WNT") %>%
  unique()

screens.and.saint <- screen %>%
  select(Prey.Gene.Name) 

screens.and.saint <-  left_join(screens.and.saint,  db.wnt, by = "Prey.Gene.Name")%>%
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

#################################################
#Dhaval addition

levels.test <- full_join(levels, db.wnt, by = "Prey.Gene.Name") 
levels.test <- mutate_at(levels.test, "Bait.Gene.Name.x", ~replace(., is.na(.), 'CSNK1G3'))

levels.test$Spec.x[is.na(levels.test$Spec.x)] <- levels.test$Spec.y[is.na(levels.test$Spec.x)] 

levels.test$source.x[is.na(levels.test$source.x)] <- "Saint.WNT"

levels <- levels.test %>%
  select(Bait.Gene.Name.x, Prey.Gene.Name, Spec.x, screen, source.x) %>%
  rename(Bait.Gene.Name = Bait.Gene.Name.x, Spec = Spec.x)

#################################################

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
         `LFQ intensity CSNK1G3-30mB-MA_1`, `LFQ intensity CSNK1G3-30mB-MA_2`, `LFQ intensity CSNK1G3-30mB-MA_3`,
         `LFQ intensity mT-CSNK1G3_1`, `LFQ intensity mT-CSNK1G3_2`)

heatmap.lfq <- left_join(extra.prey, lfq, by = c("Prey.Gene.Name" = "Gene names")) %>%
  select(-Spec, -Bait.Gene.Name) %>%
  unique()

heatmap.lfq <- heatmap.lfq %>%
  filter(rowSums(across(matches("LFQ"))) > 0)

heatmap.lfq <- heatmap.lfq[,c(1, 3:13, 2)]

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

colnames(heatmat.lfq.log) <- c("CSNK1G1_1", "CSNK1G1_2", "CSNK1G1_3", "CSNK1G2_1", "CSNK1G2_2", "CSNK1G2_3",
                               "CSNK1G3_1", "CSNK1G3_2", "CSNK1G3_3", "CSNK1G3_1", "CSNK1G3_2")

#Color for screen track
colors <- list("screen" = c("TRUE" = "#C7706F", "FALSE" = "white"))

#Screen track command
track <- ComplexHeatmap::HeatmapAnnotation(df = track.names2, which = "row", col = colors, simple_anno_size = unit(.3, "cm"))

#Color for wnt track
colors.wnt <- list("source" = c("TRUE" = "#B0D656", "FALSE" = "white"))

#Wnt track command
wnt.track <- ComplexHeatmap::HeatmapAnnotation(df = wnt.and.saint.track, which = "row", col = colors.wnt, simple_anno_size = unit(.3, "cm"))

myheatmap <- ComplexHeatmap::Heatmap(heatmat.lfq.log,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        column_order = order(as.numeric(gsub("CSNK1G", "", colnames(heatmat.lfq.log)))),
                        row_split = heatmap.levels.lfq$source.x,
                        column_split = c(rep(c("CSNK1G1_30min", "CSNK1G2_30min", "CSNK1G3_30min"), each = 3), 
                                             rep("CSNK1G3_z", 2)),
                        width = unit(6, "cm"), heatmap_height = unit(26, "cm"),
                        col = colorRamp2(c(0, 18, 19.3, 20.6, 21.9, 23.1, 14.4, 15.7, 29, 30), 
                                         c("#DDDDDD", "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b")),
                        right_annotation = c(track, wnt.track))

gb_heatmap = grid::grid.grabExpr(ComplexHeatmap::draw(myheatmap), height=20, width=6)
ggsave("C:/Users/smaranda/Documents/SmarandaSolomon/Results/CSNK1G/heatmap.pdf", height=20, width=6, compress=F, gb_heatmap)
