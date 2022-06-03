# Title     : heatmap_DB_Ck1g
# Objective : Create a heatmap using Megan & Dhaval new wnt stuff
# Created by: Smaranda Solomon
# Created on: 2/2/2022

library(devtools)
library(tidyverse)
library(circlize)

select <- get(x="select", pos = "package:dplyr")

prey <- read.csv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Results/CSNK1G/DB/Annotated_Merge_Saint_filter.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE) %>%
  dplyr::filter(grepl('CSNK1G3', Bait.Gene.Name))

just.wnt <- prey %>% 
  dplyr::filter(grepl('AAK1|AMER1|ANKRD6|APC|AXIN1|AXIN2|CELSR1|CELSR2|CSNK1A|CSNK1D|CSNK1E|CSNK1G1|CSNK1G2|CSNK1G3|CTNNB1|DVL1|DVL2|DVL3|FAT|FZD1|FZD2|FZD3|FZD4|FZD5|FZD6|FZD7|FZD8|
                      GSK3A|GSK3B|LRP5|LRP6|NOTCH|ROR1|ROR2|RYK|VANGL', Prey.Gene.Name)) %>%
  select(Experiment.ID, Bait.Gene.Name, Prey.Gene.Name, Spec) 

just.new.wnt <- prey %>% 
  dplyr::filter(grepl('FAT4|FZD3|LRP6|NOTCH3|ROR2|ROR1', Prey.Gene.Name)) %>%
  select(Experiment.ID, Bait.Gene.Name, Prey.Gene.Name, Spec) %>%
  unique()

################################################################################
#Heatmap made from LFQ intensities

lfq <- read_tsv(file='C:/Users/smaranda/Documents/SmarandaSolomon/Experimental_Data/20211014_CK1G_All_Baits/proteinGroups.txt', guess_max = 10000) %>%
  select(`Gene names`,`LFQ intensity CSNK1G1-30mB-MA_1`, `LFQ intensity CSNK1G1-30mB-MA_2`, `LFQ intensity CSNK1G1-30mB-MA_3`, 
         `LFQ intensity CSNK1G2-30mB-MA_1`, `LFQ intensity CSNK1G2-30mB-MA_2`, `LFQ intensity CSNK1G2-30mB-MA_3`,
         `LFQ intensity CSNK1G3-30mB-MA_1`, `LFQ intensity CSNK1G3-30mB-MA_2`, `LFQ intensity CSNK1G3-30mB-MA_3`,
         `LFQ intensity mT-CSNK1G3_1`, `LFQ intensity mT-CSNK1G3_2`)

heatmap.lfq <- left_join(just.wnt, lfq, by = c("Prey.Gene.Name" = "Gene names")) %>%
  select(-Spec, -Bait.Gene.Name, -Experiment.ID) %>%
  unique()

heatmap.lfq <- heatmap.lfq %>%
  filter(rowSums(across(matches("LFQ"))) > 0)

bait.order <- order(rowSums(heatmap.lfq[,-1]),decreasing=T)

heatmap.lfq <- heatmap.lfq[bait.order,]

#################################################
#Matrix and heatmap with all WNT

heatmap.mat.lfq <- as.matrix(heatmap.lfq) 

heatmap.mat.lfq <- heatmap.mat.lfq[,-1]

heatmap.mat.lfq <- matrix(as.numeric(unlist(heatmap.mat.lfq)), nrow = nrow(heatmap.mat.lfq))

rownames(heatmap.mat.lfq) <- heatmap.lfq$Prey.Gene.Name

heatmat.lfq.log <-  log2(heatmap.mat.lfq + 1)

colnames(heatmat.lfq.log) <- c("CSNK1G1_1", "CSNK1G1_2", "CSNK1G1_3", "CSNK1G2_1", "CSNK1G2_2", 
                               "CSNK1G2_3", "CSNK1G3_1", "CSNK1G3_2", "CSNK1G3_3",
                               "CSNK1G3_1", "CSNK1G3_2")

ComplexHeatmap::Heatmap(heatmat.lfq.log,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        column_order = order(as.numeric(gsub("CSNK1G", "", colnames(heatmat.lfq.log)))),
                        column_split = c(rep("30 min", 9), rep("120 min", 2)),
                        col = colorRamp2(c(0, 17, 30), c("black","white", "red")),
                        width = unit(10, "cm"), heatmap_height = unit(12, "cm"))

#################################################
#Matrix and heatmap with all NEW WNT

heatmap.lfq.new <- left_join(just.new.wnt, lfq, by = c("Prey.Gene.Name" = "Gene names")) %>%
  select(-Spec, -Bait.Gene.Name, -Experiment.ID) %>%
  unique()

heatmap.lfq.new <- heatmap.lfq.new %>%
  filter(rowSums(across(matches("LFQ"))) > 0)

bait.order <- order(rowSums(heatmap.lfq.new[,-1]),decreasing=T)

heatmap.lfq.new <- heatmap.lfq.new[bait.order,]

heatmap.mat.new <- as.matrix(heatmap.lfq.new) 

heatmap.mat.new <- heatmap.mat.new[,-1]

heatmap.mat.new <- matrix(as.numeric(unlist(heatmap.mat.new)), nrow = nrow(heatmap.mat.new))

rownames(heatmap.mat.new) <- heatmap.lfq.new$Prey.Gene.Name

heatmat.new.log <-  log2(heatmap.mat.new + 1)

colnames(heatmat.new.log) <- c("CSNK1G1_1", "CSNK1G1_2", "CSNK1G1_3", "CSNK1G2_1", "CSNK1G2_2", 
                               "CSNK1G2_3", "CSNK1G3_1", "CSNK1G3_2", "CSNK1G3_3",
                               "CSNK1G3_1", "CSNK1G3_2")

ComplexHeatmap::Heatmap(heatmat.new.log,
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        column_order = order(as.numeric(gsub("CSNK1G", "", colnames(heatmat.new.log)))),
                        column_split = c(rep("30 min", 9), rep("120 min", 2)),
                        col = colorRamp2(c(0, 18, 26), c("black","white", "red")),
                        width = unit(10, "cm"), heatmap_height = unit(12, "cm"))

