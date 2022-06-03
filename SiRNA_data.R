# Title     : SiRNA_data
# Objective : Align the SiRNA data for WNT signaling with the scored SAINT and CompPass file
# Created by: Smaranda Solomon
# Created on: 5/26/2021

sirna.data <- function(){
  
library(tidyverse)
library(org.Hs.eg.db)
library(limma)

select <- get(x="select", pos = "package:dplyr")

#Load files
st <- read.csv(file='annotations/Madan_Secondary.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

# neg <- read.csv(file='annotations/Lebenson_WNT_negative_reg.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)
# pos <- read.csv(file='annotations/Lebenson_WNT_positive_reg.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)
# att <- read.csv(file='annotations/Lebenson_WNT_attenuating_reg.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

neg <- read.csv(file='annotations/Lebenson_negative.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)
pos.high <- read.csv(file='annotations/Lebenson_pos_high_str.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)
pos.low <-  read.csv(file='annotations/Lebenson_pos_low_str.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)
att <- read.csv(file='annotations/Lebenson_attenuating.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

biechele <- read.csv(file='annotations/Biechele_et_al.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)
 
AGGF <- read.csv(file='annotations/AGGF_screens.csv', header = TRUE, sep = ",", stringsAsFactors = FALSE)

################################################################################
#Hits data sheet (Madan_Secondary)
#2/3 of the singles have to be a hit for it to be marked as hit and each hit has to be two fold (>50 or <200)

data_filter <- filter(st, siRNA.Type != "Pool")

hits.table <- data_filter %>%
  group_by(Gene.Symbol, Cell_Line) %>%
  group_modify(~ {.x %>%
      mutate(hits =  case_when(
        sum(Median.Norm.FF...Luc. > 200) >= 2 ~ "Up",
        sum(Median.Norm.FF...Luc. < 50) >= 2 ~ "Down"
        ))
      }) %>%
  summarize(Cell_Line = unique(Cell_Line), hits = unique(hits))

#Code to find Gene Symbol Alias and attach to table
names <- alias2SymbolTable(unique(hits.table$Gene.Symbol), species = "Hs") %>%
  data.frame()

to.be.joined <- hits.table %>%
  ungroup() %>%
  select("Gene.Symbol") %>%
  unique() %>%
  cbind(names) %>%
  rename(. = "Gene.Alias.Hits")

aligned.hits <- left_join(hits.table, to.be.joined, by = "Gene.Symbol") %>%
  select("Gene.Symbol", "Gene.Alias.Hits", "Cell_Line", "hits")

################################################################################
#Second data (Lebensohn)
#Up, down, att is true if significant (p<0.01), false if not significant (p>0.01)

neg.table <- neg %>%
  within(Neg <- ifelse(FDR.corrected.p.value <= 0.01, TRUE, FALSE))

#pos.table <- pos %>%
# within(Pos <- ifelse(FDR.corrected.p.value < 0.01, TRUE, FALSE))

pos.high.table <- pos.high %>%
  within(Pos.high <- ifelse(FDR.corrected.p.value <= 0.01, TRUE, FALSE))

pos.low.table <- pos.low %>%
  within(Pos.low <- ifelse(FDR.corrected.p.value <= 0.01, TRUE, FALSE))

att.table <- att %>%
  within(Att <- ifelse(FDR.corrected.p.value <= 0.01, TRUE, FALSE))

# leben.table <- neg.table %>%
#   full_join(pos.table, by = "Gene") %>%
#   full_join(att.table, by = "Gene") %>%
#   select(Gene, Neg, Pos, Att)

leben.table <- neg.table %>%
  full_join(pos.high.table, by = "Gene") %>%
  full_join(pos.low.table, by = "Gene") %>%
  full_join(att.table, by = "Gene") %>%
  select(Gene, Neg, Pos.high, Pos.low, Att)

#Code to find Gene Symbol Alias and attach to table
al <- alias2SymbolTable(unique(leben.table$Gene), species = "Hs") %>%
  data.frame()

name.join <- leben.table %>%
  select("Gene") %>%
  unique() %>%
  cbind(al) %>%
  rename(. = "Gene.Alias.Reg")

aligned.leben <- left_join(leben.table, name.join, by = "Gene")
#aligned.leben <- aligned.leben[,c(1, 5, 2:4)]
aligned.leben <- aligned.leben[,c(1, 6, 2:5)]

################################################################################
#Third data

biechele.data <- biechele %>%
  mutate(Biechele.Fold = ifelse(Fold.over.control.siRNA <= 0.5, TRUE, ifelse( Fold.over.control.siRNA >= 2, TRUE, FALSE)))%>%
  select(-Average.BAR.firefly.Luciferase..Raw., -Fold.over.control.siRNA)

################################################################################
#Fourth data, AGGF
#2/3 of the singles have to be a hit for it to be marked as hit and each hit has to be two fold (>50 or <200)

aggf_filter <- filter(AGGF, siRNA.Type != "Pool")

aggf.table <- aggf_filter %>%
  group_by(Gene.Symbol, Cell_Line) %>%
  group_modify(~ {.x %>%
      mutate(hits =  case_when(
        sum(P.val.Norm.FF...Luc. <= .05) >= 2 ~ "AGGF.Hit"
      ))
  }) %>%
  #summarize(Cell_Line = unique(Cell_Line), AGGF.Hits = unique(hits))
  group_by(Gene.Symbol) %>%
  summarize(AGGF.Hit = any(hits == "AGGF.Hit", na.rm = TRUE))

#Code to find Gene Symbol Alias and attach to table
# names.aggf <- alias2SymbolTable(unique(aggf.table$Gene.Symbol), species = "Hs") %>%
#   data.frame()
# 
# aggf.join <- aggf.table %>%
#   ungroup() %>%
#   select("Gene.Symbol") %>%
#   unique() %>%
#   cbind(names) %>%
#   rename(. = "Gene.Alias.Hits")
# 
# aggf.final <- left_join(names.aggf, aggf.join, by = "Gene.Symbol") %>%
#   select("Gene.Symbol", "Gene.Alias.Hits", "Cell_Line", "hits")


################################################################################
#Joining

combined.start <- full_join(aligned.hits, aligned.leben, by = c("Gene.Symbol" = "Gene"))

combined.middle <- full_join(combined.start, biechele.data, by = c("Gene.Symbol" = "siRNA.Target"))

combined.final <- full_join(combined.middle, aggf.table)

write.csv(combined.final, 'output/Screens_for_PIN.csv')

return(combined.final)
}