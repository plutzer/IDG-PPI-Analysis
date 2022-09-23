library(ggplot2)
library(dplyr)
library(rlang)
library(ggfortify)
library(plotly)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(pROC)
library(reshape2)

PG_path = "C:/Users/plutzer/Work/IDG_pipeline/proteinGroups.txt"
annotated_merged_path = "C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_blank_2/Annotated_Merged_Output.csv"
annotated_merged_filtered_path = "C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_blank_2/Annotated_Filtered_Merged_Output.csv"

output_dir = 'C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_blank_2'

dir.create(paste0(output_dir,"/Plots"))


# ROC curve with known interactors from biogrid
    # Adapted for colocalization

annotated = read.csv(annotated_merged_path)

saintscore = roc(annotated,"in.BioGRID","SaintScore")
wd_fdr = roc(annotated,"in.BioGRID","perm_fdrs")
wd_raw = roc(annotated,"in.BioGRID","WD")
z_score = roc(annotated,"in.BioGRID","WD")

plot_rocs = plot(saintscore,col = "black",xlim=c(1,0),main = "ROC Curve - BioGRID Interactors")
plot_rocs = plot(wd_fdr,add=T,col = "red")
plot_rocs = plot(wd_raw,add=T,col = "blue")

legend("bottomright",
       legend=c("Saint Score", "WD FDR", "WD Raw"),
       col=c("black", "red", "blue"),
       lwd=8, cex =0.8, xpd = TRUE, horiz = TRUE)
################################################################################

# PCA with raw data
raw_data = read.csv(PG_path,sep='\t')
##### Any filtering of the raw data goes here





intensity_pca_data = raw_data %>% select(contains("Intensity.") & !contains("LFQ"))

# Rename columns
names(intensity_pca_data) <- sub('^Intensity.', '', names(intensity_pca_data))

# Need to remove rows with all 0's
intensity_pca_data = intensity_pca_data[rowSums(intensity_pca_data[])>0,]

intensity_pca = prcomp(t(intensity_pca_data),center = T,scale. = T)
autoplot(intensity_pca,label = F)

# Repeat with spectral counts
spc_pca_data = raw_data %>% select(contains("MS.MS.count."))

# Rename
names(spc_pca_data) <- sub('^MS.MS.count.', '', names(spc_pca_data))

# remove 0s rows
spc_pca_data = spc_pca_data[rowSums(spc_pca_data[])>0,]

spc_pca = prcomp(t(spc_pca_data),center = T,scale. = T)
autoplot(intensity_pca,label = T, label.size = 1.8)


################################################################################

# Volcano plots
# Saint Score and SAINT BFDR
annotated_merged = read.csv(annotated_merged_path)

saint_volc_data = annotated_merged %>% select(c("Bait.Gene.Name","Prey.Gene.Name","FoldChange","BFDR")) %>%
   group_by(Bait.Gene.Name,Prey.Gene.Name) %>%
  summarize(BFDR = mean(BFDR),
            FoldChange = mean(FoldChange)) %>% ungroup() %>%
  mutate(label = paste(Bait.Gene.Name," x ",Prey.Gene.Name))

# Remove baits with no gene name
saint_volc_data = saint_volc_data %>% filter(!(is.na(Bait.Gene.Name) | Bait.Gene.Name == ""))

ggplot(saint_volc_data,aes(x = log2(FoldChange),y = -log10(BFDR))) + geom_point()


wd_volc_data = annotated_merged %>% select(c("Bait.Gene.Name","Prey.Gene.Name","WD","perm_fdrs")) %>%
  group_by(Bait.Gene.Name,Prey.Gene.Name) %>%
  summarize(perm_fdrs = mean(perm_fdrs),
            WD = mean(WD)) %>% ungroup() %>%
  mutate(label = paste(Bait.Gene.Name," x ",Prey.Gene.Name))

# Remove baits with no gene name
wd_volc_data = wd_volc_data %>% filter(!(is.na(Bait.Gene.Name) | Bait.Gene.Name == ""))

ggplot(wd_volc_data,aes(x = log2(WD),y = -log10(perm_fdrs))) + geom_point()

## Iterate the volcano plots for each bait


#### trying plotly
p <- plot_ly(data = saint_volc_data, x = log2(saint_volc_data$FoldChange), y = -log10(saint_volc_data$BFDR), text = saint_volc_data$label, mode = "markers") %>% 
  layout(title ="SAINT Volcano Plot")
p

# GSEA by bait
#################
dir.create(paste0(output_dir,"/Plots/GSEA"))

for (bait in unique(annotated_merged$Bait.Gene.Name)) {
  
  bait = "CSNK1G1"
  
  print(bait)
  interactors = annotated_merged %>% filter(Bait.Gene.Name == bait) %>%
    filter(!(is.na(Prey.Gene.Name) | Prey.Gene.Name == "")) # Remove preys that don't have a gene name
  # Should already be sorted on WD score
  genelist = log2(interactors$FoldChange)
  names(genelist) = interactors$Prey.Gene.Name
  
  genelist = sort(genelist,decreasing = TRUE)
  
  gsea_bp = gseGO(geneList=genelist, 
               ont ="BP", 
               keyType = "SYMBOL",
               nPermSimple = 10000,
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "none")
  
  dir.create(paste0(output_dir,"/Plots/GSEA/",bait))
  dot_fname_bp = paste0(output_dir,"/Plots/GSEA/",bait,"/GSEA_dot_BP.png")
  gsea_fname_bp = paste0(output_dir,"/Plots/GSEA/",bait,"/GSEA_BP.png")
  

  dotplot(gsea_bp, showCategory=10, split=".sign") + facet_grid(.~.sign)

  

  gseaplot(gsea_bp, by = "all", title = gsea_bp$Description[1], geneSetID = 1)

  
  gsea_cc = gseGO(geneList=genelist, 
                  ont ="CC", 
                  keyType = "SYMBOL",
                  nPermSimple = 10000,
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "none")
  
  dir.create(paste0(output_dir,"/Plots/GSEA/",bait))
  dot_fname_cc = paste0(output_dir,"/Plots/GSEA/",bait,"/GSEA_dot_CC.png")
  gsea_fname_cc = paste0(output_dir,"/Plots/GSEA/",bait,"/GSEA_CC.png")
  
  

  dotplot(gsea_cc, showCategory=10, split=".sign") + facet_grid(.~.sign)

  

  gseaplot(gsea_cc, by = "all", title = gsea_cc$Description[1], geneSetID = 1)

  
}  

#################
# Bait-bait similarity metrics
# Jaccard index interactors
filtered_data = read.csv(annotated_merged_filtered_path)

get_jaccard = function(filtered_data,bait_1,bait_2) {
  bait_1_set = unique(filter(filtered_data,Bait.Gene.Name == bait_1)$Prey.Gene.Name)
  bait_2_set = unique(filter(filtered_data,Bait.Gene.Name == bait_2)$Prey.Gene.Name)
  # Compute the Jaccard Index
  intersection = length(intersect(bait_1_set,bait_2_set))
  union = length(bait_1_set) + length(bait_2_set) - intersection
  return (intersection/union)
}

baits = unique(annotated_merged$Bait.Gene.Name)[! unique(annotated_merged$Bait.Gene.Name) %in% c("")]
num_baits = length(baits)

sim_matrix = matrix(nrow = num_baits,ncol = num_baits)
rownames(sim_matrix) = baits
colnames(sim_matrix) = baits
for (i in 1:num_baits) {
  for (j in 1:num_baits) {
    sim_matrix[i,j] = get_jaccard(filtered_data,baits[i],baits[j])
  }
}

heatmap(sim_matrix,symm=T)


  # Summary of correlations

testx = c(1,2,3,4,5)
testy = c(5,3,1,3,4)

test_plot <- ggplot(data.frame(testx,testy),aes(testx,testy)) + 
  geom_line()
test_plot
