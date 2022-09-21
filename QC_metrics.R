library(ggplot2)
library(dplyr)
library(rlang)
library(ggfortify)

PG_path = "C:/Users/plutzer/Work/IDG_pipeline/proteinGroups.txt"
annotated_merged_path = "C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_blank_2/Annotated_Merged_Output.csv"

# ROC curve with known interactors from biogrid
    # Adapted for colocalization
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

# GSEA by bait

# Bait-bait similarity metrics
  # Jaccard index interactors
  # Summary of correlations

testx = c(1,2,3,4,5)
testy = c(5,3,1,3,4)

test_plot <- ggplot(data.frame(testx,testy),aes(testx,testy)) + 
  geom_line()
test_plot
