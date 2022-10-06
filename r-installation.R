install.packages("devtools")
library("devtools")

devtools::install_github("smarasolo/cRomppass")
 
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

install.packages("reticulate")

devtools::install_github("IDG-Kinase/DarkKinaseTools")

