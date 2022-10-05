
repo_path = getwd()

library(reticulate)
library(cRomppass)
library(tidyverse)
library(DarkKinaseTools)
library(org.Hs.eg.db)
library(dplyr)
#source("compute_wds.R")
source("filter_interaction.R")
source("comppass-FDR.R")
source("cytoscape_outputs.R")
source("merge_annotate_filter.R")

################## ARGUMENTS ####################################################
args = commandArgs(trailingOnly = TRUE)

# Check arguments

################# Parameters ###################################################
ED_path = args[1]
PG_path = args[2]
output_dir = args[3] ## No slash at end!!!

resampling_iterations=as.integer(args[4]) #for perm_fdr calculation

quantification_method = args[5]

# Filtering parameters
BFDR_cutoff = 0.05
AvgP_cutoff = 0.7
Comppass_percent = 0.05 # Might not need?
################################################################################

############### Setting variables ##############################################
uniprot_map_path = paste0(repo_path,'/uniprot_mapping.tsv.zip')
biogrid_mv_path = paste0(repo_path,'/BIOGRID-MV-Physical-4.4.211.tab3.txt')
biogrid_all_path = paste0(repo_path,'/BIOGRID-ALL-4.4.211.tab3.txt')

# Set the SAINT path based on the quantification method
if (quantification_method == "spc") {
  SAINT_path = paste0(repo_path,"/build/SAINTexpress-spc.exe")
} else if (quantification_method == "intensity") {
  SAINT_path = paste0(repo_path,"/build/SAINTexpress-int.exe")
} else {print("Invalid argument: quantification_method")}

# This is for some reason necessary
select <- get(x="select", pos = "package:dplyr")

################################################################################

########### Generate SAINT inputs ##############################################
# For some unknown reason, this line is necessary even though it does nothing.
py_run_string("print(\"Python is Running.\")")

# Re-using the exisiting code for this because this parser works well
system(paste0("python ",repo_path,"/score_APMS_noSAINT.py", # Change filename
             " --experimentalDesign ",ED_path,
             " --proteinGroups " ,PG_path,
             " --outputPath ",output_dir,
             " --quantification-saint ",quantification_method,
             " --quantification-comppass ",quantification_method))

# Set variables for new files generated
interaction_path = paste(output_dir,'/interaction.txt',sep = '')
prey_path = paste(output_dir,'/prey.txt',sep = '')
bait_path = paste(output_dir,'/bait.txt',sep = '')

# Filter the interaction file to remove preys with 0 intensity or spectral counts in all reps
filtered_interaction_path = filter_interaction(interaction_path)

# Run SAINT with defaults
setwd(output_dir)
system(paste(SAINT_path,"-L 2",filtered_interaction_path,prey_path,bait_path))
setwd(repo_path)
################################################################################

############# Run comppass and Perm FDR ########################################
to_comp_filename = paste(output_dir,"/to_CompPASS.csv",sep='')

comp_out = run_comppass(to_comp_filename,n_iter = resampling_iterations)
################################################################################

# Merge the compass and SAINT outputs
merged = merge_scores(paste(output_dir,'/list.txt',sep=''),comp_out,output_dir)

# Annotated the merged data with Uniprot and GO info
annotated = annotate_uniprot_go(merged,uniprot_map_path)

# Check to see if a pre-processed biogrid file exists...
if (file.exists("biogrid_summary.csv")) {
  print("Found pre-processed BioGRID summary.")
  summ_biogrid = read.csv(file = "biogrid_summary.csv")
} else {
  print("No pre-processed BioGRID file found.")
  print("Processing BioGRID file. This may take several minutes.")
  summ_biogrid = preprocess_biogrid(biogrid_all_path,biogrid_mv_path,output_dir)
}

# Merge the data with the biogrid summary file.
summ_biogrid$Entrez.Gene.Interactor.A = as.character(summ_biogrid$Entrez.Gene.Interactor.A)
summ_biogrid$Entrez.Gene.Interactor.B = as.character(summ_biogrid$Entrez.Gene.Interactor.B)
annotated = left_join(annotated,summ_biogrid, by=c('Bait.GeneID'='Entrez.Gene.Interactor.A','Prey.GeneID'='Entrez.Gene.Interactor.B'))

# Fix this column for ROC analysis
annotated$in.BioGRID[is.na(annotated$in.BioGRID)] = FALSE
annotated$Multivalidated[is.na(annotated$Multivalidated)] = FALSE

#################### Filtering #################################################
# Identify test baits for filtering
test_baits = read.csv(ED_path, header = TRUE, stringsAsFactors = FALSE) %>%
  filter(Type == 'T')
annotated = annotated %>% filter(Experiment.ID %in% test_baits$Bait)

# Sort the file by WD score
annotated <- arrange(annotated, desc(WD))

# Write out the unfiltered merged and annotated file.
write.csv(annotated, paste(output_dir,"/Annotated_Merged_Output.csv",sep=''),na="",row.names=F)

# Filter based on the cutoffs specified above.
annotated_filtered = annotated %>% filter(BFDR <= BFDR_cutoff, AvgP >= AvgP_cutoff)

# Write filtered file
write.csv(annotated_filtered, paste(output_dir,"/Annotated_Filtered_Merged_Output.csv",sep=''),na="",row.names=F)
################################################################################

###################### Cytoscape Outputs #######################################
create_cytoscape_outputs(annotated,annotated_filtered,biogrid_mv_path,output_dir)
################################################################################




