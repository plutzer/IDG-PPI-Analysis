

# library("devtools")
# devtools::install_github("smarasolo/cRomppass")

library(reticulate)
library(cRomppass)


# Set paths
SAINT_path = 'C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/build/SAINTexpress-spc.exe'
ED_path = 'C:/Users/plutzer/Work/IDG_pipeline/ED_DB.csv'
PG_path = 'C:/Users/plutzer/Work/IDG_pipeline/proteinGroups.txt'
output_dir = 'C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_blank'

setwd(output_dir)

py_run_string("print(\"Python is Running.\")")

# Generate SAINT inputs from old scoreAPMS python script
system(paste("python C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/score_APMS_noSAINT.py","--experimentalDesign",ED_path,"--proteinGroups",PG_path,"--outputPath",output_dir))

interaction_path = paste(output_dir,'/interaction.txt',sep = '')
prey_path = paste(output_dir,'/prey.txt',sep = '')
bait_path = paste(output_dir,'/bait.txt',sep = '')

# Run SAINT
system(paste(SAINT_path,interaction_path,prey_path,bait_path))

# Smaranda's CompPASS script
filename = paste(output_dir,"/to_CompPASS.csv",sep='')
to_comp = read.csv(file = filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

comp_out = comppass(to_comp, stats = NULL, norm.factor = 0.98)
write.table(comp_out, file=paste(output_dir,'/compPASS.csv',sep=''), sep = "\t")







