library("devtools")

args <- commandArgs(trailingOnly = TRUE)
fn <- "output/to_CompPASS.csv"
df <- read.csv(file = fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

exit = cRomppass::comppass(df, stats = NULL, norm.factor = 0.98)

write.table(exit, file='output/compPASS.csv', sep = "\t")
