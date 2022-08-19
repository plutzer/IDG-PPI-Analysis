#source("C:/Users/plutzer/Repos/IDG-PPI-Analysis_plutzer/compPASS.R")

library(cRomppass)

to_comp_test = read.csv(file = "C:/Users/plutzer/Work/IDG_pipeline/outputs/testset_int/to_CompPASS.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

start = Sys.time()
comp_out = comppass(to_comp_test, stats = NULL, norm.factor = 0.98)
end = Sys.time()
start-end

permute = function(dataframe) {
  perm = to_comp_test[order(to_comp_test$Bait),]
  perm$new_spec = perm %>% group_by(Bait) %>% group_map(~ sample(.x$Spectral.Count,replace = TRUE)) %>% unlist()
  perm
}

# Make a list of permutation dataframes
n_perms = 5



perms = lapply(rep(to_comp_test,n_perms),permute)


# perm_test = permute(to_comp_test)

# Check to see if the permutation worked
# for (testbait in unique(perm_test$Bait)) {
#   bait_subset = perm_test %>% filter(Bait == testbait)
#   for (item in bait_subset$new_spec) {
#     if (!(item %in% bait_subset$Spectral.Count)){
#       print("OOPS")
#     }
#   }
# } LOOKS GOOD!


permuted_wd = 
for (i in c(1:10)) {
  comp_in_i = comp_out
  comp_out_i = comppass(to_comp_test, stats = NULL, norm.factor = 0.98)
}


testgroup = group_by(to_comp_test,Bait)


DqStr <- "Group   q        Dq       SD.Dq
1 -3.0 0.7351 0.0067
1 -2.5 0.6995 0.0078
1 -2.0 0.6538 0.0093
2 -3.0 0.7203 0.0081
2 -2.5 0.6829 0.0094
2 -2.0 0.6350 0.0112"
toy <- read.table(textConnection(DqStr), header=TRUE)

toy_group = toy %>% group_by(Group)

toy_group$q = sample(toy_group$q)

new_col = toy %>% group_by(Group) %>% group_map(~ sample(.x$q)) %>% unlist()


permutation = to_comp_test %>% group_by(Bait) %>% group_modify(~ sample(.x$Spectral.Count)) %>% ungroup()

to_comp_test$perm = permutation

subset1 = to_comp_test %>% filter(Bait == "O14874") # This didn't work. There are numbers here that didn't come from the correct group
