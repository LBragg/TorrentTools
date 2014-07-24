load(paste(parentDir, "screenedData/all_base_obs_filtered.RData", sep=""))

library(data.table)

myTab = tmp_dt[! is.na(FlowVal), length(ReadID), by="FlowVal,RefLen,FCYCLE,FCYCNUM"]
myTab$FCYCLE <- as.numeric(as.character(myTab$FCYCLE))
myTab$FCYCNUM <- as.numeric(as.character(myTab$FCYCNUM))
save(myTab, file=paste(parentDir, "empirical_dist_of_flow_values.RData", sep=""))



