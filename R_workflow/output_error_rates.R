load(paste(parentDir, "screenedData/all_base_obs_filtered.RData", sep=""))
nas <- which(is.na(dt_full$FlowPos))

library(plyr)
library(data.table)
setkey(dt_full, FlowPos, ReadID)



resRefLenZero = function(dt)
{
  ids = unique(dt$ReadID)   # or from the data: unique(ds$ID)
  pos = unique(dt$FlowPos[which(dt$FlowPos >= 0)])   # or from the data: unique(ds$Pos) 
  setkey(dt, ReadID,FlowPos)


  ## Trying to calculate the number of correct versus incorrect zeroes.

  #1. Number of calls at all flows (for total)
  calls_at_flows = dt[CJ(ids,pos), roll=-Inf, nomatch=0][, .N, by=FlowPos] #calls across all flow positions
  setkey(calls_at_flows, FlowPos)

  #2. Number of positive calls across flows

  pos_calls = dt[, list(pos_calls = length(ReadID)), FlowPos]
  setkey(pos_calls, FlowPos)

  #3. Subtract #2 from #1 and you've got the number of zero calls (correct)
  setkey(dt, RefLen)

  oc_zero = dt[J(0), list(count.above=length(ReadID)), FlowPos] 
  setkey(oc_zero, FlowPos)

  # Merge all these.
  by_base_zero = calls_at_flows[pos_calls[oc_zero]]
  by_base_zero[, count.at := N - pos_calls]
  return (by_base_zero)
}

## aggregated by RefLen and FlowPos.

#############################

setkey(dt_full, OOP)

res_in_phase = resRefLenZero(dt_full["False",])
res_in_phase[,  `:=`(count.below = 0, RefLen = 0),]
res_in_phase[, c("N","pos_calls") := NULL]


res_oop = resRefLenZero(dt_full["True",])
res_oop[,  `:=`(count.below = 0, RefLen = 0),]
res_oop[, c("N","pos_calls") := NULL]

#write.table(res_in_phase, file=paste(parentDir, "error_rate_for_zero_inphase.csv", sep=""))
#write.table(res_oop, file=paste(parentDir, "error_rate_for_zero_oophase.csv", sep=""))

##cheating.
byBase = function(dt)
{
    #count below count at count above
    return(dt[, list(count.below = sum(CallLen < RefLen), count.at = sum(CallLen == RefLen), count.above = sum(CallLen > RefLen)), by=list(RefLen, FlowPos)])
}

res_in_phase_non_zero = byBase(dt_full[RefLen > 0 & OOP == "False",])
res_oop_non_zero = byBase(dt_full[RefLen > 0 & OOP == "True",])

#write.table(res_in_phase, file=paste(parentDir, "error_rate_by_flowpos_reflen_in_phase.csv", sep=""))
#write.table(res_oop, file=paste(parentDir, "error_rate_by_flowpos_reflen_oop.csv", sep=""))

## combine these into one table.

res_in_phase_all = rbind(res_in_phase, res_in_phase_non_zero, use.names=TRUE)
setkey(res_in_phase_all, RefLen, FlowPos)
res_oop_all = rbind(res_oop, res_oop_non_zero, use.names=TRUE)
setkey(res_oop_all, RefLen, FlowPos)

#finally.

#write.table(res_in_phase_all, file=paste(parentDir, "error_rate_by_flowpos_inphase.csv", sep=""))
#write.table(res_oop_all, file=paste(parentDir, "error_rate_by_flowpos_oop.csv", sep=""))

## even better than above, combine these two in one data.table.

#res_in_phase_all = read.table("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/ecoli_no_calibration_test/error_rate_by_flowpos_inphase.csv")
#res_oop_all = read.table("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/ecoli_no_calibration_test/error_rate_by_flowpos_oop.csv")

res_in_phase_all$by = "OOP"
res_in_phase_all$state = "False"

res_oop_all$by = "OOP"
res_oop_all$state = "True"

res_all = rbind(res_in_phase_all, res_oop_all)

write.table(res_all, paste(parentDir, "error_rate_by_flowpos_phased.csv", sep=""))
