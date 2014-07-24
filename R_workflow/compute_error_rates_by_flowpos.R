load(paste(parentDir, "screenedData/all_base_obs_filtered.RData", sep=""))

#need to do output error rate by RLE Length
## Need to probably completed re-write this. Think about what the outcome is first.

library(plyr)
library(data.table)

runStuff <- function(dt, fname, max, oop)
{
  ## need to calculate correct zeroes again??
  
	resRefLenZero = function(dt)
	{
    print("Running refResLenZero")
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
  
  zero_res = resRefLenZero(dt[ReadPos < max & OOP == oop,])
  setkey(zero_res, FlowPos)
	write.table(zero_res, file=paste(parentDir, "error_rate_for_zero_", fname, ".csv", sep=""))

	byBase = function(dt)
	{
    print("Running by base")
	  #count below count at count above
	  return(dt[, list(count.below = sum(CallLen < RefLen), count.at = sum(CallLen == RefLen), count.above = sum(CallLen > RefLen)), by=list(RefLen, FlowPos)])
	}
  
  dt_sub = dt[ReadPos < max & OOP == oop & RefLen > 0,]
  pos_res = byBase(dt_sub)
  setkey(pos_res, RefLen, FlowPos)
	write.table(pos_res, file=paste(parentDir, "error_rate_by_flowpos_reflen_", fname, ".csv", sep=""))
}


tmp_dt = tmp_dt[which(!is.na(tmp_dt$FlowPos)),] #problem is that we have to ignore the NA flows. #May want to consider in and out of sync reads separetely.

#requires maxBounds to be defined.
for(ph in c("True", "False"))
{
	ph_str = ifelse(ph == "True", "oop", "inphase")
	for(maxV in maxBounds)
	{
	  #Readpos <= maxBounds
		ph = "True"
		label <- paste("first", maxV, ph_str, sep="")
		runStuff(tmp_dt, label, maxBounds, ph);
	}
}
