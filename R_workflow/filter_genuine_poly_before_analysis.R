library(plyr)

rep.mask1 <- unlist(apply(repeatLocations ,1,function(x) seq(x[1],x[1] + x[3] - 1)))
rep.mask2 <- unlist(apply(repeatLocations ,1,function(x) seq(x[2],x[2] + x[3] - 1)))
rep.mask <- c(rep.mask1, rep.mask2)

### data.tables
library(data.table)
tmp_dt <- fread(paste(parentDir, "rawData/flowvaluestolengths_for_overlapping_rle.txt", sep="")) ## read in all the data.
setnames(tmp_dt, c("ReadID", "FlowVal", "FlowPos", "RLEPos", "ReadPos", "RefRLEPos", "RefRealPos", "RefLen", "RefBase", "Strand", "ValidFV", "OOP", "CallLen"))

### data.tables
sbd_dt <- fread(paste(parentDir, "rawData/deletions_all.txt", sep=""), na.strings=c("N\\A", "N/A"))
sbd_dt$CallLen = 0
setnames(sbd_dt, c("ReadID", "FlowVal", "FlowPos", "RLEPos", "ReadPos", "RefRLEPos", "RefRealPos", "RefLen", "RefBase", "Strand", "ValidFV", "OOP", "CallLen"))
setkey(sbd_dt, ReadPos)
sbd_dt  = sbd_dt[!is.na(ReadPos)]
tmp_dt <- rbindlist(list(tmp_dt, sbd_dt))

rm(sbd_dt)

##can't remove NAs, but could do it in the rbindlist above.

ins_dt <- fread(paste(parentDir,"rawData/flowvalues_insertions_full.txt", sep=""))

#insertions <- read.table(paste(parentDir,"rawData/flowvalues_insertions_full.txt", sep=""), sep="\t")
setnames(ins_dt, c("ReadID","FlowVal","FlowPos","RLEPos", "ReadPos", "RefBaseBefore", "RefBaseAfter","INS_NUMBER","RefBase","Strand", "ValidFV", "OOP", "CallLen")) #this will have a call length
## alot of rearranging happens after this.
## instead of doing filter separately...
setkey(ins_dt, Strand) #make it on strand...
ins_dt[, RefRealPos :=  ifelse(Strand == 1, RefBaseBefore, RefBaseAfter)]  ###works
ins_dt[,RefRLEPos:=RefRealPos]
ins_dt[,RefLen:= 0]

#get RID of columns we don't care about
ins_dt[, c("RefBaseBefore", "RefBaseAfter", "INS_NUMBER") := NULL]
tmp_dt <- rbind(tmp_dt, ins_dt) 
rm(ins_dt)
gc()


### Manipulating the data in preparation.
tmp_dt[,ReadPos := ReadPos + 1]



# only keep rows where these are true?
# maybe set it up as a big function.
tmp_dt = tmp_dt[ReadPos >= startPos & tmp_dt$ReadPos <= truncatePos & ! is.na(FlowVal),] ###this works
setkey(tmp_dt, RefRealPos)
tmp_dt= tmp_dt[!list(rep.mask),] ### this did not work.
setkey(tmp_dt, RefRealPos)
###tmp_dt_2[J(422729),] ## this is to test whether the rep mask worked.
gc()



#tmp_dt = tmp_dt[!(tmp_dt$RefRealPos %in% rep.mask),]
#tmp_dt <- tmp_dt[tmp_dt$ReadPos <= truncatePos,]
#tmp_dt <- tmp_dt[!(tmp_dt$RefRealPos %in% rep.mask),]
#nas <- with(tmp_dt, which(is.na(FlowVal)))

##V. Old bases on reliable flow values, no longer usable.
#tmp$Error <- with(tmp, (is.na(FlowVal) | round(FlowVal,0)!=RefLen))
#tmp$Errors <- with(tmp, ifelse(is.na(FlowVal), -RefLen, round(FlowVal,0)-RefLen))

### Data.table version of this stuff
tmp_dt[,Error := CallLen != RefLen]

tmp_dt[, Errors :=  ifelse(is.na(CallLen), -RefLen, CallLen - RefLen)]
tmp_dt[, RefBase := toupper(RefBase)]
tmp_dt[, baseOnPosStrand := RefBase] #this is not great, but not much help about
setkey(tmp_dt, Strand, RefBase)
tmp_dt[J(-1, "A"), baseOnPosStrand:= "T"] #not sure if this is correct.
tmp_dt[J(-1, "T"), baseOnPosStrand:= "A"]
tmp_dt[J(-1, "C"), baseOnPosStrand:= "G"]
tmp_dt[J(-1, "G"), baseOnPosStrand:= "C"]

## Error not based on the FlowValue, but on the Called Len (i.e. What Ion Torrent called it, with phase correction)


## I think it might be safer to consider the data separately for those with flow-values and those without... or we can do things not by flow, later on.
##so what I was trying to do is work out the error-rate by flow position?

## Big Decisions to be made here! How to deal with OOP reads. What about reads that are short? 
## A comment to make in the paper is how different reads are tmp.fp versus tmp.no_fp

## There are a group short reads that are artificially upping error-rate at the 5' end. 
## There are reads with no flow-values... 
## The flow values are not so important, if only for estimating the number of correct zeroes.

## this many rows currently: 175,854,125

read_aln_dt = tmp_dt[, list(min_aln = min(ReadPos), max_aln = max(ReadPos), errors=sum(Error)), by=ReadID]
read_aln_dt[, aln_len := max_aln - min_aln + 1]#

read_aln_list = read_aln_dt[aln_len >= 30, ReadID]  ##IDs to remove.
setkey(tmp_dt, ReadID)
tmp_dt = tmp_dt[J(read_aln_list)]  ## this worked for some reason.

rm(read_aln_dt)
rm(read_aln_list)

unaligned_length = fread(paste(parentDir, "/rawData/read_lengths.tsv", sep=""))
setnames(unaligned_length, c("ReadID", "MaxReadPos", "MaxRLEPos"))

read_unalign_list = unaligned_length[MaxReadPos >= 30 & MaxReadPos <= 500, ReadID]

#setkey(tmp_dt, ReadID)

tmp_dt = subset(tmp_dt, ReadID %in% read_unalign_list) 

## We removed reads that were shorter than 30bp or longer than 500bp.
## We removed alignments that were shorter than 30bp.

rm(unaligned_length)
rm(read_unalign_list)


tmp_dt = tmp_dt[FlowPos >= 10] ##check this works

#So we've removed the non-representative reads.

## maybe just subset out the middle 

####just use the last portion for now. Must remember that memory usage gets too high if I use 800,000 reads.

#tmp_dt_sub = tmp_dt[76122248:176122248,]
#tmp_dt = tmp_dt_sub



## consider an error-prevalence by read length?
## In words, I am looking at the coverage of sites where the RefLen > 0, I want the data.table to consist of RefRealPos, num of flows aligning to that position,  the ref len for that homopolymer
## 

my_total_cover_dt = tmp_dt[RefLen > 0, list(num_flows=length(FlowPos), hp_len=unique(RefLen)), by=RefRealPos] ###this works.
adj_bases = function(r){return (list(r -1, r + 1))}

my_total_cover_dt[, c("baseBefore","baseAfter") := adj_bases(RefRealPos)] ###works

prob_by_hp_dt = tmp_dt[RefLen > 0, list(p = mean(Error)), by = RefLen] ## Works.



## To come up with a clever way to do this using data.table before I go today.

## what I want is the total number of observations at each flow
## and the total number of non-observations at each flow

ids = unique(tmp_dt$ReadID)   # or from the data: unique(ds$ID)
pos = unique(tmp_dt$FlowPos[which(tmp_dt$FlowPos >= 0)])   # or from the data: unique(ds$Pos) 
### end point will be the last FlowPos before end or OOP.

setkey(tmp_dt, ReadID,FlowPos)

##all possible flows (except for those that are negative flow-call)
calls_at_flows = tmp_dt[CJ(ids,pos), roll=-Inf, nomatch=0][, .N, by=FlowPos] #calls across all flow positions
setkey(calls_at_flows, FlowPos)
## above is number of reads at least as long as above 

setkey(tmp_dt, RefLen)

### actually called (observed) flows
called_non_zero = tmp_dt[RefLen > 0 & FlowPos > 0, list(nonzero_calls = length(ReadID)), FlowPos]
setkey(called_non_zero, FlowPos)

##Only focusing on those that have a positive FlowPos.
error_by_flow_zero = tmp_dt[J(0),][FlowPos > 0, list(num_zero_uc = sum(Error)), by=FlowPos] 
setkey(error_by_flow_zero, FlowPos)

## join all the tables by FlowPos key.
t = error_by_flow_zero[calls_at_flows[called_non_zero]] ## so we have a table for each flow position, with zero and non zero counts.

## Calculate the number of flow positions aligned with a zero RefLen (base does not exist in the reference sequence)
t[, call_w_zero_reflen := N - nonzero_calls]


## replace NAs
t[is.na(num_zero_uc), num_zero_uc:= as.integer(0)] 

### finally can calculate the probability of seeing a zero RefLen overcalled.
p = with(t,  sum(num_zero_uc)/ sum(call_w_zero_reflen))


### Location regions of the genome which seem to have an undercall of 1, and then work out if bona fide or not.

setkey(tmp_dt, RefRealPos, baseOnPosStrand)

tEl.dt = tmp_dt[FlowPos > 0 & RefLen == 0, list(errorCount = sum(Error)), by=list(RefRealPos, baseOnPosStrand)]  ## good

#tEL.df <- aggregate(Error ~ RefRealPos + baseOnPosStrand, data=singles[which(singles$FlowPos > 0),],FUN=sum) 

setnames(tEl.dt, c("InsPos", "InsertedBase", "errorCount"))

### special merge, need to link together not identical names.


###Tandom errors exist in the data.

mergedLeft <- merge.data.frame(tEl.dt, my_total_cover_dt, by.x="InsPos", by.y="RefRealPos", all.X=TRUE)[,c(1:4)]   #[,c(1:4)]

#mergedLeft <- merge(tEL.df, myTotalCover.df, by.x="InsPos", by.y="RefRealPos", all.X=TRUE)[,c(1:4)]
mergedFull <- merge.data.frame(mergedLeft, my_total_cover_dt, by.x="InsPos", by.y="baseBefore", all.X=TRUE)[,c(1:4,6)]


names(mergedFull)[4] <- "coverageBefore"
names(mergedFull)[5] <- "coverageAfter"
	
mergedFull$coverageBefore <- as.vector(mergedFull$coverageBefore)
mergedFull$coverageAfter <- as.vector(mergedFull$coverageAfter)

mergedFull$coverageBefore[which(is.na(mergedFull$coverageBefore))] <- 0
mergedFull$coverageAfter[which(is.na(mergedFull$coverageAfter))] <- 0

maxTotalCover <- apply(mergedFull[,c("coverageBefore", "coverageAfter")], MARGIN=1, FUN=max)
mergedFull$maxCover <- maxTotalCover

#what is the expected error rate for over-calling a zero
p.value <- with(mergedFull, 1 - pbinom(errorCount - 1, size=maxCover, prob=pZero))

mergedFull$p.value <- p.value

resZero <- mergedFull[,c("InsPos","p.value")] 
names(resZero) <- c("location", "p.value")

## Merged full is interesting, as it shows blocks of inserts..

write.table(resZero, file=paste(parentDir, "possible_indel_errors_in_ref_oc_zero.txt", sep=""), sep=",")


### What do I want this data.table to look like?
### Apparently I want the data frame to have location count coverage RefLen probability for the RefLen.

setkey(tmp_dt, RefRealPos, RefLen)

hpErrByLoc_dt = tmp_dt[RefLen > 0, list(count=sum(Error)), by=list(RefRealPos, RefLen)] ##Gives me RefRealpos, RefLen Count
setkey(hpErrByLoc_dt, RefRealPos)
setkey(my_total_cover_dt, RefRealPos)

mm = my_total_cover_dt[hpErrByLoc_dt]
mm[,baseAfter := NULL]
mm[,baseBefore := NULL]
mm[,hp_len := NULL]

## Above is clean enough now.
## Have position, coverage, RefLen, count (errors)

## Need to get the correct p-value

new_p <- function(coverage,count, p)
{
  return(1 - pbinom(count - 1, size=coverage, prob=p))
}

## chain together MM and it's probability
setkey(mm, RefLen)
setkey(prob_by_hp_dt, RefLen)

inter = mm[prob_by_hp_dt]
inter[, pVal := new_p(num_flows, count, p)]
inter$adj_pval = p.adjust(inter$pVal, "holm") < 0.05

toIgnore <- inter$RefRealPos[which(inter$adj_pval)] ## working apparently.

#remove these IDS from the main dataset, and also store these weird errors elsewhere
setkey(tmp_dt, RefRealPos)

weird =  subset(tmp_dt, RefRealPos %in% toIgnore)  #tmp_dt[J(toIgnore),]  ## this is the same.
tmp_dt = subset(tmp_dt, ! RefRealPos %in% toIgnore) # tmp_dt[!toIgnore] ## extra rows???



write.table(weird, file=paste(parentDir, "possible_indel_errors_in_ref.txt", sep=""), sep=",")
print(paste("Num locations masked: ", length(unique(weird$RefRealPos)))) 


tmp_dt[FlowPos > 0, FCYCLE := factor(FlowPos %% 32, levels=0:31)]
tmp_dt[FlowPos > 0, FCYCNUM := FlowPos %/% 32]

## May be unnecessary clutter
tmp_dt[, XPOS := as.numeric(sapply(strsplit(ReadID, ":"), "[", 2))]
tmp_dt[, YPOS := as.numeric(sapply(strsplit(ReadID, ":"), "[", 3))]
tmp_dt[, NUCLEOTIDE := factor(toupper(RefBase))]
tmp_dt[, GCvAT := factor(ifelse(NUCLEOTIDE %in% c("G", "C"), "GC", "AT"))]
tmp_dt[, baseOnPosStrand := NULL]


## screened data produced.

save(tmp_dt, file=paste(parentDir, "screenedData/all_base_obs_filtered.RData", sep=""))

##3 this is confusing. How can PIC 10 have a low error rate, and a high one?
## Flow positions 9,10,11,12 have higher error rates than anticipated. 
## 

tmp_dt = tmp_dt[FlowPos >= 10]
setkey(tmp_dt, FlowPos)
flowpos_tab = tmp_dt[FlowPos > 0, list(num_obs=length(ReadID)), by=FlowPos] ## some sites called more frequently than others.

flowpos_error_tab = tmp_dt[FlowPos > 0, list(error=mean(Error)), by= list(FlowPos)] ##Valid FV is no use.

readpos_error_tab = tmp_dt[, list(error=mean(Error)), by=ReadPos]

### Let's make a new column which gives 'Error' based on the Flow Val.
tmp_dt[! is.na(FlowVal), ebpc := round(FlowVal,0) != RefLen]

## for each read, compute the error rate
phase_corr_versus_raw = tmp_dt[! is.na(FlowVal), list(rawErr = sum(ebpc), phaseCorrErr = sum(Error)), by=ReadID] ## for each read, work out the number of errors before and after correction?

## generally better after phase correction. Some reads are not better after phase correction.
## Summarise the before and after corrections
write.csv(phase_corr_versus_raw[, list(leRaw=sum(rawErr < phaseCorrErr), eqRaw=sum(rawErr == phaseCorrErr), gtRaw=sum(rawErr > phaseCorrErr), p.val=t.test(rawErr, phaseCorrErr, paired=TRUE, alternative="greater")$p.val)], file=paste(parentDir,"errors_after_phase_corr.csv", sep=""))

### final outputs.
write.table(flowpos_error_tab, file=paste(parentDir, "error_rate_by_flow_position.csv", sep=""), sep=",")
write.table(readpos_error_tab, file=paste(parentDir, "error_rate_by_base_position.csv", sep=""), sep=",")
write.table(tmp_dt[,list(ReadID, CallLen, ReadPos, Errors)], file=paste(parentDir, "flow_data.csv", sep=""), sep=",")
