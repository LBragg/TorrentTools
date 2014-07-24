
setwd("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/combined_results_inhouse_400bp")

#ignore for now gag errors
dt <- fread("high_occ_indels.csv", sep=",", header=TRUE)
setnames(dt,c("chip", "desc", "kit", "species", "read", "FlowVal", "FlowPos", "RlePos", "ReadPos", "RefRLEPos", "RefRealPos", "RefLen", "RefBase", "Strand", "Valid_FV", "OOP", "CALL_LEN", "BOPS"))
#names(initdata) <- c("run", "kit", "species", "chip", "read", "FlowVal", "FlowPos", "RlePos", "ReadPos", "RefRLEPos", "RefRealPos", "RefLen", "RefBase", "Strand")

dt_old = fread("/home/bra427/Projects/IonTorrentBenchmarking/resultsIncl400bp/high_occ_indels.csv", sep=",", header=TRUE)




dt[, Error := (CALL_LEN != RefLen)]
dt[, magnitude := (RefLen - CALL_LEN)]
dt[, type := ifelse(magnitude > 0, "Deletion", ifelse(magnitude < 0, "Insertion", "Correct"))]

#for this table, you need *distinct* stuff
### What is this doing.


##



domType <- function(x)
{
  x <- x[which(x != "Correct")]
  tabType <- table(x)
  return(names(tabType[which(tabType == max(tabType))]))
}

tab1 <- with(dt, tapply(RefRealPos, list(desc, type), function(x) length(unique(x))))
test <- ddply(dt, .(desc, RefRealPos), summarize, domType=domType(type))

domPerRun <- with(test, tapply(RefRealPos, list(desc, domType), length))

write.csv(domPerRun, file="error_loc_per_run_focus_on_dominant_error_type.csv")




tab1 <- with(dt, tapply(type, list(desc, RefRealPos), domType) 
#what you would want is the dominant error at a loci?

write.csv(tab1, file="general_poly_counts.csv")

#come back to init data later.

data <- data[,c("run", "Strand", "RefRealPos", "Error") ]

#data$RefBase[which(data$RefBase == 't')] <- 'T'
#data$phase <- data$FlowPos %% 32
#data$cycle <- data$FlowPos %/% 32
#data$kitspecies <- paste(data$kit, data$species, sep=".")
#
#why this ??
data <- data[!is.na(data$phase),]


write.table(dt$read, file="all_hfi_identifiers.csv", sep=",")



#dataSub <- data
#library(plyr)

myNew <- ddply(dt, .(desc, Strand, RefRealPos), summarize, strand.count=length(Error), error.count=sum(Error))

myNew.pos <- myNew[which(myNew$Strand == 1),]
myNew.neg <- myNew[which(myNew$Strand == "-1"),]

my.merged <- merge(myNew.pos, myNew.neg, by.x=c("desc", "RefRealPos"), by.y=c("desc", "RefRealPos"), all.x=TRUE, all.y=TRUE)

ind <- is.na(my.merged)
my.merged[ind] <- rep(0, sum(ind))

my.merged <- my.merged[,c(1:2,4:5,7:8)]
names(my.merged) <- c("run", "RefRealPos", "ObsPos", "ErrPos", "ObsNeg", "ErrNeg")

my.merged$domStrand <- with(my.merged, ifelse(ObsPos >= ObsNeg, "1", "-1"))
my.merged$propDom <- with(my.merged, ifelse(ObsPos >= ObsNeg, (ObsPos / (ObsPos + ObsNeg)), 
(ObsNeg / (ObsPos + ObsNeg))))

my.merged$errorsDom <- with(my.merged, ifelse(ObsPos >= ObsNeg, ErrPos, ErrNeg))
my.merged$totalError <- with(my.merged, ErrPos + ErrNeg)

my.merged$p <- 1 - pbinom(q=my.merged$errorsDom - 1, p=my.merged$propDom, size=my.merged$totalError)
my.merged$wh <- p.adjust(my.merged$p, "holm")  < 0.05

### 2432 instances, only 11 detected as significant strand bias.. however many more examples which did not have enough coverage on both strands.



#okay, to do Glenns test, need to aggregate by run, RefRealPos, Strand
withInPhase = dt[ RefRealPos > 0 , list(strand.count = length(Error), error.count = sum(Error)),by="desc,Strand,RefRealPos,OOP"]

pos_strand_false = withInPhase[Strand == "1" & OOP == "False", list(desc, RefRealPos, strand.count, error.count),]
setkey(pos_strand_false, desc, RefRealPos)
setnames(pos_strand_false, "strand.count", "strand.count.pos.false")
setnames(pos_strand_false, "error.count", "error.count.pos.false")

pos_strand_true = withInPhase[Strand == "1" & OOP == "True", list(desc, RefRealPos, strand.count, error.count),]
setnames(pos_strand_true, "strand.count", "strand.count.pos.true")
setnames(pos_strand_true, "error.count", "error.count.pos.true")
setkey(pos_strand_true, desc, RefRealPos)
neg_strand_false = withInPhase[Strand == "-1" & OOP == "False", list(desc, RefRealPos, strand.count, error.count),]
setnames(neg_strand_false, "strand.count", "strand.count.neg.false")
setnames(neg_strand_false, "error.count", "error.count.neg.false")

setkey(neg_strand_true, desc, RefRealPos)
neg_strand_true = withInPhase[Strand == "-1" & OOP == "True", list(desc, RefRealPos, strand.count, error.count),]
setnames(neg_strand_true, "strand.count", "strand.count.neg.true")
setnames(neg_strand_true, "error.count", "error.count.neg.true")
setkey(neg_strand_false, desc, RefRealPos)

mega_merge = pos_strand_false[pos_strand_true[neg_strand_false[neg_strand_true]]]
mega_merge[is.na(mega_merge)] <- 0
mega_merge[1:10,]

## a quick look, the data looks independent of the phasing of the read... just some evidence for strand specificity.


#maybe should report p-value, and significant sites after exclusion of p-value?


#subset it

my.res <- my.merged[,c("run", "RefRealPos", "p", "wh")]


alldata <- merge(initdata, my.res, by.x=c("run", "RefRealPos"), by.y=c("run", "RefRealPos"), all.x=TRUE)


rowOfI <- alldata[which(alldata$wh),]
with(rowOfI, tapply(RefRealPos, list(run, type), function(x) length(unique(x))))

with(rowOfI, tapply(kit, list(run,type,RefRealPos), length)


#ignoring the error 'rate' we consider whether the proportion of errors observed on one strand over another is greater than expected

test <- ddply(rowOfI, .(run, RefRealPos), summarize, domType=domType(type))

domPerRun <- with(test, tapply(RefRealPos, list(run, domType), length))

write.csv(domPerRun, "strand_specific_error_per_run_by_dom_type.csv")



#dominant phase
domPhase <- function(x)
{
	tabP <- table(x)
	return(names(tabP[which(tabP == max(tabP))]))
}

domFreq <- function(x)
{
	tabP <- table(x)
	dp <- names(tabP[which(tabP == max(tabP))])
	return( max(tabP) / sum(tabP))
}



b <- ddply(rowOfI, .(run, RefRealPos), summarise, domPhase=domPhase(FlowPos %% 32), domFreq=domFreq(FlowPos %% 32))

rowOf










uniqType <- unique(data[,c("run", "RefRealPos", "RefLen")])

#deletions appear to be the most common...

#try strand bias tests for this.





#what I need for table in the paper is 

del <- which(data$type == "Deletion")
ins <- which(data$type == "Insertion")

#if you were going to do it, called tabCounts (for totals)

#join them


delByRun <- with(data[del,], tapply(RefRealPos, run, function(x){length(unique(x))}))
insByRun <- with(data[ins,], tapply(RefRealPos, run, function(x){length(unique(x))}))
byRun <- data.frame(insertions=insByRun, deletions=delByRun, run=names(insByRun))

my.merged <- merge(tabCounts, byRun, by.x="run", by.y="run", all.x=TRUE)

write.csv(my.merged, "counts_of_genomic_locations_for_snps_and_indel_poly.csv")

#okay, interested in seeing if there is a change per kit

data$ReadBase <- data$RefBase
data$ReadBase[which(data$RefBase == "A" & data$Strand == "-1")] <- "T"
data$ReadBase[which(data$RefBase == "T" & data$Strand == "-1")] <- "A"
data$ReadBase[which(data$RefBase == "C" & data$Strand == "-1")] <- "G"
data$ReadBase[which(data$RefBase == "G" & data$Strand == "-1")] <- "C"


#separate apparent 'strand bias' cases

data




#only mode 1
data <- data[which(data$magnitude == 1),]

#common locus?

numRunsWithSameLocus <- with(data, tapply(run, list(species, RefRealPos), function(x){length(unique(x))}))



myData <- aggregate(type ~ ReadBase + kit, data, FUN=length)
myData$group <- ""
myData$group[which(myData$ReadBase == "A" | myData$ReadBase == "T")] <- "AT"
myData$group[which(myData$ReadBase == "C" | myData$ReadBase == "G")] <- "GC"


#now that you've included single-base deletions, are they occurring on one strand only too?

singleBase <- data[which(data$RefLen == 0 & data$run == "31_54"),]

#need to separate by datasets

#need to work out what proportion genuinely are this strand specific issue.

myT <- with(singleBase,tapply(run, list(RefRealPos,Strand), FUN=length))

b <- data.frame(myT)



sumByAT / apply(sumByAT, MARGIN=1, sum)



myData <- aggregate(type ~ kit + ReadBase, data, FUN=length)



#it is the correct counts, minus SNPs
errCounts <- read.table("consolidated_error_counts.csv", sep=",", header=TRUE)
errCounts$totalBases <- errCounts$correctHomopolymerMinusPoly + errCounts$deletionErr + errCounts$insertionErr + errCounts$substitutionErr


#so what is my hypothesis about this.... do they have a particular flow-value distribution?
plot(density(na.omit(data$FlowVal)))
length(which(data$FlowVal == 0.51)) #few have a value of 0.51


length(with(data,which(round(data$FlowVal) == RefLen)))
#they are all errors in this file.
#break down by run

library(ggplot2)
barplot(table(data$kit))


#mean coverage for dataset 42_99 is 17.
#okay lets do just a barplot of the counts per run for these 'errors'

# assign some groups
#for each Ref Pos, by run, work out what the count is for that reference position

covByRun <- function(x)
{
	counts <- table(x)
	locations <- names(counts)
	names(counts) <- c()
	counts <- as.numeric(counts)
	return(cbind(locations, counts))
}

classes <- with(data, tapply(RefRealPos, run, covByRun)) 

res <- c()

for(i in 1:length(classes))
{
	data1 <- classes[i][[1]]
	myds <- names(classes)[i]
	nameVec <- rep(myds, dim(data1)[1])
	block <- cbind(data1, nameVec)
	res <- rbind(res, block)
}

rownames(res) <- c()
res.df <- data.frame(res)


#need to normalize for coverage...
dataMerged <- merge(data, res.df, by.x=c("RefRealPos", "run"), by.y=c("locations", "nameVec"), all.x=TRUE) 
dataMerged$covGroup <- cut(as.numeric(as.character(dataMerged$counts)), breaks=c(0,5,10,15,20,25,30,50)) 

#time to barplot

dataMerged$KitSpecies <- paste(dataMerged$kit, dataMerged$species, sep=".")
merge2 <- merge(dataMerged, errCounts, by.x="run", by.y="run", all.x=TRUE)


library(ggplot2)
myPlot <- ggplot(dataMerged, aes(x=covGroup, y=..density.., fill=kitspecies, group=run)) +
geom_bar(position = "dodge")
myPlot

#kit level errors

table(data$magnitude) #by far the most common magnitude is 1.

distinctRows <- unique(merge2[, c("run", "RefRealPos", "counts", "totalBases")])


errorCounts <- with(distinctRows, tapply(as.numeric(as.character(counts)), run, sum))
counts <- errorCounts
rownames(counts) <- c()
runs <- rownames(errorCounts)

error.df <- data.frame(run=runs, count=counts)

run.df <- merge(error.df, errCounts, by.x="run", by.y="run", all.x=TRUE)
run.df$kitSpecies <- paste(run.df$kit, run.df$species, sep=".")
run.df$kitChip <- paste(run.df$kit, run.df$chip, sep=".")

ggplot(run.df, aes(x=run, y=count, group=kitChip, colour=kitChip)) + geom_point(size=5) + theme_bw()



library(plyr)

processPhase <- function(x)
{
	tabledX <- table(na.omit(x))
	maxCount <- as.numeric(as.character(max(tabledX)))	
	modePhase <- as.numeric(as.character(names(tabledX[which(tabledX == maxCount)])))
	totalObs <- as.numeric(as.character(length(x)))	
	return(data.frame("modeCount"=maxCount, "modePhase"=modePhase, "totalObs"=totalObs))
}


modePhase <- function(x,strand)
{	
#	print("Mode Phase")
#	print(x)
#	print(strand)
	strandOcc <- table(strand)
	domStrand <- names(strandOcc[which(strandOcc == max(strandOcc))][1])
	x <- x[which(strand == domStrand)]
	myT <- table(x)
	retVal = names(myT[which(myT == max(myT))][1])

#	print(paste("Returning", retVal))
	return(retVal)
}

#x <- rep(8,10)
#strand <- rep(1,10)

modeCount <- function(x,strand)
{	
#	print("Mode count")
#	print(x)
#	print(strand)
	strandOcc <- table(strand)
	domStrand <- names(strandOcc[which(strandOcc == max(strandOcc))][1])
	x <- x[which(strand == domStrand)]
	myMax <- max(table(x))[1]

#	print(paste("Returning", myMax))
	return(myMax)
}

domStrand <- function(strand)
{
#	print("Dom Strand")
#	print(strand)
	strandOcc <- table(strand)
	domStrand <- names(strandOcc[which(strandOcc == max(strandOcc))][1])

#	print(paste("Returning ", domStrand))
	return (domStrand)
}


strandDomination <- function(strand)
{
#	print("Dom Strand Prop")
#	print(strand)
	strandOcc <- table(strand)
	domStrandCount <- strandOcc[which(strandOcc == max(strandOcc))][1]
	
#	print(paste("Returning ", domStrandCount, " / ", sum(strandOcc)))

	return(domStrandCount / sum(strandOcc))
}

strandObs <- function(strand)
{
	strandOcc <- table(strand)
	domStrandCount <- strandOcc[which(strandOcc == max(strandOcc))][1]
	return(domStrandCount)
}






myRes <- ddply(data, .(RefRealPos,run), summarise,
mode.phase = modePhase(phase,Strand),
mode.count = modeCount(phase,Strand),
mode.strand = domStrand(Strand),
strand.dom.freq = strandDomination(Strand),
strand.count = strandObs(Strand)
)

#this is saved
#save(myRes, file="aggregated_by_strand_dominance_full_ins_dels.RData")
load("aggregated_by_strand_dominance_full_ins_dels.RData")

#okay I wanted to say definitively how many of these look like the 'real' deal.

library(scatterplot3d)

#try merging myRes with uniqType.
uniqType


my.merged <- merge(myRes, uniqType, by.x=c("RefRealPos", "run"), by.y=c("RefRealPos", "run"), all.x=TRUE)





a200bpMan <- c("43_69", "44_70", "47_75", "33_76", "34_77", "48_78", "35_79")
a200bpOne <- c("41_94",  "42_95", "58_97", "59_99")

#interested in the 200bp kits only

myResMan <- my.merged[which(my.merged$run %in% a200bpMan),]
myResOne <- my.merged[which(my.merged$run %in% a200bpOne),]




numInstancesMan <- myResMan[which(myResMan$strand.count >= 5 & myResMan$strand.dom.freq >= 0.80),]
numInstancesOne <- myResOne[which(myResOne$strand.count >= 5 & myResOne$strand.dom.freq >= 0.80),]

with(numInstancesMan, tapply(RefRealPos, run, length))

> with(numInstancesMan, tapply(RefRealPos, run, length))
30_53 31_54 33_76 34_77 35_79 38_60 40_64 41_94 42_95 43_69 44_70 47_75 48_78 58_97 59_99 
   NA    NA   359   286   969    NA    NA    NA    NA  2079  2462   724   840    NA    NA 

with(numInstancesOne, tapply(RefRealPos, run, length))

> with(numInstancesOne, tapply(RefRealPos, run, length))
30_53 31_54 33_76 34_77 35_79 38_60 40_64 41_94 42_95 43_69 44_70 47_75 48_78 58_97 59_99 
   NA    NA    NA    NA    NA    NA    NA  2597  3862    NA    NA    NA    NA  1050  2207 

length(which(myResMan$strand.count > 5)

smoothScatter(myResMan$strand.count, myResMan$mode.count / myResMan$strand.count)

