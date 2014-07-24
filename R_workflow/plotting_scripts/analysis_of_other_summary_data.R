#analysis of additional files (substitution rates and indel rates mainly.

#what files do I have?

label100 <- "Ion OneTouch \n Template Kit"
label200old <- "Ion Xpress \n Template 200 kit"
label200new <- "Ion OneTouch 200 \n Template kit"



library(ggplot2)
theme_set(theme_bw())

data <- read.table("consolidated_error_counts.csv", sep=",", header=TRUE)


tabCounts <- data[, c("run", "kit", "species", "chip", "correctHomopolymerMinusPoly", "subSNPCount", "subSNPLocs")]

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(2,3,4,8)]

data$group <- ""
data$group[which(data$kit == "200")] <- label200old
data$group[which(data$kit == "100")] <- label100
data$group[which(data$kit == "200Plus")] <- label200new
data$group <- factor(data$group,levels=c(label100, label200old, label200new))


data$prop51 <- data$num51SNPerr / data$substitutionErr
summary(data$prop51)

#> summary(data$prop51)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1635  0.2093  0.2495  0.2506  0.2877  0.3306 



names(data)[4] <- c("chip")
names(data)[5] <- c("CoverageHPMinusPoly")


#insertion rate

data$delRate <- with(data, deletionErr / (CoverageHPMinusPoly + substitutionErr))
data$insRate <- with(data, insertionErr / (CoverageHPMinusPoly + substitutionErr))

#to make this work in the boxplot, need to expand crap

dataSubs <- data[, c("run", "kit", "species", "group", "substitutionErrRate")]
dataSubs$errType <- "Substitutions"
names(dataSubs)[5] <- "rate" 
dataIns <- data[, c("run", "kit", "species", "group", "insRate")]
dataIns$errType <- "Insertions"

names(dataIns)[5] <- "rate"

dataDels <- data[,c("run", "kit", "species", "group", "insRate")]
dataDels$errType <- "Deletions"

names(dataDels)[5] <- "rate"


combData <- rbind(dataSubs, dataIns, dataDels)

library(scales)

m <- ggplot(combData, aes(x = errType, y=rate, fill=errType)) + facet_wrap(~ group,nrow = 1,scale = "free_y")   + 
geom_boxplot() + coord_cartesian(ylim=c(0.0001,0.05)) + scale_fill_manual(breaks="errType", values=cbPalette)  + scale_y_log10(labels=percent_format()) + 
xlab("Error type") + ylab("Error rate")+ opts(legend.position = "none") + 
opts(axis.text.x=theme_text(angle=-90))  + opts(strip.text.x = theme_text(size = 11, colour = "black"))



postscript("avg_rates_of_error_across_datasets.eps",horiz=TRUE,onefile=FALSE,width=8.5,height=11,paper="A4")
m
dev.off()


#done

data <- read.table("consolidated_sub_locations.csv", sep=",", header=TRUE)

#would appear that the majority are a single case in each dataset
errByRefLoc <- with(data, tapply(run, referenceLocation, function(x){length(unique(x))}))

#majority are 1-off... possible that some still coincide with repetitive regions.


data <- read.table("consolidated_snp_rate_by_base_pos.csv", sep=",", header=TRUE)

data$group <- ""
data$group[which(data$kit == "200")] <- label200old
data$group[which(data$kit == "100")] <- label100
data$group[which(data$kit == "200Plus")] <- label200new
data$group <- factor(data$group,levels=c(label100, label200old, label200new))
data$comboStr <- paste(data$species, ".", data$chip, sep="")



m <- ggplot(data, aes(x = BasePos, y=SNPRate, colour=comboStr, group=run)) + facet_wrap(~ group,nrow = 1,scale = "free_y")+ 
geom_line(data = subset(data, group=label100)) +  
geom_line(data = subset(data, group=label200old)) +
geom_line(data = subset(data, group=label200new)) + coord_cartesian(ylim=c(0,0.01), xlim=c(0,100))+
xlab("Base position") + ylab("Error rate") + opts(title="Substitution rate by base position(bases 1-100)") +
scale_colour_manual(values=cbPalette) + opts(strip.text.x = theme_text(size = 11, colour = "black")) 


n <- ggplot(data, aes(x = BasePos, y=SNPRate, colour=comboStr, group=run)) + facet_wrap(~ group,nrow = 1,scale = "free_y")+ 
geom_line(data = subset(data, group=label100)) +  
geom_line(data = subset(data, group=label200old)) +
geom_line(data = subset(data, group=label200new)) + coord_cartesian(ylim=c(0,0.05), xlim=c(0,200))+
xlab("Base position") + ylab("Error rate")+ opts(title="Substitution rate by base position(bases 1-200)")+
scale_colour_manual(values=cbPalette) + opts(strip.text.x = theme_text(size = 11, colour = "black"))




library(gridExtra)
postscript("substitution_rate_by_base_position.eps",horiz=TRUE,onefile=FALSE,width=8.5,height=11,paper="A4")
grid.arrange(m, n, ncol=1)
dev.off()

###
data <- read.table("consolidated_substitution_rates_between_bases.csv",sep=",", header=TRUE)
data$group <- ""
data$group[which(data$kit == "200")] <- label200old
data$group[which(data$kit == "100")] <- label100
data$group[which(data$kit == "200Plus")] <- label200new
data$group <- factor(data$group,levels=c(label100, label200old, label200new))

data$comboStr <- paste(data$species, ".", data$chip, sep="")
data$subType <- paste(data$RefBase,"->", data$ReadBase, sep="")

plot1 <- ggplot(data, aes(x = subType, y=Freq, fill=RefBase)) + facet_wrap(~ group,nrow = 1,scale = "free_y") + 
geom_boxplot(data = subset(data, group=label100)) +  
geom_boxplot(data = subset(data, group=label200old)) +
geom_boxplot(data = subset(data, group=label200new)) + opts(axis.text.x=theme_text(angle=-45)) + 
xlab("Substitution type") + ylab("Proportion of substitutions") + scale_fill_manual(values=cbPalette) + opts(strip.text.x = theme_text(size = 11, colour = "black"))

postscript("substitution_preferences.eps",horiz=TRUE,onefile=FALSE,width=8.5,height=11,paper="A4")
plot1
dev.off()

#despite all these plots, how about I try and visualise read length distributions per datasets, the mean may not be that informative

readLen <- read.table("consolidated_read_lengths.csv", sep=",")
names(readLen) <- c("run", "kit", "species", "chip", "ReadLength", "Count", "TotalReads")

readLen$group <- ""
readLen$group[which(readLen$kit == "200")] <- label200old
readLen$group[which(readLen$kit == "100")] <- label100
readLen$group[which(readLen$kit == "200Plus")] <- label200new
readLen$group <- factor(readLen$group,levels=c(label100, label200old, label200new))

readLen$comboStr <- paste(readLen$species, ".", readLen$chip, sep="")
readLen$freq <- readLen$Count / readLen$TotalReads

readPlot <- ggplot(readLen, aes(x = ReadLength, y=freq, colour=comboStr,group=run)) + facet_wrap(~ group,nrow = 1,scale = "free_y")+ 
geom_line(data = subset(readLen, group="100")) +  
geom_line(data = subset(readLen, group="200-old")) +
geom_line(data = subset(readLen, group="200-onetouch"))  + coord_cartesian(ylim=c(0,0.05), xlim=c(0,800)) + 
ylab("Frequency") + xlab("Read length")+ scale_colour_manual(values=cbPalette) + opts(strip.text.x = theme_text(size = 11, colour = "black"))


postscript("read_length_distribution.eps",horiz=TRUE,onefile=FALSE,width=8.5,height=11,paper="A4")
readPlot
dev.off()

#ignore for now gag errors
data <- read.table("high_occ_indels.csv", sep=",")
names(data) <- c("run", "kit", "species", "chip", "read", "FlowVal", "FlowPos", "RlePos", "ReadPos", "RefRLEPos", "RefRealPos", "RefLen", "RefBase", "Strand")

data$phase <- data$FlowPos %% 32
data$cycle <- data$FlowPos %/% 32
data$kitspecies <- paste(data$kit, data$species, sep=".")
data$magnitude <- data$RefLen - round(data$FlowVal)
data$type <- ifelse(data$magnitude > 0, "Deletion", "Insertion") 
data$RefBase[which(data$RefBase == 't')] <- 'T'
data <- data[!is.na(data$phase),]


uniqType <- unique(data[,c("run", "RefRealPos", "RefLen")])

#deletions appear to be the most common...





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


if(0)
{

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


#for luca, generates the reference?

myRes <- ddply(megaStructure, .(RefRealPos), summarise,
mode.phase = modePhase(FCYCNUM,Strand),
mode.count = modeCount(FCYCNUM,Strand),
mode.strand = domStrand(Strand),
strand.dom.freq = strandDomination(Strand),
strand.count = strandObs(Strand))

write.table(myRes, file="bias_in_phase_in_normal_data.csv", sep=",")
)

#okay what would I like to plot if I could plot something.

#save(myRes, file="dominance_of_strand_error.RData")
load("dominance_of_strand_error.RData")
load("combined_normal_data_dom_strand.RData")

library(ggplot2)

library(vcd)
pal <- heat_hcl(7)

myRes$niceStrand <- ""
myRes$niceStrand[which(myRes$mode.strand == "-1")] <- "Antisense"
myRes$niceStrand[which(myRes$mode.strand == "1")] <- "Sense"
myRes$niceStrand <- factor(myRes$niceStrand, levels=c("Sense", "Antisense"))



res$niceStrand <- ""
res$niceStrand[which(res$mode.strand == "-1")] <- "Antisense"
res$niceStrand[which(res$mode.strand == "1")] <- "Sense"

res$niceStrand <- factor(res$niceStrand, levels=c("Sense", "Antisense"))

res$category <- "No context-specific errors"
myRes$category <- "Context-specific errors"

subsetRes <- res[sample(1:630000, 43753),]

combData <- rbind(myRes[, intersect(names(myRes), names(res))], subsetRes[, intersect(names(myRes), names(res))])

combData$category <- factor(combData$category, levels=c("No context-specific errors", "Context-specific errors"))


plot1 <- ggplot(combData, aes(x = strand.dom.freq, y =(mode.count / strand.count)), main="Phase bias in common errors") + facet_wrap(~ niceStrand + category) +
stat_density2d(geom="tile", aes(fill=..density..), contour=FALSE) +
scale_fill_gradient(low="red", high="yellow",breaks=c(0,2, 10,15,20)) + 
xlab("Dominance of strand") +  ylab("Dominance of phase (in reads on dominant strand)") + 
coord_cartesian(xlim = c(0.5, 1), ylim=c(0,1)) + theme_bw() + opts(axis.text.x=theme_text(angle=-90))  + opts(strip.text.x = theme_text(size = 11, colour = "black")) +
 scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + opts(panel.margin = unit(0.15, "inches"))

svg("context_specific_error_strand_and_phase_dominance_versus_normal.svg", height=8, width=10)
plot1
dev.off()



ofInt <- myRes$RefRealPos[which(myRes$strand.dom.freq > 0.80 & (myRes$mode.count / myRes$strand.count) > 0.5)]
mySub <- myRes[which(myRes$RefRealPos %in% ofInt),]


plot2 <- ggplot(subsetRes, aes(x = strand.dom.freq, y =(mode.count / strand.count)), main="General phase bias") + facet_wrap(~ niceStrand) +
stat_density2d(geom="tile", aes(fill=..density..), contour=F, n=c(30, 30)) +
xlab("Dominance of strand") + scale_fill_gradient(low="red", high="yellow",breaks=c(0,2,5,8,10,15)) +  ylab("Dominance of phase (in reads on dominant strand)") + coord_cartesian(xlim = c(0.5, 1), ylim=c(0,1)) 

#removed breaks ,) 

library(gridExtra)

postscript("comparison_of_phase_of_strand_biased_errors_versus_regular_data.eps")
grid.arrange(plot1, plot2, ncol=1)
dev.off()

s <- unique(data[,c("run", "species", "RefRealPos", "RefBase", "Strand")])

kitErrorOcc <- table(s$RefBase , s$run)

csum <- apply(kitErrorOcc, MARGIN=2, sum)
kitErrorOcc[,1:15] / csum[1:15]

