load(paste(parentDir, "screenedData/all_base_obs_filtered.RData", sep=""))

genomeFileSul <- "~/IT/SulfobusTokodaii_bases.gc"
genomeFileBac <- "~/IT/BacillusAmyloliquefaciens_bases.gc"

gcFile <- read.table(genomeFile, sep="\t");
names(gcFile) <- c("RefRealPos", "Base")
genomeLength <- max(gcFile$RefRealPos) + 1 #take care of RHS

rep.mask1 <- unlist(apply(repeatLocations ,1,function(x) seq(x[1],x[1] + x[3] - 1)))
rep.mask2 <- unlist(apply(repeatLocations ,1,function(x) seq(x[2],x[2] + x[3] - 1)))
rep.mask <- c(rep.mask1, rep.mask2)

## is there a better way to do this? I take the rep mask, and assign it to a genomic bin

binsWhichRepeatsFallInto <- cut(x=rep.mask, breaks=seq(from=0, to=genomeLength, by=windowSize), include.lowest = TRUE, include.highest=FALSE)

#I work out the levels of the genomic bins
myLevels <- levels(binsWhichRepeatsFallInto)
df.myLevels <- data.frame(myLevels)

#then I count 
countsPerLevel <- apply(df.myLevels, MARGIN=1, FUN=function(x){length(which(binsWhichRepeatsFallInto == x))}) 


bin.tab = table(binsWhichRepeatsFallInto)
toUse = factor(names(bin.tab[bin.tab == 0]), levels=myLevels)

library(data.table)
#what is dt_full
setkey(dt_full, CallLen)
dt_full = dt_full[CallLen > 0,]

dt_full = dt_full[which(!(is.na(dt_full$FlowVal))),]

posStrand <- which(dt_full$Strand == 1) #is this necessary
negStrand <- which(dt_full$Strand == -1) #is this necessary

setkey(dt_full, Strand)


flPos <- as.vector(dt_full[list(1), c("CallLen"), with=FALSE]) # needs to be a vector.
nucPos <- with(dt_full[list(1),], rep(RefBase, CallLen))
refPos <- with(dt_full[list(1),], rep(RefRealPos, CallLen) + unlist(sapply(CallLen - 1, function(x) seq.int(0,x))))

#Why does this have an invalid times argument!

flNeg <- dt_full[list(-1), c("CallLen"), with=FALSE]
nucNeg <- with(dt_full[list(-1),], rep(RefBase, CallLen))
refNeg <- with(dt_full[list(-1),], rep(RefRealPos, CallLen) - unlist(sapply(CallLen - 1, function(x) seq.int(0,x))))

## finally move on to the next bits.
combnuc <- c(as.character(nucPos), as.character(nucNeg))
combref <- c(refPos, refNeg)

expanded2 = data.table("RefRealPos" = combref, "Nucleotide"=combnuc)


#expanded2 <- data.frame("RefRealPos" = combref, "NUCLEOTIDE"=combnuc)

coverage = expanded2[, list("num_reads_covering" = length(Nucleotide), "nucleotide" = unique(Nucleotide)), by=list(RefRealPos)]

gcBins <-  cut(x=gcFile$RefRealPos, breaks=seq(from=0, to=genomeLength, by=windowSize), include.lowest = TRUE, include.highest=FALSE)
coverageBins <- cut(x=coverage$RefRealPos, breaks=seq(from=0, to=genomeLength, by=windowSize), include.lowest = TRUE, include.highest=FALSE)
numBins <- floor(genomeLength / windowSize)

gcContent <- function(x)
{
        num <-  sum( as.numeric(x == "G" | x == "C" | x == "c" | x == "g"))
        total <- length(x)
        prop <- num / total
}

gcForBin <- tapply(gcFile$Base, gcBins, gcContent)

#instead of being the names, just take the middle position of the bin

myNames <- (windowSize / 2) + windowSize * (0 : (numBins - 1))

namesForBin <-  myNames

coverageForBin <- tapply(coverage$num_reads_covering, coverageBins, mean)

namesForBin <- namesForBin[toUse]
coverageForBin <- coverageForBin[toUse]
gcForBin <- gcForBin[toUse]

gcBreaks <- cut(x=gcForBin, breaks=seq(from=0, to=1.01, by=0.10),include.lowest = TRUE, include.highest=FALSE)

coverageForBin[which(is.na(coverageForBin))] <- 0

bigDF <- data.frame(binName = namesForBin, coverage = coverageForBin, gcForBin = gcBreaks)

write.table(bigDF, file=paste(parentDir,"dataframe_gc_and_breaks_", windowSize, ".csv", sep=""), sep=",")

