rep.mask1 <- unlist(apply(repeatLocations ,1,function(x) seq(x[1],x[1] + x[3] - 1)))
rep.mask2 <- unlist(apply(repeatLocations ,1,function(x) seq(x[2],x[2] + x[3] - 1)))
rep.mask <- c(rep.mask1, rep.mask2)

#complete substitution information needs to go here
#also need complete coverage information, we store that in coverageRLE

#coverage RLE does not exist

#load screened data.
load(paste(parentDir, "screenedData/all_base_obs_filtered.RData", sep=""))

tmp <- tmp[!(is.na(tmp$FlowVal)),]
tmp <- tmp[ tmp$FlowVal >= 0.5,]

posStrand <- which(tmp$Strand == 1)
negStrand <- which(tmp$Strand == -1)

flPos <- round(tmp$FlowVal[posStrand])
nucPos <- with(tmp[posStrand,], rep(NUCLEOTIDE, flPos))
refPos <- with(tmp[posStrand,], rep(RefRealPos,flPos) + unlist(sapply(flPos-1, function(x) seq.int(0,x))))


flNeg <- round(tmp$FlowVal[negStrand])

nucNeg <- with(tmp[negStrand,], rep(NUCLEOTIDE, flNeg))
refNeg <- with(tmp[negStrand,], rep(RefRealPos,flNeg) - unlist(sapply(flNeg-1, function(x) seq.int(0,x))))
combnuc <- c(nucPos, nucNeg)
combref <- c(refPos, refNeg)

expanded2 <- data.frame("RefRealPos" = combref, "NUCLEOTIDE"=combnuc)

myCover <- table(expanded2$RefRealPos)
RefRealPos <- names(myCover)
names(myCover) <- c()
df.cov <- data.frame(RefRealPos, "NUM_READS_COVERING"=myCover)
nuc.df <- unique(expanded2)
my.merged <- merge(df.cov, nuc.df, by.x="RefRealPos", by.y="RefRealPos", all.x=TRUE)

coverage <- my.merged[,c(1,3,4)]
coverage$RefRealPos <- as.numeric(coverage$RefRealPos)
coverage <- coverage[!(coverage$RefRealPos %in% rep.mask),]

#there was some stuff here about GC bias, but removecd it for now.

#just load the substitutions for now, some datasets will have an extra column
data <- read.table(paste(parentDir, "rawData/substitution_information_complete.txt.screened", sep=""), sep="\t", header=TRUE)
names(data) <- c("READ_ID", "RLE_POS", "READ_BASE_POS", "FLOW_POSITION","REF_NUCLEOTIDE", "READ_NUCLEOTIDE", "FLOW_VALUE", "REF_BASE_POS", "STRAND")

#we need to mask reads which fall in repeat locations
#we need to mask the repeat locations themselves.
#we need to mask polymorphisms
#had to increment the read base position for substitutions too.
data$READ_BASE_POS <- data$READ_BASE_POS + 1;

data <- data[data$READ_BASE_POS >= startPos,]
data <- data[data$READ_BASE_POS <= truncatePos,]
data2 <- data[!(data$REF_BASE_POS %in% rep.mask),]

#aggregate the substitution information
polyByReference <- table(data2$REF_BASE_POS)

polyPos <- as.numeric(names(polyByReference))
names(polyByReference) <- c()
polyCount <- polyByReference

df.poly <- data.frame(snpCount=polyCount, position=polyPos)
df.poly2 <- data.frame(snpCount=df.poly$snpCount.Freq, position=df.poly$position)
df.cov <- data.frame(covCount=coverage$NUM_READS_COVERING, position=coverage$RefRealPos, nucleotide=coverage$NUCLEOTIDE)

df.merged <- merge(df.cov, df.poly2,by.x="position", by.y="position", all.x=TRUE, all.y=TRUE)
df.merged$snpCount[is.na(df.merged$snpCount)] <- 0
df.merged$covCount[is.na(df.merged$covCount)] <- 0

df.merged$totalCov <- df.merged$covCount + df.merged$snpCount

#if the 200bp kit, cutoff of 220, 100 bp kit, cutoff of 110.

p <- sum(df.merged$snpCount) / sum(df.merged$totalCov)

#looking for significant substitutions
p.value <- with(df.merged, 1 - pbinom(snpCount - 1, size=totalCov, prob=p))

#pdf("subs_pvalue_hist.pdf")
#hist(p.value, breaks=100)
#dev.off()

#assuming a rate 0.01, only get 11 locations which look like substitutions

table(wh <- (p.adjust(p.value, "holm") <0.05))
polyLoc <- df.merged$position[which(wh == TRUE)]

write.table(polyLoc, file=paste(parentDir, "polymorphism_locations.csv", sep=""), sep=",")

###new
numSNPInstances <- dim(data2[(data2$REF_BASE_POS %in% polyLoc),])[1]
numSNPLocations <- length(polyLoc)

#after removal of polymorphisms

data2 <- data2[!(data2$REF_BASE_POS %in% polyLoc),]
numLocationsSubErrors <- length(unique(data2$REF_BASE_POS))
numSubsErrors <- dim(data2)[1]

errorsPerRead <- table(table(as.character(data2$READ_ID)))

#write the trend in the substitution rate
subOcc <- data.frame(table(data2$READ_BASE_POS))
load(paste(parentDir, "flow_data_expanded.RData", sep=""),header=TRUE) #loads an object called expandedFlow.
totalCov <- nrow(expandedFlow)
covByReadPos <- data.frame(table(expandedFlow$BasePosition))
names(subOcc) <- c("ReadPos", "SNPFreq")
names(covByReadPos) <- c("ReadPos", "HPcov")
df.snpMerged <- merge(y=subOcc, x=covByReadPos, by.x="ReadPos", by.y="ReadPos", all.x=TRUE, all.y=TRUE)
df.snpMerged$rate <- df.snpMerged$SNPFreq / (df.snpMerged$HPcov + df.snpMerged$SNPFreq)
write.table(df.snpMerged, file=paste(parentDir, "SNP_error_rate_by_read_pos.csv", sep=""), sep=",")
rate <- numSubsErrors / (totalCov + numSubsErrors)
subRes <- c()
subRes$numSubs <- numSubsErrors
subRes$covHPNoSNP <- totalCov
subRes$subRate <- rate
subRes$snpLocs <- numSNPLocations
subRes$snpinstances <- numSNPInstances
subRes$point51 <- length(which(data2$FLOW_VALUE == 0.51))

#many cases of 0.51 as a sub

write.table(subRes, file=paste(parentDir, "substitution_rates.csv", sep=""), sep=",")

subMatrix <- table(data2$REF_NUCLEOTIDE, toupper(data2$READ_NUCLEOTIDE))
subMatrixProp <- subMatrix / (sum(subMatrix))

write.table(subMatrixProp, file=paste(parentDir, "proportion_of_substitutions.csv", sep=""), sep=",")

