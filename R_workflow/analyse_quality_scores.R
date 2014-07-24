library(plyr)

library(data.table)

rep.mask1 <- unlist(apply(repeatLocations ,1,function(x) seq(x[1],x[1] + x[3] - 1)))
rep.mask2 <- unlist(apply(repeatLocations ,1,function(x) seq(x[2],x[2] + x[3] - 1)))
rep.mask <- c(rep.mask1, rep.mask2)


tmp <- read.table(paste(parentDir, "rawData/errors2basepos.txt", sep=""), header=TRUE, sep="\t")
names(tmp) <- c("Read","BasePosition","FlowPosition","RLEPosition", "RefPosition","Quality", "Strand", "ValidFV", "OOP", "Type")
tmp = tmp[!is.na(tmp$BasePosition),]
tmp$BasePosition <- as.numeric(as.character(tmp$BasePosition)) + 1
tmp = tmp[with(tmp, BasePosition >= startPos),]
tmp = tmp[with(tmp, BasePosition <= truncatePos),]
tmp <- tmp[with(tmp, Type != "DELETION"),]

indelLocs <- read.table(paste(parentDir, "possible_indel_errors_in_ref.txt", sep=""), sep=",")
snpLocs <- read.table(paste(parentDir, "polymorphism_locations.csv", sep=""), sep=",") #should be a single vector
indelLocsUniq <-  unique(indelLocs$RefRealPos)
snpLocsUniq <-  unique(unlist(snpLocs))
tmp <- tmp[with(tmp, !(RefPosition %in% rep.mask)),]
tmp <- tmp[with(tmp, !(RefPosition %in% indelLocsUniq)),]
tmp <- tmp[with(tmp, !(RefPosition %in% snpLocsUniq)),]
tmp$Quality <- as.numeric(tmp$Quality)

## finally in the interesting bit.
## Might want to see how indicative the quality is of error rate, by Type.
## save tmp perhaps.
## does this tabled qual score even work?
#instead, let's try and get confirmation that what I am seeing is wrong.

#which(tmp$Quality == 5 & tmp$by == "OO")

tableQualScores = function(d)  
{
  
  tabledQual <- with(d, table(Quality,Type))
  accuracy <- tabledQual[,1] / (tabledQual[,1] + tabledQual[,3])
  errorProb <- 1 - accuracy
  qualityScores <- log(errorProb, base=10) * -10
  itQuals <- names(qualityScores)
  names(qualityScores) <- c()
  toSave <- cbind(qualityScores, itQuals);
  return(toSave)
}

res_valid = tableQualScores(tmp[which(tmp$ValidFV == "True"),])
res_invalid = tableQualScores(tmp[which(tmp$ValidFV == "False"),])
res_in_phase = tableQualScores(tmp[which(tmp$OOP == "False"),])
res_oo_phase = tableQualScores(tmp[which(tmp$OOP == "True"),])

## needs to be backed up
#save(res_valid, file=paste(parentDir, "quality_scores_valid_fv.csv", sep=""));
#save(res_invalid, file=paste(parentDir, "quality_scores_invalid_dv.csv", sep=""));
#save(res_in_phase, file=paste(parentDir, "quality_scores_inphase.csv", sep=""));
#save(res_oo_phase, file=paste(parentDir, "quality_scores_oop.csv", sep=""));

res_valid = as.data.frame(res_valid)
res_invalid = as.data.frame(res_invalid)
res_in_phase = as.data.frame(res_in_phase)
res_oo_phase = as.data.frame(res_oo_phase)


##Friday morning. Load all the above, and make a combined data.frame, stupid to have so many outfiles
#res_valid = read.table("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/ecoli_no_calibration_test/quality_scores_valid_fv.csv", header=TRUE, sep=",")
#res_invalid = read.table("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/ecoli_no_calibration_test/quality_scores_invalid_dv.csv", header=TRUE, sep=",")
#res_in_phase = read.table("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/ecoli_no_calibration_test/quality_scores_inphase.csv", header=TRUE, sep=",")
#res_oo_phase = read.table("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/ecoli_no_calibration_test/quality_scores_oop.csv", header=TRUE, sep=",")

#make one mega table (not really mega)
res_valid$by = "ValidFV"
res_valid$state = "True"
res_invalid$by = "ValidFV"
res_invalid$state = "False"
res_in_phase$by = "OOP"
res_in_phase$state = "False"
res_oo_phase$by = "OOP"
res_oo_phase$state = "True"

res_quality_scores = rbindlist(list(res_valid, res_invalid, res_in_phase, res_oo_phase))
write.table(res_quality_scores, file=paste(parentDir, "quality_scores_by_phase_fv.csv"), sep=",")


