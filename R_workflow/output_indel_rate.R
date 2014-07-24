
tmp <- read.csv((paste(parentDir, "flow_data_expanded.csv", sep=""), header=TRUE) #is there a header?
res <- c()
res$insertErrors <- with(tmp, length(which(tmp$Type == "Insertion")))
res$deletionErrors <- with(tmp, length(which(tmp$Type == "Deletion")))
res$correctHPs <- with(tmp, length(which(tmp$Type == "Correct")))
res$total <- res$insertError + res$deletionErrors + res$correctHPs
res$insertRate <- res$insertErrors / res$total
res$deletionRate <- res$deletionErrors / res$total
write.table(res, paste(parentDir, "number_insertion_deletion_errors.csv", sep=""), sep=",");

#only part of it.

# Next -- we need the error-rate by base

tmp$Error <- tmp$Type != "Correct";
library(data.table)

dt = data.table(tmp)

errorRateByBasePosPhase = dt[, list(error_rate = mean(Error)), by=c("BasePosition", "OOP")]
errorRateByBasePosValidFV = dt[, list(error_rate = mean(Error)), by=c("BasePosition", "ValidFV")]
errorsPerReadPhase = dt[, list(error_rate = sum(Error)), by=c("ReadID", "OOP")]
errorsPerReadValidFV = dt[, list(error_rate = sum(Error)), by=c("ReadID", "ValidFV")]

#write.table(errorRateByBaseposPhase, file=paste(parentDir, "error_rate_by_base_position_phase.csv", sep="")); #what name does it need again
#write.table(errorRateByBasePosValidFV, file=paste(parentDir, "error_rate_by_base_position_fv.csv", sep="")); #what name does it need again

#write.table(errorsPerReadPhase, file=paste(parentDir, "num_errors_per_reads_phase.csv", sep=""));
#write.table(errorsPerReadValidFV, file=paste(parentDir, "num_errors_per_reads_fv.csv", sep=""));

#errorsPerReadPhase = read.table("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/ecoli_no_calibration_test/num_errors_per_reads_phase.csv")
#errorsPerReadFV = read.table("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/ecoli_no_calibration_test/num_errors_per_reads_fv.csv")

setnames(errorsPerReadPhase, "OOP", "State")
errorsPerReadPhase$by = "OOP"
setnames(errorsPerReadValidFV, "ValidFV", "State")
errorsPerReadValidFV$by = "ValidFV"

comb_errors_per_read = rbind(errorsPerReadPhase, errorsPerReadValidFV)

write.table(comb_errors_per_read, file=paste(parentDir, "num_errors_per_read_with_phase_fv.csv", sep=""))

errorRateByBasePosPhase[, by := "Phase"]
setnames(errorRateByBasePosPhase, "OOP", "State")

errorRateByBasePosValidFV[, by := "ValidFV"]
setnames(errorRateByBasePosValidFV, "ValidFV", "State")

comb = rbind(errorRateByBasePosPhase, errorRateByBasePosValidFV)

write.table(comb, file=paste(parentDir, "error_rate_by_base_position.csv", sep=""))

