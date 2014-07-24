#Figure 1

## Up here, would set up the directory for plotting (maybe this should be done from a configuration file)
setwd("/home/bra427/Projects/IonTorrentBenchmarking/full_run_files/combined_results_inhouse_400bp/")

#so first, we look at consolidated error counts
data <- read.table("consolidated_error_counts.csv", sep=",", header=TRUE) #new format.

data$KIT = "New400bp"
data$DESC = "400bp sequencing kit TS 3.4"

#old data
data2 <- read.table("/home/bra427/Projects/IonTorrentBenchmarking/resultsIncl400bp/consolidated_error_counts.csv", sep=",", header=TRUE)
data2 = data2[, c("CHIP", "DESC", "SEQUENCING_KIT", "SPECIES", "correctHomopolymerMinusPoly", "correctHomopolymers", 
                  "deletionErr", "insertionErr", "num51SNPerr", "subSNPCount", "subSNPLocs", "substitutionErr", "substitutionErrRate")]

names(data2)[3] = "KIT"

combined_data = rbind(data[, names(data2),], data2)


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(1,2,3,4,7,8)]




#Error rates (global)
## should probably look at how correctHomopolymer is calculated.

## this is not working the way I would expect?

combined_data$delRate <- with(combined_data, deletionErr / (correctHomopolymers  + insertionErr))
combined_data$insRate <- with(combined_data, insertionErr / (correctHomopolymers  + deletionErr))
combined_data$subRate <- with(combined_data, substitutionErr / (correctHomopolymerMinusPoly + substitutionErr))


#to make this work in the boxplot, need to expand crap

dataSubs <- combined_data[, c("CHIP", "DESC", "KIT", "SPECIES", "substitutionErrRate")]
dataSubs$errType <- "Substitutions"
names(dataSubs)[5] <- "rate" 
dataIns <- combined_data[, c("CHIP", "DESC", "KIT", "SPECIES", "insRate")]
dataIns$errType <- "Insertions"
names(dataIns)[5] <- "rate"

dataDels <- combined_data[,c("CHIP", "DESC", "KIT", "SPECIES", "delRate")]
dataDels$errType <- "Deletions"
names(dataDels)[5] <- "rate"

combErrData <- rbind(dataSubs, dataIns, dataDels)

combErrData$KIT[which(combErrData$KIT == "New400bp")] = "400bp (new template and sequencing kits)"
combErrData$KIT[which(combErrData$KIT == "New300bp")] = "300bp (new 200bp template kit)"
combErrData$KIT[which(combErrData$KIT == "OneTouch200bp")] = "200bp (old template kit - OneTouch200bp)"
combErrData$KIT[which(combErrData$KIT == "Old200bp")] = "200bp (new 200bp template kit, old 200bp sequencing kit)"

#meanBlock <- with(combData, tapply(rate, list(kit,errType), mean))
#testing <- with(combData[which(combData$errType != "Substitution"),], tapply(rate, list(kit), mean))
#meanBlock
#meanBlock * 100
#apply(meanBlock, MARGIN=1, FUN=sum) * 100
#sdBlock <- with(combData, tapply(rate, list(kit,errType), sd))
#percentOfMean <- sdBlock/ meanBlock
#apply(percentOfMean, MARGIN=2, FUN=mean)

library(scales)
library(ggplot2)

## to be used when I have lots of new datasets.
if(0)
{
  meanErrRateBox <- ggplot(combErrData, aes(x = errType, y=rate, fill=errType)) + facet_wrap(~ DESC,nrow = 1,scale = "free_y")   + 
  geom_boxplot() + coord_cartesian(ylim=c(0.0001,0.07)) + scale_fill_manual(breaks="errType", values=cbPalette)  + scale_y_log10(labels=percent_format()) + 
  xlab("Error type") + ylab("Error rate")  +  theme_bw() + scale_colour_manual(values=cbPalette) +
  theme(strip.text.x = element_text(size = 11, colour = "black")) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + theme(legend.position = "none") + 
  theme(axis.text.x=element_text(angle=-90))
}

## boxplot is not really useful for these datasets. 

## Good enough
svg("avg_rates_of_error_across_datasets.svg",width=8.5,height=6)
ggplot(combErrData, aes(x=errType, y=rate, shape=KIT, col=KIT)) + xlab("Error type") + ylab("Error rate")  +  theme_bw() + 
  scale_colour_manual(values=cbPalette) + geom_point() + scale_y_log10(labels=percent_format()) + theme(axis.text.x=element_text(angle=-90))
+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()

#combining the plotting output of both indels and subs into one file
#Figure 2

## Problems here:



tmpHpErr_phase = read.table("error_rate_by_base_value.csv", sep=",", header=TRUE)
tmpHpErr_phase$run = "Ecoli_new"
tmpHpErr_no_phase = read.table("/home/bra427/Projects/IonTorrentBenchmarking/resultsIncl400bp/error_rate_by_base_value.csv", sep=",", header=TRUE)
tmpHpErr_no_phase$by = "Phase"
tmpHpErr_no_phase$State = "False" #not OOP.
names(tmpHpErr_no_phase)[8] = "run"
names(tmpHpErr_no_phase)[7] = "KIT"
names(tmpHpErr_no_phase)[1] = "RowNumber"
names(tmpHpErr_no_phase)[3] = "BasePosition"
names(tmpHpErr_no_phase)[2] = "error_rate"
tmpHpErr_no_phase$run = as.character(tmpHpErr_no_phase$run)

tmpHpErr_no_phase$run[which(tmpHpErr_no_phase$run == "59")] = "59_99"
tmpHpErr_no_phase$run[which(tmpHpErr_no_phase$run == "82")] = "82_01"
tmpHpErr_no_phase$run[which(tmpHpErr_no_phase$run == "83")] = "83_01"
tmpHpErr_no_phase$run[which(tmpHpErr_no_phase$run == "84")] = "84_01"
tmpHpErr_no_phase$run[which(tmpHpErr_no_phase$run == "Ecoli")] = "Ecoli (IT supplied)"


tmpHpErr_combined = rbind(tmpHpErr_phase[, -1], tmpHpErr_no_phase[, c("BasePosition", "error_rate", "State", "CHIP", "DESC", "KIT", "SPECIES", "by","run")])

library(ggplot2)
theme_set(theme_bw())

tmpHpErr_combined$type <- "Homopolymer"

tmpHpErr_combined[is.na(tmpHpErr_combined)] <- 0

hpErr <- tmpHpErr_combined[,c("error_rate", "BasePosition", "SPECIES", "CHIP", "KIT", "DESC", "by", "State", "run")]
names(hpErr) <- c("rate", "ReadBasePos", "species", "chip", "kit", "desc", "by", "State", "run")

hpErr$error_type = "Homopolymer"

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(1,3,5,8)]

tmpSubErr_phase <- read.table("consolidated_snp_rate_by_base_pos.csv", sep=",", header=FALSE)
tmpSubErr_phase[is.na(tmpSubErr_phase)] <- 0
tmpSubErr_phase$run  = "Ecoli_new"

tmpSubErr_no_phase = read.table("/home/bra427/Projects/IonTorrentBenchmarking/resultsIncl400bp/consolidated_snp_rate_by_base_pos.csv", sep=",", header=FALSE)
tmpSubErr_no_phase[is.na(tmpSubErr_no_phase)] <- 0
names(tmpSubErr_no_phase) = c("CHIP", "DESC", "run", "KIT", "SHORT_RUN", "SPECIES", "SEQUENCING_KIT", "ReadBasePos", "Correct", "Error", "Rate")
tmpSubErr_no_phase$run = as.character(tmpSubErr_no_phase$run)

tmpSubErr_no_phase$run[which(tmpSubErr_no_phase$run == "Ecoli400")] = "Ecoli (IT supplied)"

names(tmpSubErr_phase) = c("CHIP", "DESC",  "KIT", "SPECIES", "ReadBasePos", "Correct", "Error", "Rate", "run")


tmpSubErr = rbind(tmpSubErr_phase, tmpSubErr_no_phase[, c("CHIP", "DESC", "KIT", "SPECIES", "ReadBasePos", "Correct", "Error", "Rate", "run")])
#tmpSubErr <- tmpSubErr[which((tmpSubErr$BasePos <= 99 & tmpSubErr$kit == "100") | (tmpSubErr$BasePos <= 198 & tmpSubErr$kit != "100")),]
#tmpSubErr$species <- as.character(tmpSubErr$species)

tmpSubErr$type <- "Substitution"
subErr <- tmpSubErr[,c("Rate", "ReadBasePos", "SPECIES", "CHIP", "KIT", "DESC", "run")]
subErr$by = "Phase"
subErr$State = "False" ##no OOP
subErr$error_type = "Substitutions"
names(subErr) <- c("rate", "ReadBasePos", "species", "chip", "kit", "desc", "run", "by", "State", "error_type")


res <- rbind(hpErr, subErr)

library(scales)
library(grid)

censored_res = subset(res, (kit == "New400bp" & ReadBasePos < 400) | 
                        (kit == "New300bp" & ReadBasePos < 300) | 
                        (kit == "OneTouch200bp" & ReadBasePos < 200) | 
                        (kit == "Old200bp" & ReadBasePos < 200) | 
                        (kit == "400" & ReadBasePos < 400))

censored_res$kit[which(censored_res$kit == "New400bp")] = "400bp (new template and sequencing kits)"
censored_res$kit[which(censored_res$kit == "New300bp")] = "300bp (new 200bp template kit)"
censored_res$kit[which(censored_res$kit == "OneTouch200bp")] = "200bp (old template kit - OneTouch200bp)"
censored_res$kit[which(censored_res$kit == "Old200bp")] = "200bp (new 200bp template kit, old 200bp sequencing kit)"
censored_res$kit[which(censored_res$kit == "400")] = "400bp (new template and sequencing kits)"


censored_res$kit = factor(censored_res$kit, levels = c("200bp (old template kit - OneTouch200bp)", "200bp (new 200bp template kit, old 200bp sequencing kit)",
                                                       "300bp (new 200bp template kit)", "400bp (new template and sequencing kits)"))

svg("error_rate_hp_versus_substitution.svg",onefile=FALSE,width=10,height=7)

##Phasing is causing problems.#

phase_corr = subset(censored_res, by == "Phase" & State == "False")

hpVersusSub <- ggplot(phase_corr, aes(x=ReadBasePos, y=rate, colour=kit, group=run)) + facet_grid( ~ error_type,scale = "free_y") + 
xlab("Base position") + ylab("Error rate") + coord_cartesian(xlim = c(0, 400), ylim=c(0,0.075)) +
geom_line(size=0.85) + labs(colour = "Species (Chip)")+  theme_bw() + scale_y_continuous(labels=percent_format())+
theme(strip.text.x = element_text(size = 11, colour = "black")) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + theme(panel.margin = unit(0.85, "lines"))

##

hpVersusSub
dev.off()

## What about phase corrected versus un-phase corrected?
new_data_only_phase = subset(censored_res, run =="Ecoli_new" & by == "Phase" &  error_type == "Homopolymer")


hpVersusSub <- ggplot(new_data_only_phase, aes(x=ReadBasePos, y=rate, colour=State)) +
  xlab("Base position") + ylab("Error rate") + coord_cartesian(xlim = c(0, 400), ylim=c(0,0.1)) +
  geom_line(size=0.85) + labs(colour = "Species (Chip)")+  theme_bw() + scale_y_continuous(labels=percent_format())+
  theme(strip.text.x = element_text(size = 11, colour = "black")) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + theme(panel.margin = unit(0.85, "lines"))

svg("error_rate_hp_oop_versus_in_phase.svg",onefile=FALSE,width=10,height=7)
hpVersusSub
dev.off()
#Figure 3


###
data <- read.table("consolidated_substitution_rates_between_bases.csv",sep=",", header=TRUE)
data$run = "Ecoli_new"
##old data

data2 = read.table("/home/bra427/Projects/IonTorrentBenchmarking/resultsIncl400bp/consolidated_substitution_rates_between_bases.csv", header=TRUE, sep=",")
names(data2) = c("CHIP", "DESC", "run", "KIT", "SHORT_NAME", "SPECIES", "TEMPLATE_IT", "RefBase", "ReadBase", "Freq")

comb_data = rbind(data, data2[, c("DESC", "CHIP", "KIT", "run", "SPECIES", "RefBase", "ReadBase", "Freq")])
comb_data$subType = with(comb_data, paste(RefBase, "=>", ReadBase, sep=""))
subFreq <- ggplot(comb_data, aes(x = subType, y=Freq, col=RefBase, shape=run)) + geom_point() + facet_wrap( ~ KIT) +
  theme_bw() + scale_fill_manual(values=cbPalette) +
  theme(strip.text.x = element_text(size = 11, colour = "black")) + 
  theme(panel.grid.minor=element_blank()) +
  theme(axis.text.x=element_text(angle=-45))

svg("substitution_preferences.svg", onefile=FALSE,width=10,height=5)
subFreq
dev.off()
#Can I do a boxplot? Probab.y not

  
# Not this.
if(0)
{
geom_boxplot(data = subset(data, group=run))  
geom_boxplot(data = subset(data, group=label100)) +  
geom_boxplot(data = subset(data, group=label200old)) +
geom_boxplot(data = subset(data, group=label200new)) + 
theme_bw() + scale_fill_manual(values=cbPalette) +
opts(strip.text.x = theme_text(size = 11, colour = "black")) + 
opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank()) +
opts(axis.text.x=theme_text(angle=-45)) + 
xlab("Substitution type (Reference base -> Read base)") + ylab("Proportion of substitutions") + labs(fill = "Reference \nBase") 

svg("substitution_preferences.svg", onefile=FALSE,width=10,height=5)
subFreq
dev.off()
}

#Supp Figure
#despite all these plots, how about I try and visualise read length distributions per datasets, the mean may not be that informative


#dont need to run this again.
readLen <- read.table("consolidated_read_lengths.csv", sep=",", header=TRUE)
names(readLen) <- c("chip", "desc", "kit", "species", "ReadLength", "Count", "TotalReads")

readPlot <- ggplot(readLen, aes(x = ReadLength, y= Count / TotalReads)) +  #+ facet_wrap(~ group,nrow = 1,scale = "free_y")+ 
geom_line()   + coord_cartesian(ylim=c(0, 0.03), xlim=c(0,650)) + 
ylab("Frequency") + xlab("Read length")+ theme_bw() + scale_colour_manual(values=cbPalette) +
theme(strip.text.x = element_text(size = 11, colour = "black")) + 
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

svg("read_length_distribution.svg",onefile=FALSE,width=9,height=6)
readPlot
dev.off()




library(data.table)

# Supp Figure
res <- read.table("error_rate_by_flow_values.csv", sep=",", header=TRUE)
res = res[,-1]


res$KIT = "400bp (new 400bp sequencing kit and new 200bp template kit)"

#two problems here, if I want to include the old data, got to subset it.
dt = data.table(res)

#below works, merge this with old data.
err_by_flow_pos_only = dt[, list(rate = (sum(count.above) + sum(count.below)) / sum(count.at)), by="FlowPos,CHIP,DESC,KIT,SPECIES"]
err_by_flow_pos_only$RUN = "Ecoli 400bp (in-house)"

library(ggplot2)
library(scales)

err_by_flow_pos_only_phase = dt[, list(rate = (sum(count.above) + sum(count.below)) / sum(count.at)), by="FlowPos,CHIP,DESC,KIT,SPECIES,state"]

## plot by OOP versus INPHASE>
errByFlowPos_phase <- ggplot(err_by_flow_pos_only_phase, aes(x = FlowPos)) +
  xlab("Flow position") + 
  ylab("Error rate")  +
  coord_cartesian(xlim = c(0, 800), ylim=c(0,0.20)) + 
  geom_line(data = subset(err_by_flow_pos_only_phase, state="False"), aes(x=FlowPos, y=rate, group=state, colour=state)) + 
  geom_line(data = subset(err_by_flow_pos_only_phase, state="True"), aes(x=FlowPos, y=rate, group=state, colour=state)) + 
  theme_bw() + 
  opts(strip.text.x = theme_text(size = 11, colour = "black")) + 
  opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank()) +
  labs(colour = "Out-of-phase?") + guides(fill=guide_legend(title=NULL)) + opts(panel.margin = unit(0.75, "lines"))

svg("error_rate_by_flow_position_and_phase.svg",onefile=FALSE,width=9,height=6)
errByFlowPos_phase 
dev.off()




data2 = read.table("/home/bra427/Projects/IonTorrentBenchmarking/resultsIncl400bp/error_rate_by_flow_values.csv", header=TRUE, sep=",")
## Old error rate.
data2 = data2[,c("position","CHIP", "DESC", "TEMPLATE_KIT", "SPECIES", "x", "SHORT_NAME")]
names(data2) = c("FlowPos", "CHIP", "DESC", "KIT", "SPECIES", "rate", "RUN")
data2$KIT = as.character(data2$KIT)
data2$RUN = as.character(data2$RUN)



res = rbind(err_by_flow_pos_only, data2)


res$KIT[which(res$DESC == "300bp sequencing kit with 200bp template kit")] = "300bp (new 300bp sequencing kit, new 200bp template kit)"
res$KIT[which(res$DESC == "200bp OneTouch sequencing kit with old 200bp template kit")] = "200bp (old 200bp sequencing kit, old 200bp template kit)"
res$KIT[which(res$DESC == "400bp sequencing kit from IonCommunity Full plate of new 200bp template kit with old 200bp sequencing kit.")] = "400bp (new 400bp sequencing kit and new 200bp template kit)"
res$KIT[which(res$DESC == "Half plate of New 200bp template and old 200bp sequencing kit.")] = "200bp (old 200bp sequencing kit, new 200bp template kit)"
res$KIT[which(res$DESC == "400bp sequencing kit from IonCommunity")] = "400bp (new 400bp sequencing kit and new 200bp template kit)"
res$KIT[which(res$DESC == "Full plate of new 200bp template kit with old 200bp sequencing kit.")] = "200bp (old 200bp sequencing kit, new 200bp template kit)"

res$RUN[which(res$DESC == "400bp sequencing kit from IonCommunity")] = "Ecoli 400bp (Ion Community)"

##Above fixed all the different descriptions


table(res$KIT)

res$KIT = factor(res$KIT, levels=c("200bp (old 200bp sequencing kit, old 200bp template kit)", "200bp (old 200bp sequencing kit, new 200bp template kit)",
                                   "300bp (new 300bp sequencing kit, new 200bp template kit)", "400bp (new 400bp sequencing kit and new 200bp template kit)"))



#res <- res[which(res$position >= 11),]
#res$rate = (res$count.above + res$count.below) / res$count.at

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(2,3,4,6,7,8)]

#things to fix, position should be 'flow position', x= flow error rate"
# the title should be null in the legend
# legend shuold be on done once? 
library(scales)
require(grid)
svg("error_rate_by_flow_pos.svg",onefile=FALSE,width=12,height=6)

errByFlowPosFull <- ggplot(res,aes(x = FlowPos)) +  facet_wrap(~ KIT, scale = "free_y", nrow=2)  +
xlab("Flow position") + 
ylab("Error rate")  +
coord_cartesian(xlim = c(0, 800), ylim=c(0,0.30)) + 
 geom_line(data = subset(res, kit="400bp (new 400bp sequencing kit and new 200bp template kit)"), aes(x=FlowPos, y=rate, group=RUN, colour=RUN)) + 
 geom_line(data = subset(res, kit="300bp (new 300bp sequencing kit, new 200bp template kit)"), aes(x=FlowPos, y=rate, group=RUN, colour=RUN)) + 
 geom_line(data = subset(res, group="200bp (old sequencing kit, old 200bp template kit)"), aes(x=FlowPos, y=rate, group=RUN, colour=RUN)) + 
 geom_line(data = subset(res, group="200bp (old 200bp sequencing kit, new 200bp template kit)"), aes(x=FlowPos, y=rate, group=RUN, colour=RUN)) + 
  theme_bw() + scale_colour_manual(values=cbPalette) +
opts(strip.text.x = theme_text(size = 11, colour = "black")) + 
opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank()) +
labs(colour = "Run ID") + guides(fill=guide_legend(title=NULL)) + opts(panel.margin = unit(0.75, "lines"))

errByFlowPosFull 
dev.off()

## too much data for my computer!
# Supplementary figure also.

if(0)
{
  
#this is the confusing one eep.  
res <- read.table("error_rate_by_error_distribution.csv", sep=",", header=TRUE)

library(data.table)
dt = data.table(res)

dt_sum_by_error_rate = dt[error_rate <= 75, list(num_reads = length(ReadID)), by="by,State,error_rate"]
dt_num_rows_per_treat = dt[error_rate <= 75, list(num_reads_total = length(ReadID)), by="by,State"]
setkey(dt_sum_by_error_rate, by, State)
setkey(dt_num_rows_per_treat, by, State)
comb = dt_sum_by_error_rate[dt_num_rows_per_treat][, freq := num_reads / num_reads_total] 

## probably be better as a percentage of data,
#might be better to go another scale.
library(scales)

ticks <- unique( c(seq(from=0, to=15, by=1), seq(from=15, to=75, by=15)))

##somehow fix this.

m <- ggplot(subset(comb, by =="OOP"), aes(x=error_rate, colour=State, fill=State, group=State,)) 
m = m + geom_bar(aes(y=freq), stat = "identity", position = position_dodge()) +
  scale_y_continuous(labels = percent_format()) + scale_x_continuous(breaks=ticks)

svg("errors_per_read_dist_phased.svg",onefile=FALSE,width=10,height=6)
m
dev.off()


data2 = read.table("/home/bra427/Projects/IonTorrentBenchmarking/resultsIncl400bp/error_rate_by_error_distribution.csv", sep=",", header=TRUE)
dt_2 = data.table(data2)

## grab the other data.
dt_sum_by_error_rate2 = dt_2[x <= 45, list(num_reads = length(X)), by="SHORT_NAME,DESC,x"]
dt_num_rows_per_treat2 = dt_2[x <= 45, list(num_reads_total = length(X)), by="SHORT_NAME,DESC"]
setkey(dt_sum_by_error_rate2, SHORT_NAME, DESC)
setkey(dt_num_rows_per_treat2, SHORT_NAME, DESC)
comb2 = dt_sum_by_error_rate2[dt_num_rows_per_treat2][, freq := num_reads / num_reads_total] 

setnames(comb2, "x", "error_rate")

## reformat the first data, only looking at OOP and error rate < 75 

new_data_num_reads =  dt[by=="OOP" & error_rate <= 45, list(num_reads = length(ReadID)), by="error_rate"]
new_data_num_reads[, num_reads_total := sum(num_reads)]
new_data_num_reads[,freq := num_reads / num_reads_total]
new_data_num_reads[, DESC := "Ecoli 400bp in-house"]

comb2$SHORT_NAME = NULL

setcolorder(comb2, c("error_rate", "num_reads","num_reads_total", "freq", "DESC"))


all_data = rbindlist(list(comb2, new_data_num_reads))

ticks <- unique( c(seq(from=0, to=15, by=1), seq(from=15, to=45, by=15)))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(2,3,4,5,6,8)]

m <- ggplot(all_data, aes(x=error_rate, colour=DESC, fill=DESC)) 
m = m + geom_bar(aes(y=freq), stat = "identity", position = position_dodge()) +
  scale_y_continuous(labels = percent_format()) + scale_x_continuous(breaks=ticks) + theme_bw()

m

svg("errors_per_read_dist_phased_versus_old.svg",onefile=FALSE,width=10,height=6)
m
dev.off()

