load("res_removed_key_column.RData")


library(ggplot2)
theme_set(theme_bw())

label100 <- "Ion OneTouch \n Template Kit"
label200old <- "Ion Xpress \n Template 200 kit"
label200new <- "Ion OneTouch 200 \n Template kit"

#res <- res[sample(1:dim(res)[1], 500000),]


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(2,3,4,8)]

 load("combined_qual_by_hp.RData")

res <- res[which(res$position > 1),]
res <- res[which(res$position < 5),]

randomCorr <- sample(which(res$type == "CORRECT"), 1000000) 
randomIns <- sample(which(res$type == "INSERTION"), 1000000) 

mySub <-res[c(randomCorr, randomIns),]

mySub$group <- ""
mySub$group[which(mySub$kit == "200")] <- label200old
mySub$group[which(mySub$kit == "100")] <- label100
mySub$group[which(mySub$kit == "200Plus")] <- label200new
mySub$group <- factor(mySub$group,levels=c(label100, label200old, label200new))

names(mySub)[1] <- "FoldChange"

mySub$type <- factor(mySub$type, levels=c("CORRECT", "INSERTION"))

m <- ggplot(mySub, aes(x=position, y = FoldChange, fill=type, colour=type)) + facet_wrap(~ group,nrow = 3,scale = "free_y")  + 
stat_sum_df("mean_cl_normal", geom = "smooth") +  scale_colour_manual( values=cbPalette)  + scale_fill_manual( values=cbPalette) +
xlab("Error type") + ylab("Error rate")+ opts(axis.text.x=theme_text(angle=-90))  + 
opts(strip.text.x = theme_text(size = 11, colour = "black"))

pdf("quality_scores_decrease_regardless_of_accuracy.pdf", height=6, width=8)
m
dev.off()