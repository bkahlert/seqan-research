#!/usr/bin/env Rscript
table <- read.table("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/per_pos_qual.tsv", sep="\t",header=TRUE)
x<-boxplot(table$Median~table$Base,plot=FALSE)
x$stats[2,]<-table$Lower.Quartile 
x$stats[1,]<-table$X10th.Percentile
x$stats[3,]<-table$Median
x$stats[5,]<-table$X90th.Percentile
x$stats[4,]<-table$Upper.Quartile
x$names<-table$Base
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_quality.png")
bxp(x,medlwd=1,medcol="red",boxfill = "yellow", border =1, ylim=c(0,max(table$X90th.Percentile)),xlim=c(0,length(table$Median)),main="Qualityscore across all bases")
lines(table$Mean,lwd=1,col="blue")
dev.off()
table <- read.table("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/seq_qual.tsv", sep="\t",header=TRUE)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_sequence_quality.png")
plot(table$Quality,table$Count,type="l",col="red",xlab="Mean Sequence Quality",ylab="",main="Qualityscore distribution over all sequences")
legend("topright",legend="Average Quality per read",cex=.5,text.col="red")
dev.off()
table <- read.table("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/pos_gc_content.tsv", sep="\t",header=TRUE)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_gc_content.png")
table$Base
plot(table$Base,table$GC,type="l",col="red",xlab="Position in read (bp)",ylab="",main="GC content across all bases",ylim=c(0,100))
legend("topright",legend="%GC",cex=.5,text.col="red")
dev.off()
