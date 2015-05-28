#!/usr/bin/env Rscript
table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/per_pos_qual.tsv", sep="\t",header=TRUE,skip=1)
x<-boxplot(table$Median~table$Base,plot=FALSE)
x$stats[2,]<-table$Lower.Quartile 
x$stats[1,]<-table$X10th.Percentile
x$stats[3,]<-table$Median
x$stats[5,]<-table$X90th.Percentile
x$stats[4,]<-table$Upper.Quartile
x$names<-table$Base
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_quality.png")
bxp(x,medlwd=1,medcol="red",boxfill = "yellow", border =1, ylim=c(0,max(table$X90th.Percentile)),xlim=c(0,length(table$Median)),main="Qualityscore across all bases")
rect(0,0,max(table$Base)+2,20,col="lightcoral",border="lightcoral")
rect(0,20,max(table$Base)+2,28,col="lightgoldenrod",border="lightgoldenrod")
rect(0,28,max(table$Base)+2,max(table$X90th.Percentile)+3,col="lightgreen",border="lightgreen")
bxp(x,medlwd=1,medcol="red",boxfill = "yellow", border =1, ylim=c(0,max(table$X90th.Percentile)),xlim=c(0,length(table$Median)),main="Qualityscore across all bases",add=TRUE)
lines(table$Mean,lwd=1,col="blue")
grid(nx=max(table$Base),ny=0,lty=1)
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/seq_qual.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_sequence_quality.png")
plot(table$Quality,table$Count,type="l",col="red",xlab="Mean Sequence Quality",ylab="",main="Qualityscore distribution over all sequences")
legend("topright",legend="Average Quality per read",cex=.5,text.col="red")
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/pos_gc_content.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_gc_content.png")
plot(table$X.Base,table$GC,type="l",col="red",xlab="Position in read (bp)",ylab="",main="GC content across all bases",ylim=c(0,100))
legend("topright",legend="%GC",cex=.5,text.col="red")
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/pos_nuc_dist.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_sequence_content_old.png")
plot(table$X.Base,table$G,type="l",col="red",ylim=c(0,100),xlab="Mean Sequence Quality",ylab="",main="Sequence content accross all bases")
lines(table$A,lwd=1,col="blue")
lines(table$T,lwd=1,col="green")
lines(table$C,lwd=1,col="black")
legend("topright",legend=c("%G","%A","%T","%C"),cex=.5,text.col=c("red","blue","green","black"))
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/seq_gc_content.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_sequence_gc_content.png")
plot(table$X.GC.Content,table$Count,type="l",col="red",xlab="Mean GC content (%)",ylab="",main="GC distribution over all sequences")
legend("topright",legend=c("GC count per read"),cex=.5,text.col=c("red"))
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/pos_n_content.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_n_content.png")
plot(table$X.Base,table$GC,type="l",col="red",xlab="Position in read (bp)",ylab="",main="N content across all bases",ylim=c(0,100))
legend("topright",legend="%N",cex=.5,text.col="red")
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/readlength.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/sequence_length_distribution.png")
x <- c(min(table$X.Length)-1,table$X.Length,max(table$X.Length)+1)
y <- c(0,table$Count,0)
z <- seq(0,max(y),by=max(y)/10)
plot(x,y,xaxt="n",yaxt="n",type="l",col="red",xlab="Sequence length (bp)",ylab="",main="Distribution of sequences length over all sequences",panel.first = grid(nx=0,ny=5,lty="dotted"))
axis(2, at=z,labels=z, col.axis="black", las=1)
axis(1, at=x,labels=x, col.axis="black", las=1)
legend("topright",legend="Sequence Length",cex=.5,text.col="red")
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/kmerContent.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/kmer_content.png")
plot(table$X.Base,table$GC,type="l",col="red",xlab="Position in read (bp)",ylab="",main="N content across all bases",ylim=c(0,100))
legend("topright",legend="%N",cex=.5,text.col="red")
dev.off()

