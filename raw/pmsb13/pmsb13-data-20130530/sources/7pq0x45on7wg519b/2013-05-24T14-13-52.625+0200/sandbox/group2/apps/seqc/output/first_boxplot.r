#!/usr/bin/env Rscript
table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/per_pos_qual.tsv", sep="\t",header=TRUE,skip=1)
x<-boxplot(table$Median~table$Base,plot=FALSE)
x$stats[2,]<-table$Lower.Quartile 
x$stats[1,]<-table$X10th.Percentile
x$stats[3,]<-table$Median
x$stats[5,]<-table$X90th.Percentile
x$stats[4,]<-table$Upper.Quartile
x$names<-table$Base
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_quality.png",width=800,height=600)
bxp(x,medlwd=1,medcol="red",boxfill = "yellow", border =1, ylim=c(0,max(table$X90th.Percentile)),xlim=c(0,length(table$Median)),main="Qualityscore across all bases")
rect(0,0,max(table$Base)+2,20,col="lightcoral",border="lightcoral")
rect(0,20,max(table$Base)+2,28,col="lightgoldenrod",border="lightgoldenrod")
rect(0,28,max(table$Base)+2,max(table$X90th.Percentile)+3,col="lightgreen",border="lightgreen")
bxp(x,medlwd=1,medcol="red",boxfill = "yellow", border =1, ylim=c(0,max(table$X90th.Percentile)),xlim=c(0,length(table$Median)),main="Qualityscore across all bases",add=TRUE)
lines(table$Mean,lwd=1,col="blue")
grid(nx=max(table$Base)/2,ny=0,lty=1)
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/seq_qual.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_sequence_quality.png",width=800,height=600)
plot(table$X.Quality,table$Count,type="l",col="red",xlab="Mean Sequence Quality",ylab="",main="Qualityscore distribution over all sequences")
legend("topright",legend="Average Quality per read",text.col="red")
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/pos_gc_content.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_gc_content.png",width=800,height=600)
plot(table$Base,table$GC,type="l",col="red",xlab="Position in read (bp)",ylab="",main="GC content across all bases",ylim=c(0,100))
legend("topright",legend="%GC",text.col="red")
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/pos_nuc_dist.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_sequence_content_old.png",width=800,height=600)
plot(table$X.Base,table$G,type="l",col="red",ylim=c(0,100),xlab="Mean Sequence Quality",ylab="",main="Sequence content accross all bases")
lines(table$A,lwd=1,col="blue")
lines(table$T,lwd=1,col="green")
lines(table$C,lwd=1,col="black")
lines(table$N,lwd=1,col="orange")
legend("topright",legend=c("%G","%A","%T","%C","%N"),text.col=c("red","blue","green","black","orange"))
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/seq_gc_content.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_sequence_gc_content.png",width=800,height=600)
plot(table$X.GC,table$Count,type="l",col="red",xlab="Mean GC content (%)",ylab="",main="GC distribution over all sequences")
legend("topright",legend=c("GC count per read"),text.col=c("red"))
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/pos_n_content.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/per_base_n_content.png",width=800,height=600)
plot(table$X.Base,table$X.GC,type="l",col="red",xlab="Position in read (bp)",ylab="",main="N content across all bases",ylim=c(0,100))
legend("topright",legend="%N",text.col="red")
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/readlength.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/sequence_length_distribution.png",width=800,height=600)
x <- c(min(table$X.Length)-1,table$X.Length,max(table$X.Length)+1)
y <- c(0,table$Count,0)
z <- seq(0,max(y),by=max(y)/10)
plot(x,y,xaxt="n",yaxt="n",type="l",col="red",xlab="Sequence length (bp)",ylab="",main="Distribution of sequences length over all sequences",panel.first = grid(nx=0,ny=5,lty="dotted"))
axis(2, at=z,labels=z, col.axis="black", las=1)
axis(1, at=x,labels=x, col.axis="black", las=1)
legend("topright",legend="Sequence Length",text.col="red")
dev.off()

table <- read.csv("/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/expected/kmerContent.tsv", sep="\t",header=TRUE,skip=1)
png(file="/home/dkersting/Development/seqan-trunk/sandbox/group2/apps/seqc/output/Images/kmer_content.png",width=800,height=600)
ylimvar<-c(table[,2],table[,3],table[,4],table[,5],table[,6],table[,7])
plot(table$Position,table[,2],lwd=3,type="l",col="red",xlab="Position in read (bp)",ylab="",main="Absolute enrichment over read length",ylim=c(0,max(ylimvar)))
lines(table[,3],lwd=3,col="blue")
lines(table[,4],lwd=3,col="green")
lines(table[,5],lwd=3,col="black")
lines(table[,6],lwd=3,col="orange")
lines(table[,7],lwd=3,col="magenta")
ltv<-colnames(table)
legendtext<-c(ltv[2],ltv[3],ltv[4],ltv[5],ltv[6],ltv[7])
legend("topright",legend=legendtext,text.col=c("red","blue","green","black","orange","magenta"),bg="gray95",horiz=TRUE)
dev.off()

