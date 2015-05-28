loadtable <- function(x){ 
  daten <- read.table(x,sep="\t",header=TRUE)
  return(daten)
}

vergleich <- function(x,y){
  
}


main<- function(){
  
  statistik <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/statistik.txt")
  blosum_30 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/outputBlosum30verify.txt")
  
}