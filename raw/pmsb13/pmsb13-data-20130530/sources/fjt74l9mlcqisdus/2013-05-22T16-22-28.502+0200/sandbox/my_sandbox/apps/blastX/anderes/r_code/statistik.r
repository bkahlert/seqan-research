loadtable <- function(x){ 
  daten <- read.table(x,sep="\t",header=TRUE)
  return(daten)
}

vergleich <- function(x,y){
  x_read <- x[,1]
  x_prot <- x[,2]
  x_begin <- x[,3]
  x_end <-x[,4]
  
  y_read <- y[,1]
  y_prot <- y[,3]
  y_begin <- y[,4]
  y_end <-y[,5]
  
  richtig <- c()
  for (i in 1:length(x_read)){
    for (j in 1:length(y_read)){
      if (x_read[i]==y_read[j]){
        if (x_prot[i]==y_prot[j]){
          if (x_end[i]==y_end[j] & x_begin[i]==y_begin[j]){
            richtig <- c(richtig,i)
          }
          else{
            print("fehler")
          } 
        }
        else print("falsches_prot")
      }
    }
  }
  return (richtig)
}


main<- function(){
  
  statistik <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/statistik.txt")
  blosum_30 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/matchBlosum30.txt")
  blosum_45 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/matchBlosum45.txt")
  blosum_80 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/matchBlosum80.txt")
  pam120 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/matchPam120.txt")
  pam200 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/matchPam200.txt")
  pam250 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/matchPam250.txt")
  pam40 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/matchPam40.txt")
  Vtml200 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/matchVtml200.txt")
  length(vergleich(statistik,blosum_30))
  length(vergleich(statistik,blosum_45))
  length(vergleich(statistik,blosum_80))
  length(vergleich(statistik,pam120))
  length(vergleich(statistik,pam200))
  length(vergleich(statistik,pam250))
  length(vergleich(statistik,pam40))
  length(vergleich(statistik,Vtml200))
}
main()