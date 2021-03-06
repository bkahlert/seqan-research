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
  for (i in 1:10){
    for (j in 1:10){
      if (x_read[i]==y_read[j]){
        if (x_prot[i]==y_prot[j]){
          if (x_begin[i]==y_begin[j]){
            richtig <- c(richtig,i)
          }
          else{
            print("falscher anfang")
            print(x_begin[i])
            print(y_begin[j])
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
  blosum_30 <- loadtable("D:/seqan_dev/Development/seqan-trunk-build/vs9/sandbox/my_sandbox/apps/blastX/outputBlosum30verify.txt")
  vergleich(statistik,blosum_30)
}
main()