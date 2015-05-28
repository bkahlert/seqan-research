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
          
        }
        
      }
    }
  }
  return (richtig)
}

vereinigung <- function(x,y){
  vereinigung<-c(y)
  count <- 0
  for (i in 1:length(x)){
    for (j in 1:length(y)){
      if (x[i]==y[j]) count <- 1
    }
    if (count == 0) vereinigung <- c(vereinigung,x[i])
    count <- 0
  }
  return(vereinigung)
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
  
  found_blosum30 <- vergleich(statistik,blosum_30)
  found_blosum45 <- vergleich(statistik,blosum_45)
  found_blosum80 <-vergleich(statistik,blosum_80)
  found_pam120 <-vergleich(statistik,pam120)
  found_pam200 <-vergleich(statistik,pam200)
  found_pam250 <-vergleich(statistik,pam250)
  found_pam40 <-vergleich(statistik,pam40)
  found_vtml200 <-vergleich(statistik,Vtml200)
  
  two <- vereinigung(found_blosum30,found_blosum45)
  three <- vereinigung(two,found_blosum80)
  four <- vereinigung(three,found_pam120)
  fife <- vereinigung(four,found_pam200)
  six <- vereinigung(fife,found_pam250)
  seven<- vereinigung(six,found_pam40)
  eight <- vereinigung(seven,found_vtml200)
  
  print(length(two))
  print(length(three))
  print(length(four))
  print(length(fife))
  print(length(six))
  print(length(seven))
  print(length(eight))
 
  
  
}
main()