loadtable <- function(x){ 
  daten <- read.table(x,sep="\t",header=TRUE)
  return(daten)
}

vergleich <- function(x,y,bin){
  x_read <- x[,1]
  x_prot <- x[,2]
  x_begin <- x[,3]
  x_end <-x[,4]
  
  y_read <- y[,1]
  y_prot <- y[,3]
  y_begin <- y[,4]
  y_end <-y[,5]
  
  richtig <- c()
  false <- c()
  for (i in 1:length(x_read)){
    for (j in 1:length(y_read)){
      if (x_read[i]==y_read[j]){
        if (x_prot[i]==y_prot[j]){
          if (abs(x_end[i]-y_end[j])<=2 & abs(x_begin[i]-y_begin[j])<=2){
            richtig <- c(richtig,i)
          }
          else{
            print("new")
            print(abs(x_end[i]-y_end[j]))
            print(abs(x_begin[i]-y_begin[j]))
            false <- c(false,i)
          }
        }
      }
    }
  }
  if (bin==1) return (richtig)
  else return (false)
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
  
  statistik <- loadtable("D:/seqan_dev/Development/seqan-trunk/sandbox/my_sandbox/apps/blastX/anderes/r_code/data/k2/statistik.txt")
  blosum_30 <- loadtable("D:/seqan_dev/Development/seqan-trunk/sandbox/my_sandbox/apps/blastX/anderes/r_code/data/k2/matchBlosum30.txt")
  blosum_45 <- loadtable("D:/seqan_dev/Development/seqan-trunk/sandbox/my_sandbox/apps/blastX/anderes/r_code/data/k2/matchBlosum45.txt")
  blosum_80 <- loadtable("D:/seqan_dev/Development/seqan-trunk/sandbox/my_sandbox/apps/blastX/anderes/r_code/data/k2/matchBlosum80.txt")
  pam120 <- loadtable("D:/seqan_dev/Development/seqan-trunk/sandbox/my_sandbox/apps/blastX/anderes/r_code/data/k2/matchPam120.txt")
  pam200 <- loadtable("D:/seqan_dev/Development/seqan-trunk/sandbox/my_sandbox/apps/blastX/anderes/r_code/data/k2/matchPam200.txt")
  pam250 <- loadtable("D:/seqan_dev/Development/seqan-trunk/sandbox/my_sandbox/apps/blastX/anderes/r_code/data/k2/matchPam250.txt")
  pam40 <- loadtable("D:/seqan_dev/Development/seqan-trunk/sandbox/my_sandbox/apps/blastX/anderes/r_code/data/k2/matchPam40.txt")
  Vtml200 <- loadtable("D:/seqan_dev/Development/seqan-trunk/sandbox/my_sandbox/apps/blastX/anderes/r_code/data/k2/matchVtml200.txt")
  
  found_blosum30 <- vergleich(statistik,blosum_30,1)
  found_blosum45 <- vergleich(statistik,blosum_45,1)
  found_blosum80 <-vergleich(statistik,blosum_80,1)
  found_pam120 <-vergleich(statistik,pam120,1)
  found_pam200 <-vergleich(statistik,pam200,1)
  found_pam250 <-vergleich(statistik,pam250,1)
  found_pam40 <-vergleich(statistik,pam40,1)
  found_vtml200 <-vergleich(statistik,Vtml200,1)
  
  
  
  not_found_blosum30 <- vergleich(statistik,blosum_30,0)
  not_found_blosum45 <- vergleich(statistik,blosum_45,0)
  not_found_blosum80 <-vergleich(statistik,blosum_80,0)
  not_found_pam120 <-vergleich(statistik,pam120,0)
  not_found_pam200 <-vergleich(statistik,pam200,0)
  not_found_pam250 <-vergleich(statistik,pam250,0)
  not_found_pam40 <-vergleich(statistik,pam40,0)
  not_found_vtml200 <-vergleich(statistik,Vtml200,0)
  
  two <- vereinigung(found_blosum30,found_blosum45)
  three <- vereinigung(two,found_blosum80)
  four <- vereinigung(three,found_pam120)
  fife <- vereinigung(four,found_pam200)
  six <- vereinigung(fife,found_pam250)
  seven<- vereinigung(six,found_pam40)
  eight <- vereinigung(seven,found_vtml200)
  
  #ntwo <- vereinigung(not_found_blosum30,not_found_blosum45) 
  #nthree <- vereinigung(ntwo,not_found_blosum80)
  #nfour <- vereinigung(nthree,not_found_pam120)
  #nfife <- vereinigung(nfour,not_found_pam200)
  #nsix <- vereinigung(nfife,not_found_pam250)
  #nseven<- vereinigung(nsix,not_found_pam40)
  #neight <- vereinigung(nseven,not_found_vtml200)
  
  
  result <- c(length(found_blosum30),length(two),length(three),length(four),length(fife),length(six),length(seven),length(eight))
  print(result)
  #result2<- c(length(not_found_blosum30),length(ntwo),length(nthree),length(nfour),length(nfife),length(nsix),length(nseven),length(neight))
  
}

main()