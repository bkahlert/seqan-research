#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include <fstream>

using namespace seqan;
using namespace std;
/*
String<char> getMasterId(String<char> RefSeq){
  String<char> MasterId;
  bool RefRef=0; 
  for (int i=0;i<length(RefSeq);i++){
    if (toChar(RefSeq[i])=="$")
      RefRef++;
    if (RefRef==1)
      appendValue(MasterId,RefSeq[i]);
  }
  return MasterId;
}
*/  
int main(){
  String<char> SampleString=">seq1 $seq1$";
  cout << SampleString[2];
  return 0;
  
}
