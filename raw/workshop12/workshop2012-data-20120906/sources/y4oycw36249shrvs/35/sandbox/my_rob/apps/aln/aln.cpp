#include <iostream>
#include <seqan/align.h>



using namespace seqan;

int main (){
  typedef String<char> TString;
  typedef Align<TString,ArrayGaps> TAlign;
  
  TString seq1 = "CDFGTFG";
  TString seq1 = "CDATGDTFG";
  TAlign align;
  resize (rows(align),2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  
  std::cout << align << std::endl;
  
}