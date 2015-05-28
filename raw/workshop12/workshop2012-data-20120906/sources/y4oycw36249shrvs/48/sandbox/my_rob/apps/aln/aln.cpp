#include <iostream>
#include <seqan/align.h>



using namespace seqan;



void msaDemo(){
  typedef String<AminoAcid> TString;
  typedef StringSet<TString> TStringSet;
  typedef Graph<Alignment<StringSet<TString,Dependent<>>>> TAlingGraph;
  
  TStringSet seq;
    appendValue(seq,"DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE");
    appendValue(seq,"RVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK");
    appendValue(seq,"FPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK");
    appendValue(seq,"HIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK");
  TAlignGreph alignG(seq);
  globalMsaAlignment(alinG,Blosum62(-1,-11));
  

  std::cout << "MSA alignment:" << std::endl;
  std::cout << alignG << std::endl;
}

void globalAln(){
  typedef String<char> TString;
  typedef Align<TString> TAlign;
  TString seq1 = "blablubalu";
  TString seq2 = "abba";
  
  TAlign align;
  resize(rows(align),2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  
  Score<int,Simple> scoring(1,-1,-1,-1);
  int score=globalAlignment(align,scoring,AlignConfig<true,true,true,true>());
  std::cout << "alignment:" << std::endl;
  std::cout << score << std::endl;
  std::cout << align << std::endl;
}

void localAln(){
  typedef String<Dna> TString;
  typedef Align<TString> TAlign;
  TString seq1 = "ataaGcgtctcg";
  TString seq2 = "tcatagagttgc";
  
  TAlign align;
  resize(rows(align),2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  
  Score<int,Simple> scoring(2,-1,-2,0);
  int score=localAlignment(align,scoring);
  std::cout << "Local_alignment:" << std::endl;
  std::cout << score << std::endl;
  std::cout << align << std::endl;
  LocalAlignmentEnumerator<Score<int>,Unbanded> enumerator(scoring,5);
  while(nextLocalAlignment(align,enumerator)){
	std::cout << "Local_alignment:" << std::endl;
	std::cout << getScore(enumerator) << std::endl;
	std::cout << align << std::endl;
  }
}



int main (){
  typedef String<char> TString;
  typedef Align<TString,ArrayGaps> TAlign;
  
  TString seq1 = "CDFGTFG";
  TString seq2 = "CDATGDTFG";
  TAlign align;
  resize (rows(align),2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2); 	
  
  insertGap(row(align,0),2);
  insertGaps(row(align,0),5,2);
  
  std::cout << align << std::endl;
  
  std::cout << "Source" << std::endl;
  for (unsigned i=0;i<length(row(align,0));i++){
	  std::cout << toSourcePosition(row(align,0),i) << " ";
  }
  std::cout << std::endl;
  for (unsigned i=0;i<length(row(align,1));i++){
	  std::cout << toSourcePosition(row(align,1),i) << " ";
  }
  std::cout << std::endl;
  std::cout << "Source" << std::endl;
  for (unsigned i=0;i<length(source(row(align,0)));i++){
	  std::cout << toViewPosition(row(align,0),i) << " ";
  }
  std::cout << std::endl;
  for (unsigned i=0;i<length(source(row(align,1)));i++){
	  std::cout << toViewPosition(row(align,1),i) << " ";
  }
  std::cout << std::endl;
  
  globalAln();
  localAln();
  mseDemo();
  return 1;
  
}