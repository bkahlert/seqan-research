#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
  typedef String<char> TSequence;		//sequence type
  typedef Align<TSequence,ArrayGaps> TAlign;	//align Type
  
  TSequence seq1="AAATGACGGATTG";
  TSequence seq2="AGTCGGATCTACTG";
  
  TAlign align;
  resize(rows(align),2);
  assignSource(row(align,0),seq1);
  assignSource(row(align,1),seq2);
  
  int score = globalAlignment(align, Score<int,Simple>(4,-2,-4,-2));
    ::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << align << ::std::endl;

    return 0;
}

  
  