#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<Dna> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "AGTCGGATCTACTG";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);
     int score = globalAlignment(align, Score<int,Simple>(4,-2,-2,-4));
    //::std::cout << "Score: " << score << ::std::endl;
    //for(int i=0;i<length(rows(align)[0]);i++)
     //if((rows(align)[0][i])=='-'){
     ::std::cout << align << ::std::endl;
  //  ::std::cout << begin(align) << ::std::endl;
    return 0;
}