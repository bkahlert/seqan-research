#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    
  
  
  SimpleScore scoringScheme(2, -1, -1, -2);
LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, 5, -3, 0);
 
Dna5String seqH = "CGAGAGAGACCGAGA";
Dna5String seqV = "TTCTGAGATCCGTTTTT";
 
Align<Dna5String> align;
resize(rows(align), 2);
assignSource(rows(align)[0], seqH);
assignSource(rows(align)[1], seqV);
 
int i = 0;
while (nextLocalAlignment(align, enumerator))
{
    std::cout << i << "-th alignment:\n";
    std::cout << align << "\n\n";
    std::cout << "score == " << getScore(enumerator) << "\n";
}
  
  
  
  /*typedef String<Dna> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "AGTCGGATCTACTG";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);
     int score = globalAlignment(align, Score<int,Simple>(4,-2,-2,-4));
    //::std::cout << "Score: " << score << ::std::endl;
    for(int i=0;i<length(rows(align)[0]);i++)
     if((rows(align)[0][i])=='-'){
     ::std::cout << align[0][i] << ::std::endl;}*/
    return 0;
}