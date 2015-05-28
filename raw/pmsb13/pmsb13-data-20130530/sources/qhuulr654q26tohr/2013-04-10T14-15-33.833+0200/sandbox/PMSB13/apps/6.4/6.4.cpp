#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<AminoAcid> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
    
    TSequence seq1 = "PNCFDAKQRTASRPL";
    TSequence seq2 = "CFDKQKNNRTATRDTA";
    
    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);
    
    Score<int> scoring(3,-2,-1,-5);
    LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(scoring);
    
    for (int i=0;i<3;++i)
    {
	nextLocalAlignment(align, enumerator)
	std::cout << "Score = " << getScore(enumerator) << std::endl;
	std::cout << align;
	std::cout << "Aligns Seq1[" << clippedBeginPosition(row(align, 0)) << ":" << (clippedEndPosition(row(align, 0))-1) << "]";
	std::cout << " and Seq2[" << clippedBeginPosition(row(align, 1)) << ":" <<  (clippedEndPosition(row(align, 1))-1) << "]" << std::endl << std::endl;
    }
    
    
}