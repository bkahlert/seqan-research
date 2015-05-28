#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
    
    TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "AGTCGGATCTACTG";
    
    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);
    
    int score=globalAlignment(align,Score<int, Simple>(4,-2,-2,-4));
    cout << align << endl;
    cout << score << endl;
    
}