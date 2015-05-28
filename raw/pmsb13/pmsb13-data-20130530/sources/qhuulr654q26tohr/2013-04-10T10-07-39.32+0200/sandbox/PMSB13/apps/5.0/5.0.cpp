#include <iostream>
#include <seqan/align.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
    
    TSequence seq1 = "ACGTCACCTC";
    TSequence seq2 = "ACGGGCCTATC";
    
    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);
    
    std::cout << align;
    TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);
    insertGaps(row1,2,2);
    insertGap(row1,5);
    insertGaps(row2,9,2);
    std::cout << align;
    
    typedef Iterator<TRow>::Type TIt1;
    typedef Iterator<Rows<TAlign>::Type>::Type TIt2;
    TIt1 it = begin(row(align,0));
    
    /*while(!atEnd(it))
	goNext(it);
    
    
    cout << countGaps(it) << endl;*/
    
}