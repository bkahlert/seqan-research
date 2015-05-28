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
    
    std::cout << align;
    TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);
    insertGaps(row1,2,2);
    insertGap(row1,5);
    insertGaps(row2,9,2);
    std::cout << align;
    
    
}