#include <iostream>
#include <seqan/align.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
    
    TSequence seq1 = "CDFGDC";
    TSequence seq2 = "CDEFGAHGC";
    
    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);
    
    std::cout << align;
    TRow &row1 = rows(align)[0];
    TRow &row2 = row(align,1);
    insertGap(row1,2);
    insertGaps(row1,5,2);
    std::cout << align;
    std::cout << length(row(align,0)) << std::endl;
}