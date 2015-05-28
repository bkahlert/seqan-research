#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type

    TSequence seq1 = "ACGTCACCTC";
    TSequence seq2 = "ACGGGCCTATC";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);
    std::cout << align;
    TRow &row1 = row(align, 0);
    TRow &row2 = row(align, 1);
    insertGaps(row1, 2, 2);
    insertGap(row1, 7);
    insertGaps(row2, 9, 4);
    std::cout << align;
    typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it1 = begin(row1);
    TRowIterator it1End = end(row1);
    TRowIterator it2 = begin(row2);
    TRowIterator it2End = end(row2);
    unsigned noOfGaps = 0;
    for(; it1 != it1End && it2 != it2End; ++it1, ++it2)
    {
        if(isGap(it1)&&isGap(it2))
            noOfGaps+=2;
        if(isGap(it1)||isGap(it2))
            ++noOfGaps;

    }
    std::cout << noOfGaps << std::endl;
}
