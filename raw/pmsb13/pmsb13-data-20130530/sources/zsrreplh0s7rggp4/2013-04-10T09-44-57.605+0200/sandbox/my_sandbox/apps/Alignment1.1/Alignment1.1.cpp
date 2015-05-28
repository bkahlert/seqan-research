#include <iostream>
#include <seqan/align.h>

using namespace seqan;

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
    insertGap(row1,5);
    insertGaps(row1,2,2);
    insertGaps(row2,9,2);
    std::cout << align;

    // Initialize the row iterators.
    TRowIterator itRow1 = begin(row1);
    TRowIterator itEndRow1 = end(row1);
    TRowIterator itRow2 = begin(row2);
    // Iterate over both rows simultaneously.
    int gapCount = 0;
    for(;itRow1 != itEndRow1; ++itRow1, ++itRow2)
    {
        if(isGap(itRow1))
        {
            gapCount += countGaps(itRow1);  // Count the number of consecutive gaps from the current position in row1.
            itRow1 += countGaps(itRow1);    // Jump to next position to check for gaps.
            itRow2 += countGaps(itRow1);    // Jump to next position to check for gaps.
        }
        if(isGap(itRow2))
        {
            gapCount += countGaps(itRow2);  // Count the number of consecutive gaps from the current position in row2.
            itRow1 += countGaps(itRow2);    // Jump to next position to check for gaps.
            itRow2 += countGaps(itRow2);    // Jump to next position to check for gaps.
        }
    }
    // Print the result.
    std::cout << "Number of gaps: " << gapCount << std::endl;
}