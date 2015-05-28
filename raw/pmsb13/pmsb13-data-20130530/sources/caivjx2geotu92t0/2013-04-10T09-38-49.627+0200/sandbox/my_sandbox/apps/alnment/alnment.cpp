#include <iostream>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main() {
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
	insertGap(row2,9);
    std::cout << align;

	int count = 0;
	typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it = begin(row1);
    TRowIterator itEnd = end(row1);
    for(; it != itEnd; ++it)
    {
        if(isGap(it)) ++count;
    }
	TRowIterator it = begin(row2);
    TRowIterator itEnd = end(row2);
    for(; it != itEnd; ++it)
    {
        if(isGap(it)) ++count;
    }
    std::cout << count << std::endl;



	return 0;
}