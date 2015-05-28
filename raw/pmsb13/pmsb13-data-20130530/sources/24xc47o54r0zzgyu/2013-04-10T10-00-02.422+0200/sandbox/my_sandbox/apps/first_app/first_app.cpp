#include <iostream>
#include <seqan/align.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>

using namespace seqan;

typedef String<char> TSequence;
typedef Align<TSequence,ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;



int main()
{
   
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

	unsigned gaps;

	typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it = begin(row1);
    TRowIterator itEnd = end(row1);
	TRowIterator ita = begin(rows(align));
	TRowIterator itaEnd = end(rows(align));
	
	for (; ita != itaEnd; ++ita)
	{
		for(; it != itEnd; ++it)
		{
			if(isGap(it))
				++gaps;
		}
	}
	std::cout<<gaps;
    return 0;
}