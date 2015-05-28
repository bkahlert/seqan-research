#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
	typedef String<char> TSequence;
	typedef Align<TSequence,ArrayGaps> TAlign;
	typedef Row<TAlign>::Type TRow;

	TSequence seq1= "ACGTCACCTC";
	TSequence seq2= "ACGGGCCTATC";
	TAlign align;
	resize(rows(align),2);
	assignSource(row(align,0),seq1);
	assignSource(row(align,1),seq2);

	TRow row1=row(align,0);
	TRow row2=row(align,1);

	std::cout << align << std::endl;

	insertGaps(row1,2,3);

	std::cout << align << std::endl;

}