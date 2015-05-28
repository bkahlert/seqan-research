#include <iostream>
#include <seqan/align.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>

using namespace seqan;

typedef String<Rna> TSequence;
typedef StringSet<TSequence> TStringSet;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow; 

int main()
{
	TSequence seq1 = "AAGUGACUUAUUG ";
	TSequence seq2 = "AGUCGGAUCUACUG ";
	TAlign align;
    resize(rows(align),2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);

	int score = globalAlignment(align,MyersHirschberg());
    ::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << align << ::std::endl;


	TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);

	typedef Iterator<TRow>::Type TRowIterator;
    
	for (unsigned j=0;j<2;++j)
	{
		TRowIterator it = begin(row(align,j));
		TRowIterator itEnd = end(row(align,j));
		for(unsigned i=0; it != itEnd; ++it)
		{
			if(isGap(it))
				std::cout << "Gap in seq."<<j+1<<" at pos." << i <<"."<<std::endl;
			++i;
		}
	}
	std::cout << std::endl;

	return 0;
}