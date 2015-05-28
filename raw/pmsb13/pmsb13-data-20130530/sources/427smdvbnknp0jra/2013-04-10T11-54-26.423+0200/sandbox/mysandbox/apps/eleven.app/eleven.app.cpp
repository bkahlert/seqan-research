#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<Rna> TSequence;
    typedef Align<TSequence, ArrayGaps> TAlign;

    TSequence seq1 = "AAGUGACUUAUUG";
    TSequence seq2 = "AGUCGGAUCUACUG";

    TAlign align;
    resize(rows(align),2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);

	
	int score = globalAlignment(align,Score<int,Simple>(1,-1,-1),MyersHirschberg());
    ::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << align << ::std::endl;

    return 0;
}