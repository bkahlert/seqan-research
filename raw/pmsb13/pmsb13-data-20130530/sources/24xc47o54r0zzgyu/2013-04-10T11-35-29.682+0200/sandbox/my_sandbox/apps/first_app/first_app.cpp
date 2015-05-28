#include <iostream>
#include <seqan/align.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>

using namespace seqan;

typedef String<Rna> TSequence;
typedef StringSet<TSequence> TStringSet;
typedef Align<TSequence, ArrayGaps> TAlign;

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

	return 0;
}