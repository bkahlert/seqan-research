#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{

	typedef String<Dna> TSequence;                 
    typedef StringSet<TSequence> TStringSet;
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;
	typedef Graph<Alignment<TDepStringSet> > TAlignGraph;

	TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "TGGGA";

    TStringSet sequences;
    appendValue(sequences,seq1);
    appendValue(sequences,seq2);

    TAlignGraph alignG(sequences);

	int score = globalAlignment(alignG, Score<int,Simple>(0,-1,-1), AlignConfig<true, false, false, true>());
    ::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << alignG << ::std::endl;

    return 0;
}