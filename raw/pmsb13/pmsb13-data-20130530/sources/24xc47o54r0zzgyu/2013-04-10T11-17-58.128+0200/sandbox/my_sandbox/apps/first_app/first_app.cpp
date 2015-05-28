#include <iostream>
#include <seqan/align.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>

using namespace seqan;

typedef String<Dna> TSequence;
typedef StringSet<TSequence> TStringSet;
typedef StringSet<TSequence, Dependent<> > TDepStringSet;  
typedef Graph<Alignment<TDepStringSet> > TAlignGraph;

int main()
{
	TSequence seq1 = "AAATGACGGATTG";
	TSequence seq2 = "TGGGA";
	TStringSet sequences;
    appendValue(sequences,seq1);
    appendValue(sequences,seq2);

	TAlignGraph alignG(sequences);

	int score = globalAlignment(alignG, Score<int,EditDistance>(), AlignConfig<false, true, false, true>());
    ::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << alignG << ::std::endl;

	return 0;
}