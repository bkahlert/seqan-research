#include <iostream>
#include <seqan/align.h>

using namespace seqan;
int main()
{
	typedef String<char>            TSequence;  // sequence type
    typedef Align<TSequence,ArrayGaps>  TAlign;     // align type

    TSequence seq1 = "blablablu";
    TSequence seq2 = "abab";
	
	StringSet<TSequence> sequences;
    appendValue(sequences,seq1);
    appendValue(sequences,seq2);
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;
    TAlignGraph alignG(sequences);

	AlignConfig<true,true,true,true> ac;
    int score = globalAlignment(alignG, Score<int>(1,-1,-1,-1), ac, Gotoh());
    ::std::cout << "Score = " << score << ::std::endl;
    ::std::cout << alignG;

	return 0;
}