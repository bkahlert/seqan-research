#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
    typedef String<char> TSequence;                 // sequence type
    typedef StringSet<TSequence> TStringSet;       // container for strings
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;   // dependent string set
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;   // alignment graph
TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "TGGGACCCCCCCCCC";

    TStringSet sequences;
    appendValue(sequences,seq1);
    appendValue(sequences,seq2);

    TAlignGraph alignG(sequences);
    int score=-1;
    if(length(seq1)<=length(seq2)) {
    score = globalAlignment(alignG, Score<int,Simple>(1,-1,-1), AlignConfig<true,false,false,true>());
    } else {
    score = globalAlignment(alignG, Score<int,Simple>(1,-1,-1), AlignConfig<false,true,true,false>()); 
    }
    ::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << alignG << ::std::endl;

    return 0;
}

