#include <iostream>
#include <seqan/align.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<Dna> TSequence;
    typedef StringSet<TSequence> TStringSet;
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;
    TSequence seq1 = "TTGT";
    TSequence seq2 = "TTAGT";
    
    TStringSet strings;
    appendValue(strings, seq1);
    appendValue(strings, seq2);
    
    TAlignGraph alignG(strings);
    globalMsaAlignment(alignG, EditDistance);
    cout << alignG << endl;
    
}