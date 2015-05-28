#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;
using namespace std;

int main()
{
    typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
    
    typedef StringSet<TSequence> TStringSet;
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;
    
    TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "TGGGA";
    
    TStringSet strings;
    appendValue(strings, seq1);
    appendValue(strings, seq2);
    appendValue(strings, seq3);
    
    TAlignGraph alignG(strings);
    
    int score=globalAlignment(alignG,Score<int, Simple>(0,-1,-1),AlignConfig<true,false,false,true);
    cout << alignG << endl;
    cout << score << endl;
    
}