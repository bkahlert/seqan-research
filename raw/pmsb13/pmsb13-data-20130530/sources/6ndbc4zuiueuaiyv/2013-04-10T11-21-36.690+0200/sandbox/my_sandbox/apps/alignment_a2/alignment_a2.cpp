#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main()
{
    typedef String<Dna> TSequence;
    typedef StringSet<TSequence> TStringSet;
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;

    TSequence seq1 = "GARFIELDTHECAT";
    TSequence seq2 = "GARFIELDTHEBIGCAT";
    TSequence seq3 = "THEBIGCAT";

    TStringSet strings;
    appendValue(strings, seq1);
    appendValue(strings, seq2);
    appendValue(strings, seq3);

    //TAlignGraph alignG(strings);
    TAlignGraph globalMsaAlignment(strings, 4);
   
    ::std::cout << globalMsaAlignment << ::std::endl;
/*
    addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG),0), 0, 2),
                    addVertex(alignG, positionToId(stringSet(alignG),1), 0, 2),
                    addVertex(alignG, positionToId(stringSet(alignG),2), 0, 2));

    addVertex(alignG, positionToId(stringSet(alignG), 1), 2, 1);

    addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG),0), 2, 2),
                    addVertex(alignG, positionToId(stringSet(alignG),1), 2, 2),
                    addVertex(alignG, positionToId(stringSet(alignG),2), 3, 2));

    ::std::cout << alignG << ::std::endl;
*/
    return 0;
}

