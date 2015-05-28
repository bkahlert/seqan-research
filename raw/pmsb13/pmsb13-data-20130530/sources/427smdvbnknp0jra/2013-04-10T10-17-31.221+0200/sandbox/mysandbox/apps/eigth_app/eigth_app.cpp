#include <iostream>
#include <seqan/align.h>

using namespace seqan;
/*
int main()
{
	typedef String<char> TSequence;
	typedef Align<TSequence,ArrayGaps> TAlign;
	typedef Row<TAlign>::Type TRow;

	TSequence seq1= "ACGTCACCTC";
	TSequence seq2= "ACGGGCCTATC";
	TAlign align;
	resize(rows(align),2);
	assignSource(row(align,0),seq1);
	assignSource(row(align,1),seq2);

	TRow &row1=row(align,0);
	TRow &row2=row(align,1);

	std::cout << align << std::endl;

	insertGap(row1,5);
	insertGaps(row1,2,2);

	insertGaps(row2,9,2);

	std::cout << align << std::endl;

	int count=0;
	typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it = begin(row1);
    TRowIterator itEnd = end(row1);
    for(; it != itEnd; ++it)
    {
        if(isGap(it))
            count++;
        
    }
    std::cout << "total count of gaps: " << count << std::endl;


}*/

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

    TAlignGraph alignG(strings);

	std::cout << alignG << std::endl; 

	return 0;

}