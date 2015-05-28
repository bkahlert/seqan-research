#include <iostream>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main() {
	typedef String<char> TSequence;                 // sequence type
    typedef Align<TSequence,ArrayGaps> TAlign;      // align type
    typedef Row<TAlign>::Type TRow;                 // gapped sequence type
    typedef StringSet<TSequence> TStringSet;       // container for strings
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;   // dependent string set
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;   // alignment graph

	TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "AGTCGGATCTACTG";
	 
    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);

	int score = globalAlignment(align, Score<int,Simple>(4,-2,-2,-4));
	::std::cout << "Score: " << score << ::std::endl;
    ::std::cout << align << ::std::endl;

	// Assignment 2
	TSequence seq3 = "AAATGACGGATTG";
	TSequence seq4 = "TGGGA";

	TStringSet sequences;
    appendValue(sequences,seq3);
    appendValue(sequences,seq4);

    TAlignGraph alignG(sequences);

	int score2 = globalAlignment(alignG, Score<int,Simple>(1,-1,-1), AlignConfig<true, false, false, true>());
    ::std::cout << "Score: " << score2 << ::std::endl;
    ::std::cout << alignG << ::std::endl;

	// Assignment 3
	typedef String<Rna> TRnaSequence;	// sequence type
	TRnaSequence seq5 = "AAGUGACUUAUUG";
	TRnaSequence seq6 = "AGUCGGAUCUACUG";
	TStringSet sequences2;
	appendValue(sequences2,seq5);
	appendValue(sequences2,seq6);

	int score3 = globalAlignmentScore(sequences2,MyersBitVector());
    ::std::cout << "Score: " << score3 << ::std::endl;
    ::std::cout << sequences2 << ::std::endl;


	/*
	std::cout << align;
    TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);
    insertGaps(row1,2,2);
    insertGap(row1,5);
	insertGaps(row2,9,2);
    std::cout << align;

	int count = 0;
	typedef Iterator<TRow>::Type TRowIterator;
    TRowIterator it = begin(row1);
	TRowIterator it2 = begin(row2);
    TRowIterator itEnd = end(row1);
    for(; it != itEnd; ++it, ++it2)
    {
        if(isGap(it)) ++count;
		if(isGap(it2)) ++count;
    }
	
    TRowIterator it2End = end(row2);
    for(; it2 != it2End; ++it2)
    {
        if(isGap(it2)) ++count;
    }
    std::cout << count << std::endl;
	int count2 = countGaps(it);
	std::cout << count2 << std::endl;*/



	return 0;
}