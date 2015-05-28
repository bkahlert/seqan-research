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
	typedef Align<TRnaSequence,ArrayGaps> TRAlign;      // align type
	TRnaSequence seq5 = "AAGUGACUUAUUG";
	TRnaSequence seq6 = "AGUCGGAUCUACUG";
	TRAlign align2;
	resize(rows(align2), 2);
    assignSource(row(align2,0),seq5);
    assignSource(row(align2,1),seq6);

	int score3 = globalAlignment(align2,MyersHirschberg());
    ::std::cout << "Score: " << score3 << ::std::endl;
    ::std::cout << align2 << ::std::endl;

	//Assignment 4
	typedef String<AminoAcid> TASSequence;	// sequence type
	typedef Align<TASSequence,ArrayGaps> TASAlign;      // align type
	TASSequence seq7 = "PNCFDAKQRTASRPL";
	TASSequence seq8 = "CFDKQKNNRTATRDTA";
	TASAlign align3;
	resize(rows(align3), 2);
    assignSource(row(align3,0),seq7);
    assignSource(row(align3,1),seq8);

	Score<int> scoring(3, -2, -5, -1);
    LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(scoring, 3);
	int k = 0;
    while (k<2)
    {	nextLocalAlignment(align3, enumerator);
        std::cout << "Score = " << getScore(enumerator) << std::endl;
        std::cout << align3;
		++k;
        //std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali2, 0)) << ":" << (clippedEndPosition(row(ali2, 0))-1) << "]";
        //std::cout << " and Seq2[" << clippedBeginPosition(row(ali2, 1)) << ":" <<  (clippedEndPosition(row(ali2, 1))-1) << "]" << std::endl << std::endl;
    }



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