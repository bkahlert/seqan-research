#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/align.h>
#include <seqan/seeds2.h>

using namespace seqan;

// set parameters

unsigned const minimalLength = 12;
float const errorRate = 0.05;
//unsigned const qGramLength = 8;
CharString const filename = "exastellar.out";

typedef int TScoreValue;
TScoreValue const matchScore = 1;
TScoreValue const mismatchScore = -2;
TScoreValue const gapScore = -2;
TScoreValue const minScore = 9;
TScoreValue const xDrop = 3;

template<typename TSeed, typename TSeq>
void
writeSeed(TSeed & seed, TSeq const & seq0, TSeq const & seq1) {
    std::cout << "Seed from position " << leftPosition(seed, 0);
    std::cout << " to " << rightPosition(seed, 0) << ": ";
    std::cout << infix(seq0, leftPosition(seed, 0), rightPosition(seed, 0)+1) << std::endl;
    std::cout << "Seed from position " << leftPosition(seed, 1);
    std::cout << " to " << rightPosition(seed, 1) << ": ";
    std::cout << infix(seq1, leftPosition(seed, 1), rightPosition(seed, 1)+1) << std::endl;
}

int main(int argc, char const ** argv)
{
	if (argc != 3)
	{
		std::cerr << "ERROR: Invalid argument count!" << std::endl;
		std::cerr << "USAGE: ministallar IN1.fa IN2.fa" << std::endl;
		return 1;
	}
	
	Score< TScoreValue > judgement_of_carrion( matchScore, mismatchScore, gapScore );
	
	// read sequences from fasta files
	typedef Dna5String TSequence;
	StringSet<CharString> id1, id2;
	StringSet<TSequence> seq1, seq2;

	std::ifstream streamInPrimus(argv[1], std::ios::in | std::ios::binary);
	RecordReader<std::ifstream, SinglePass<> > readerUno(streamInPrimus);
	if (read2(id1, seq1, readerUno, Fasta()) != 0)
	{
		std::cerr << "Could not read file " << argv[1] << std::endl;
		return 1;
	}
    
    std::ifstream streamInDeux(argv[2], std::ios::in | std::ios::binary);
	RecordReader<std::ifstream, SinglePass<> > readerDeux(streamInDeux);
	if (read2(id2, seq2, readerDeux, Fasta()) != 0)
	{
		std::cerr << "Could not read file " << argv[2] << std::endl;
		return 1;
	}

	// prepare output file
	std::ofstream outFile;
	outFile.open(toCString(filename));

	// define finder and pattern
	// TODO: Define types for finder and pattern using SeqAn's find interface from the index
	//       module. Both should be specialized for the approproiate swift filter algorithm.
	//       The pattern will need the definition of a q-gram-index. Use q-grams of length 8
	//       for now, and open addressing.
	// HINT: The tutorial on Pattern Matching describes the find interface. At the end of that
	//       tutorial page, you will find a link to a Swift HowTo page.
	
	typedef Index<TSequence, IndexQGram<Shape<Dna5, UngappedShape<8> >, OpenAddressing> > TIndex;
	
	TIndex mt_doom( seq2[1] );
	
	Finder< TSequence, Swift< SwiftLocal > > gollum( seq1[0] );
	Pattern< TIndex, Swift< SwiftLocal > > the_ring( mt_doom );

	// TODO: Repeat the following steps for all hits of the finder.

	// call the function find to obtain a swift hit
	// TODO: uncomment the following line

	while( find(gollum, the_ring, errorRate, minimalLength) ){
	
	    typedef Infix<TSequence>::Type TInfix;
	    TInfix infix1 = infix(gollum, seq1[0]);
	    TInfix infix2 = infix(the_ring, seq2[0]);
	    
	    Align< TSequence > the_showdown;
	    appendValue( rows( the_showdown ), infix1 );
	    appendValue( rows( the_showdown ), infix2 );
	    
	    TScoreValue the_score = localAlignment( the_showdown, judgement_of_carrion, minScore );
    	if( the_score < minScore ){
	        std::cerr << "DerTooooooooooooD!" << std::endl;
	        return -1;
	    }
	    
	    Seed< Simple > wind( clippedBeginPosition( row( the_showdown, 0 ) ), clippedBeginPosition( row( the_showdown, 1 ) ), clippedEndPosition( row( the_showdown, 0 ) ), clippedEndPosition( row( the_showdown, 1 ) ) );
	    writeSeed( wind, seq1[0], seq2[0] );	    	    
	    extendSeed( wind, seq1[0], seq2[0], 2, xDrop, judgement_of_carrion, GappedXDrop() );
	    std::cout << "After extension:" << std::endl;
	    writeSeed( wind, seq1[0], seq2[0] );
	    std::cout << "---" << std::endl;
	}	

	// create a seed for the local alignment, and conduct gapped X-drop extension
	// TODO: Define a seed on the original sequences but for the subsequences of the local
	//       alignment. Extend this seed in both directions using gapped X-drop extension.
	//       Check that the extended seed has a length of at least minimalLength.
	// HINT: Have a look at the Seed-and-Extend tutorial.

	// find best global alignment of extended seed
	// TODO: Compute the best global alignment of extended seed and its score.
	// HINT: Create an align object on infixes of the sequences.

	// TODO: Output the alignment as a match to the output file.

	outFile.close();

	return 0;
}