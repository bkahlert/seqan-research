#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/align.h>
#include <seqan/seeds2.h>

using namespace seqan;

// set parameters

unsigned const minimalLength = 12;
float const errorRate = 0.05;
//unsigned const qGramLength = 8;
CharString const filename = "ministellar.out";

typedef int TScoreValue;
TScoreValue const matchScore = 1;
TScoreValue const mismatchScore = -2;
TScoreValue const gapScore = -2;
TScoreValue const minScore = 9;
TScoreValue const xDrop = 3;


int main(int argc, char const ** argv)
{
	if (argc != 3)
	{
		std::cerr << "ERROR: Invalid argument count!" << std::endl;
		std::cerr << "USAGE: ministallar IN1.fa IN2.fa" << std::endl;
		return 1;
	}

    Score<int, Simple> scoring = Score<int, Simple>(matchScore, mismatchScore, gapScore);
	// DONE: Define a scoring scheme 'score' with linear gap costs using the above defined
	//       score values for scoring matches, mismatches, and gaps.
	// HINT: You can find a section on Schoring Schemes in the Alignments tutorial.

	// read sequences from fasta files
	typedef Dna5String TSequence;
	StringSet<CharString> ids1, ids2;
	StringSet<TSequence> seqs1, seqs2;

	std::ifstream streamIn1(argv[1], std::ios::in | std::ios::binary);
	RecordReader<std::ifstream, SinglePass<> > reader1(streamIn1);
	if (read2(ids1, seqs1, reader1, Fasta()) != 0)
	{
		std::cerr << "Could not read file " << argv[1] << std::endl;
		return 1;
	}
	// DONE: Do the same for the second file.

    std::ifstream streamIn2(argv[2], std::ios::in | std::ios::binary);
    RecordReader<std::ifstream, SinglePass<> > reader2(streamIn2);
	if (read2(ids2, seqs2, reader2, Fasta()) != 0)
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

    typedef Index<TSequence, IndexQGram<UngappedShape<8> > > TQIndex;
    //typedef Finder<TQIndex, Swift<SwiftLocal> > TFinder;
    //typedef Pattern<TQIndex, Swift<SwiftLocal> > TPattern;

	// TODO: Define variables of your pattern and finder types, and initialize them with the
	//       first sequences in the sequence sets seqs1 and seqs2.
    //TFinder finder = TFinder(seqs1[0]);
    //TPattern pattern = TPattern(seqs2[0]);
	// TODO: Repeat the following steps for all hits of the finder.

	// call the function find to obtain a swift hit
	// TODO: uncomment the following line

	//find(finder, pattern, errorRate, minimalLength);

	// extract infixes from hit
	// TODO: uncomment the following three lines

	//typedef Infix<TSequence>::Type TInfix;
	//TInfix infix1 = infix(finder, seqs1[0]);
	//TInfix infix2 = infix(pattern, seqs2[0]);

	// find best local alignment of infixes, and check for minimal score
	// TODO: Define an align object and initialize it with the infixes. Find the best local
	//       alignment of the infixes. Use the above defined scoring scheme. Make sure that
	//       it reaches the minimal score minScore.

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