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

	// TODO: Define a scoring scheme 'score' with linear gap costs using the above defined
	//       score values for scoring matches, mismatches, and gaps.
	// HINT: You can find a section on Schoring Schemes in the Alignments tutorial.
	Score<TScoreValue, Simple> score;
	setScoreGap(score, gapScore);
	setScoreMismatch(score, mismatchScore);
	setScoreMatch(score, matchScore);

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
	streamIn1.close();

	// TODO: Do the same for the second file.

	std::ifstream streamIn2(argv[2], std::ios::in | std::ios::binary);
	RecordReader<std::ifstream, SinglePass<> > reader2(streamIn2);
	if (read2(ids2, seqs2, reader2, Fasta()) != 0)
	{
		std::cerr << "Could not read file " << argv[2] << std::endl;
		return 1;
	}
	streamIn2.close();

	// prepare output file
	std::ofstream outFile;
	outFile.open(toCString(filename));

	// Print some info output
	std::cerr << "read " << length(ids1) << " reads from " << argv[1] << std::endl;
	std::cerr << "read " << length(ids2) << " reads from " << argv[2] << std::endl;

	// define finder and pattern
	// TODO: Define types for finder and pattern using SeqAn's find interface from the index
	//       module. Both should be specialized for the approproiate swift filter algorithm.
	//       The pattern will need the definition of a q-gram-index. Use q-grams of length 8
	//       for now, and open addressing.
	// HINT: The tutorial on Pattern Matching describes the find interface. At the end of that
	//       tutorial page, you will find a link to a Swift HowTo page.

	typedef Index<TSequence, IndexQGram<Shape<Dna, UngappedShape<8> >, OpenAddressing> > TIndex;
	typedef Pattern<TIndex, Swift<SwiftLocal> > TPattern;
	typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;

	// TODO: Define variables of your pattern and finder types, and initialize them with the
	//       first sequences in the sequence sets seqs1 and seqs2.

	TIndex index(seqs2[0]);
	TPattern pattern(index);

	TFinder finder(seqs1[0]);

	// TODO: Repeat the following steps for all hits of the finder.

	// call the function find to obtain a swift hit
	// TODO: uncomment the following line

	while (find(finder, pattern, errorRate, minimalLength)) {
	    // extract infixes from hit
	    // TODO: uncomment the following three lines

	    typedef Infix<TSequence>::Type TInfix;
	    TInfix infix1 = infix(finder, seqs1[0]);
	    TInfix infix2 = infix(pattern, seqs2[0]);

	    std::cerr << "infix 1 " << infix1 << std::endl;
	    std::cerr << "infix 2 " << infix2 << std::endl;

	    // find best local alignment of infixes, and check for minimal score
	    // TODO: Define an align object and initialize it with the infixes. Find the best local
	    //       alignment of the infixes. Use the above defined scoring scheme. Make sure that
	    //       it reaches the minimal score minScore.

	    Align<String<Dna> > align;
	    appendValue(rows(align), infix1);
	    appendValue(rows(align), infix2);

	    TScoreValue score = localAlignment(align, score);
	    if (score < minScore)
		continue;


	    // create a seed for the local alignment, and conduct gapped X-drop extension
	    // TODO: Define a seed on the original sequences but for the subsequences of the local
	    //       alignment. Extend this seed in both directions using gapped X-drop extension.
	    //       Check that the extended seed has a length of at least minimalLength.
	    // HINT: Have a look at the Seed-and-Extend tutorial.

	    std::cerr << "beginPosition() " <<  beginPosition(infix1) << std::endl;
	    Seed<Simple> seed(beginPosition(infix1)+clippedBeginPosition(row(align,0)), 
		    beginPosition(infix2)+clippedBeginPosition(row(align,1)), 
		    beginPosition(infix1)+clippedEndPosition(row(align,0)),
		    beginPosition(infix2)+clippedEndPosition(row(align,1)));

	    std::cerr << "seed before extension: [" << getBeginDim0(seed) << "," << getEndDim0(seed) << "] [" 
		<< getBeginDim1(seed) << "," << getEndDim1(seed) << "]\n";

	    extendSeed(seed, seqs2[0], seqs1[0], EXTEND_BOTH, score, xDrop, GappedXDrop());
	    
	    std::cerr << "seed after extension: [" << getBeginDim0(seed) << "," << getEndDim0(seed) << "] [" 
		<< getBeginDim1(seed) << "," << getEndDim1(seed) << "]\n";

	    // Check if the match is long enough after it was extended
	    if (getEndDim1(seed)-getBeginDim1(seed) < minimalLength)
	       continue;	

	    // find best global alignment of extended seed
	    // TODO: Compute the best global alignment of extended seed and its score.
	    // HINT: Create an align object on infixes of the sequences.

	    TInfix extendedInfix1 = infix(seqs1[0], getBeginDim0(seed), getEndDim0(seed));
	    TInfix extendedInfix2 = infix(seqs2[0], getBeginDim1(seed), getEndDim1(seed));
	
	    Align<String<Dna> > globAlign;
	    appendValue(rows(globAlign), extendedInfix1);
	    appendValue(rows(globAlign), extendedInfix2);

	    globalAlignment(globAlign, score);

	    std::cerr << globAlign << std::endl;



	    // TODO: Output the alignment as a match to the output file.

	}

	outFile.close();

	return 0;
}