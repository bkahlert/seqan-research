#include "lagan_functions.h"
#include "stellar/stellar.h"

using namespace seqan;

int main(int argc, char const ** argv){


	seqan::CharString idOne;
	seqan::DnaString seqOne;
	seqan::CharString idTwo;
	seqan::DnaString seqTwo;

	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq1;
	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq2;

	if (!read_FASTA(argv[1], idOne, seqOne)){
		std::cerr << "error: unable to read first sequence";
		return 1;
	}
	if (!read_FASTA(argv[2], idTwo, seqTwo)){
		std::cerr << "error: unable to read second sequence";
		return 1;
	}
	unsigned minRepeatLength = 1000;
	unsigned maxRepeatPeriod = 1;
	double epsilon = 0.1;
	int minLength = 100;
	double xDrop = 5.0;
	unsigned compactThresh = 500;
	unsigned disableThresh = 50;
	unsigned numMatches = 50;

	seqan::StringSet<DnaString> setTwo;
	seqan::appendValue(setTwo, seqTwo);

	StringSet<QueryMatches<StellarMatch<DnaString, CharString> > > matches;
	resize(matches, 1);

	typedef Index<StringSet<DnaString, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
	TQGramIndex qgramIndex(setTwo);


	resize(indexShape(qgramIndex), 8);
	Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);

	indexRequire(qgramIndex, QGramSADir());

	typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
	TFinder swiftFinder(seqOne, minRepeatLength, maxRepeatPeriod);

	stellar(swiftFinder,
			swiftPattern,
			epsilon,
			minLength,
			xDrop,
			disableThresh,
			compactThresh,
			numMatches,
			0,
			idOne,
			false,
			matches,
			AllLocal());

	if (!write_seed_positions(seedsPosSeq1, seedsPosSeq2, argv[3])){
		std::cerr << "error: STELLAR output file not found";
		return 1;
	}
	typedef StellarMatch<DnaString, CharString> TMatch;
	typedef typename Iterator<String<TMatch> >::Type TIterator;
	QueryMatches<TMatch> &queryMatches = value(matches, 0);
	TIterator it = begin(queryMatches.matches);
	TIterator itEnd = end(queryMatches.matches);

	while (it != itEnd){
		std::cout << "beginPosition: " << seqan::beginPosition((*it).row1) << std::endl;
		++it;
	}
	typedef Seed<Simple> SSeed;
	typedef SeedSet<Simple> SSeedSet;

	SSeedSet seedSet;

	//writing seeds and adding them to a SeedSet
	for (unsigned i = 0; i < seedsPosSeq1.size(); ++i){
		SSeed seed(seedsPosSeq1[i].first, seedsPosSeq2[i].first, seedsPosSeq1[i].second, seedsPosSeq2[i].second);
		addSeed(seedSet, seed, Single());
	}
	//std::cout << "Trennlinie" << std::endl;
	/*typedef Iterator<SSeedSet >::Type SetIterator;

	for (SetIterator it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it){
		std::cout << *it;
		std::cout << std::endl;
	}
	typedef Iterator<String<SSeed> > StringIterator;
	*/
	seqan::String<SSeed> seedChain;

	unsigned q = 12;
	chainSeedsGlobally(seedChain, seedSet, seqan::SparseChaining());

	for(unsigned i = 0; i < length(seedChain); ++i){
		std::cout << "Initial: seedChain[" << i << "] = " << seedChain[i] << std::endl;
	}

	std::vector<std::pair<unsigned, unsigned> > startPos;
	std::vector<std::pair<unsigned, unsigned> > endPos;

	if (beginPositionH(seedChain[0]) > q && beginPositionV(seedChain[0]) > q){
		startPos.push_back(std::make_pair(0, 0));
		endPos.push_back(std::make_pair(beginPositionH(seedChain[0]), beginPositionV(seedChain[0])));
	}
	for (unsigned i = 0; i < length(seedChain) - 1; ++i){
		unsigned startSeq1 = endPositionH(seedChain[i]);
		unsigned startSeq2 = endPositionV(seedChain[i]);
		unsigned endSeq1 = beginPositionH(seedChain[i + 1]);
		unsigned endSeq2 = beginPositionV(seedChain[i + 1]);

		startPos.push_back(std::make_pair(startSeq1, startSeq2));
		endPos.push_back(std::make_pair(endSeq1, endSeq2));
	}
	if (endPositionH(back(seedChain)) < length(seedChain) - q && endPositionV(back(seedChain)) < length(seedChain) - q){
		startPos.push_back(std::make_pair(endPositionH(back(seedChain)), endPositionV(back(seedChain))));
		endPos.push_back(std::make_pair(length(seqOne), length(seqTwo)));
	}

	if(!insert_seeds_in_initial_seedChain(seedChain, startPos, endPos, seqOne, seqTwo, q)){
		std::cerr << "error: couldn't insert new found seeds!";
		return 1;
	}
	//erase(seedChain, 8);
	//erase(seedChain, 16);
	for(unsigned i = 0; i < length(seedChain); ++i){
		std::cout << "seedChain[" << i << "] = " << seedChain[i] << std::endl;
	}

	Align<DnaString, ArrayGaps> alignment;
	resize(seqan::rows(alignment), 2);
	Score<int, Simple> scoringFunction(2, -1, -2);
	seqan::assignSource(seqan::row(alignment, 0), seqOne);
	seqan::assignSource(seqan::row(alignment, 1), seqTwo);

	int result = bandedChainAlignment(alignment, seedChain, scoringFunction, 2);
	const char* filename = "lagan.txt";

	output_alignment(filename, idOne, idTwo, alignment, result);
	return 0;
}
