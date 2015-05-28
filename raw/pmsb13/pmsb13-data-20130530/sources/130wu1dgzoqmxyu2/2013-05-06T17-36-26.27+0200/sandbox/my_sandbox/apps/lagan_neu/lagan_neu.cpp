#include "lagan_functions.h"
#include "stellar/stellar.h"
#include <seqan/arg_parse.h>

using namespace seqan;

struct LaganOptions
{
	//LAGAN parameters
	char * seqOne;
	char * seqTwo;
	unsigned q;
	char * outputFile;

	//stellar parameters
	unsigned stellarKmer;

	unsigned repeatLength;
	unsigned repeatPeriod;
	double epsilon;
	int minLength;
	double xDrop;
	unsigned sortThresh;
	unsigned disableThresh;
	unsigned numMatches;
	//bool databaseStrand; forward and reverse search of the database
	/*LaganOptions() {
			outputFile = "laganBlub.txt";

			stellarKmer = seqan::maxValue<unsigned>();
			epsilon = 0.05;
			minLength = 100;
			xDrop = 5;

			disableThresh = seqan::maxValue<unsigned>();
			numMatches = 50;
			repeatPeriod = 1;
			repeatLength = 1000;
		} */

};

seqan::ArgumentParser::ParseResult
parseCommandLine(LaganOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("LAGAN");
    setShortDescription(parser, "Reimplementation of LAGAN (lite version) , a global aligner");
	setDate(parser, __DATE__);
	setVersion(parser, "0.5 beta");
	setCategory(parser, "Global Alignment");

	addUsageLine(parser, "<\\fIFASTA FILE 1\\fP> <\\fIFASTA FILE 2\\fP> [\\fIOPTIONS\\fP]");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA FILE 1"));
	setValidValues(parser, 0, "fa fasta");  // allow only fasta files as input
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA FILE 2"));
	setValidValues(parser, 1, "fa fasta");  // allow only fasta files as input

    // Define Options
	addSection(parser, " LAGAN options");

    addOption(parser, seqan::ArgParseOption(
    		"q", "qGram", "qGram length to use for the index.",
    		seqan::ArgParseArgument::INTEGER, "UNSIGNED"));
    setDefaultValue(parser, "q", "12");

    addOption(parser, seqan::ArgParseOption(
    		"o", "ouputFilename", "name of the output file.",
    		seqan::ArgParseArgument::OUTPUTFILE));
    setDefaultValue(parser, "o", "laganBlub.txt");

    addSection(parser, " Stellar: Main Options");
    addOption(parser, ArgParseOption("e", "epsilon", "Maximal error rate (max 0.25).", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "e", "0.05");
    setMinValue(parser, "e", "0.0000001");
    setMaxValue(parser, "e", "0.25");
    addOption(parser, ArgParseOption("l", "minLength", "Minimal length of epsilon-matches.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "l", "100");
    setMinValue(parser, "l", "0");

    addSection(parser, " Stellar: Filtering Options");

	addOption(parser, ArgParseOption("k", "stellarKmer", "Length of the q-grams in stellar (max 32).", ArgParseArgument::INTEGER));
	setMinValue(parser, "k", "1");
	setMaxValue(parser, "k", "32");
	setDefaultValue(parser, "k", "8");
	addOption(parser, ArgParseOption("rp", "repeatPeriod",
									 "Maximal period of low complexity repeats to be filtered.", ArgParseArgument::INTEGER));
	setDefaultValue(parser, "rp", "1");
	addOption(parser, ArgParseOption("rl", "repeatLength",
									 "Minimal length of low complexity repeats to be filtered.", ArgParseArgument::INTEGER));
	setDefaultValue(parser, "rl", "1000");

	addSection(parser, "Verification Options");

	addOption(parser, ArgParseOption("x", "xDrop", "Maximal x-drop for extension.", ArgParseArgument::DOUBLE));
	setDefaultValue(parser, "x", "5");
	addOption(parser, ArgParseOption("dt", "disableThresh",
	                                     "Maximal number of verified matches before disabling verification for one query "
	                                     "sequence (default infinity).", ArgParseArgument::INTEGER));
	setMinValue(parser, "dt", "0");
	setDefaultValue(parser, "dt", "50");
	addOption(parser, ArgParseOption("n", "numMatches",
	                                     "Maximal number of kept matches per query and database. If STELLAR finds more matches, "
	                                     "only the longest ones are kept.", ArgParseArgument::INTEGER));
	setDefaultValue(parser, "n", "50");
	addOption(parser, ArgParseOption("s", "sortThresh",
									 "Number of matches triggering removal of duplicates. Choose a smaller value for saving "
									 "space.", ArgParseArgument::INTEGER));
	setDefaultValue(parser, "s", "500");

	std::cout << "ficken" << std::endl;
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    //Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
    	return res;

    // Extract option values.
    std::cout << "blasen" << std::endl;
    seqan::getArgumentValue(options.seqOne, parser, 0);
    seqan::getArgumentValue(options.seqTwo, parser, 1);
    getOptionValue(options.q, parser, "q");
    getOptionValue(options.outputFile, parser, "o");
    getOptionValue(options.stellarKmer, parser, "k");
    getOptionValue(options.repeatLength, parser, "rl");
    getOptionValue(options.repeatPeriod, parser, "rp");
    getOptionValue(options.epsilon, parser, "e");
    getOptionValue(options.minLength, parser, "l");
    getOptionValue(options.xDrop, parser, "x");
    getOptionValue(options.disableThresh, parser, "dt");
    getOptionValue(options.numMatches, parser, "n");
    getOptionValue(options.sortThresh, parser, "s");
    std::cout << "69" << std::endl;
    //options.toUppercase = isSet(parser, "uppercase");
    //options.toLowercase = isSet(parser, "lowercase");

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv){
	/*
	//LAGAN parameters
	unsigned q = 12;
	const char* outputFilename = "lagan.txt";

	//stellar parameters
	unsigned stellarKmer = 8;

	unsigned minRepeatLength = 1000;
	unsigned maxRepeatPeriod = 1;
	double epsilon = 0.1;
	int minLength = 50;
	double xDrop = 5.0;
	unsigned compactThresh = 500;
	unsigned disableThresh = 50;
	unsigned numMatches = 50;
	bool databaseStrand = false; //forward and reverse search of the database
	*/
	// Parse the command line.
	std::cout << "blastaa" << std::endl;
	LaganOptions options;
	std::cout << "2. teil" << std::endl;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
	std::cout << "3. teil" << std::endl;
	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		std::cout << "blub" << std::endl;
		return res == seqan::ArgumentParser::PARSE_ERROR;
	}
	std::cout << options.q << std::endl;
	seqan::CharString idOne;
	seqan::DnaString seqOne;
	seqan::CharString idTwo;
	seqan::DnaString seqTwo;

	typedef Seed<Simple> SSeed;
	typedef SeedSet<Simple> SSeedSet;

	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq1;
	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq2;

	//read-in of FASTA sequences
	if (!read_FASTA(options.seqOne, idOne, seqOne)){
		std::cerr << "error: unable to read first sequence";
		return 1;
	}
	if (!read_FASTA(options.seqTwo, idTwo, seqTwo)){
		std::cerr << "error: unable to read second sequence";
		return 1;
	}

	//stellar call to obtain initial seeds
	seqan::StringSet<DnaString> setTwo;
	seqan::appendValue(setTwo, seqTwo);

	StringSet<QueryMatches<StellarMatch<DnaString, CharString> > > matches;
	resize(matches, 1);

	typedef Index<StringSet<DnaString, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
	TQGramIndex qgramIndex(setTwo);
	resize(indexShape(qgramIndex), options.stellarKmer);
	Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);

	indexRequire(qgramIndex, QGramSADir());

	typedef Finder<DnaString, Swift<SwiftLocal> > TFinder;
	TFinder swiftFinder(seqOne, options.repeatLength, options.repeatPeriod);

	stellar(swiftFinder,
			swiftPattern,
			options.epsilon,
			options.minLength,
			options.xDrop,
			options.disableThresh,
			options.sortThresh,
			options.numMatches,
			0, //argument for forward/reverse search
			idOne,
			false,
			matches,
			AllLocal());

	/*if (!write_seed_positions(seedsPosSeq1, seedsPosSeq2, argv[3])){
		std::cerr << "error: STELLAR output file not found";
		return 1;
	}
	*/
	//extracting positions from stellar matches (seeds) and compute best chain
	typedef StellarMatch<DnaString, CharString> TMatch;
	typedef typename Iterator<String<TMatch> >::Type TIterator;
	QueryMatches<TMatch> &queryMatches = value(matches, 0);
	TIterator it = begin(queryMatches.matches);
	TIterator itEnd = end(queryMatches.matches);

	int z = 0;
	SSeedSet seedSet;

	while (it != itEnd){
		std::cout << "beginPositionDatabase[" << z << "]: " << seqan::beginPosition((*it).row1) << std::endl;
		std::cout << "endPositionDatabase[" << z << "]: " << seqan::endPosition((*it).row1) << std::endl;
		std::cout << "beginPositionQuery[" << z << "]: " << seqan::beginPosition((*it).row2) << std::endl;
		std::cout << "endPositionQuery[" << z << "]: " << seqan::endPosition((*it).row1) << std::endl;

		SSeed stellarSeed(seqan::beginPosition((*it).row1),
						seqan::beginPosition((*it).row2),
						seqan::endPosition((*it).row1),
						seqan::endPosition((*it).row2));

		addSeed(seedSet, stellarSeed, seqan::Single());
		++it;
 		++z;
	}

	seqan::String<SSeed> seedChain;
	seqan::chainSeedsGlobally(seedChain, seedSet, seqan::SparseChaining());



	//writing seeds and adding them to a SeedSet
	/*for (unsigned i = 0; i < seedsPosSeq1.size(); ++i){
		SSeed seed(seedsPosSeq1[i].first,
				seedsPosSeq2[i].first,
				seedsPosSeq1[i].second,
				seedsPosSeq2[i].second);

		addSeed(seedSet, seed, Single());
	}*/
	//std::cout << "Trennlinie" << std::endl;
	/*typedef Iterator<SSeedSet >::Type SetIterator;

	for (SetIterator it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it){
		std::cout << *it;
		std::cout << std::endl;
	}
	typedef Iterator<String<SSeed> > StringIterator;
	*/
	//seqan::String<SSeed> seedChain;

	//chainSeedsGlobally(seedChain, seedSet, seqan::SparseChaining());

	for(unsigned i = 0; i < length(seedChain); ++i){
		std::cout << "Initial: seedChain[" << i << "] = " << seedChain[i] << std::endl;
	}

	std::vector<std::pair<unsigned, unsigned> > startPos;
	std::vector<std::pair<unsigned, unsigned> > endPos;

	//creating a container to define the gaps between the initial stellar matches
	if (beginPositionH(seedChain[0]) > options.q && beginPositionV(seedChain[0]) > options.q){
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
	if (endPositionH(back(seedChain)) < length(seedChain) - options.q && endPositionV(back(seedChain)) < length(seedChain) - options.q){
		startPos.push_back(std::make_pair(endPositionH(back(seedChain)), endPositionV(back(seedChain))));
		endPos.push_back(std::make_pair(length(seqOne), length(seqTwo)));
	}
	//finding and ordered insertion of the new seeds into the global seedChain
	if(!insert_seeds_in_initial_seedChain(seedChain, startPos, endPos, seqOne, seqTwo, options.q)){
		std::cerr << "error: couldn't insert new found seeds!";
		return 1;
	}
	//erase(seedChain, 8);
	//erase(seedChain, 16);
	for(unsigned i = 0; i < length(seedChain); ++i){
		std::cout << "seedChain[" << i << "] = " << seedChain[i] << std::endl;
	}

	//global alignment with seedChain
	Align<DnaString, ArrayGaps> alignment;
	resize(seqan::rows(alignment), 2);
	Score<int, Simple> scoringFunction(2, -1, -2);
	seqan::assignSource(seqan::row(alignment, 0), seqOne);
	seqan::assignSource(seqan::row(alignment, 1), seqTwo);

	int result = bandedChainAlignment(alignment, seedChain, scoringFunction, 2);

	//output
	output_alignment(options.outputFile, idOne, idTwo, alignment, result);
	return 0;
}
