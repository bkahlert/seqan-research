#include "lagan_functions.h"
#include <seqan/arg_parse.h>

using namespace seqan;

struct LaganOptions
{
	//LAGAN parameters
	seqan::CharString seqOne;
	seqan::CharString seqTwo;
	unsigned q;
	seqan::CharString outputFile;

	//stellar parameters
	unsigned repeatLength;
	unsigned repeatPeriod;
	double epsilon;
	int minLength;
	double xDrop;
	unsigned stellarKmer;
	unsigned sortThresh;
	unsigned disableThresh;
	unsigned numMatches;


	LaganOptions() {
		outputFile = "lagan.txt";
		q = 12;

		repeatLength = 1000;
		repeatPeriod = 1;
		epsilon = 0.05;
		minLength = 100;
		xDrop = 5;

		sortThresh = 500;
		disableThresh = 500;
		numMatches = 50;
	}

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
	setValidValues(parser, 0, "fa fasta");
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "FASTA FILE 2"));
	setValidValues(parser, 1, "fa fasta");

	// Define Options
	addSection(parser, " LAGAN options");

	addOption(parser, seqan::ArgParseOption(
			"q", "qGram", "qGram length to use for the index.",
			seqan::ArgParseArgument::INTEGER, "UNSIGNED"));
	setMinValue(parser, "q", "3");
	setMaxValue(parser, "q", "14");

	addOption(parser, seqan::ArgParseOption(
			"o", "ouputFilename", "name of the output file.",
			seqan::ArgParseArgument::OUTPUTFILE));

	addSection(parser, " Stellar: Main Options");
	addOption(parser, ArgParseOption("e", "epsilon", "Maximal error rate (max 0.25).", ArgParseArgument::DOUBLE));
	setMinValue(parser, "e", "0.0000001");
	setMaxValue(parser, "e", "0.25");
	addOption(parser, ArgParseOption("l", "minLength", "Minimal length of epsilon-matches.", ArgParseArgument::INTEGER));
	setMinValue(parser, "l", "0");

	addSection(parser, " Stellar: Filtering Options");

	addOption(parser, ArgParseOption("k", "stellarKmer", "Length of the q-grams in stellar (max 32).", ArgParseArgument::INTEGER));
	setMinValue(parser, "k", "1");
	setMaxValue(parser, "k", "32");
	addOption(parser, ArgParseOption("rp", "repeatPeriod",
			"Maximal period of low complexity repeats to be filtered.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("rl", "repeatLength",
			"Minimal length of low complexity repeats to be filtered.", ArgParseArgument::INTEGER));

	addSection(parser, "Stellar: Verification Options");

	addOption(parser, ArgParseOption("x", "xDrop", "Maximal x-drop for extension.", ArgParseArgument::DOUBLE));
	addOption(parser, ArgParseOption("dt", "disableThresh",
			"Maximal number of verified matches before disabling verification for one query "
			"sequence (default infinity).", ArgParseArgument::INTEGER));
	setMinValue(parser, "dt", "0");
	addOption(parser, ArgParseOption("n", "numMatches",
			"Maximal number of kept matches per query and database. If STELLAR finds more matches, "
			"only the longest ones are kept.", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("s", "sortThresh",
			"Number of matches triggering removal of duplicates. Choose a smaller value for saving "
			"space.", ArgParseArgument::INTEGER));

	// Parse command line.
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	//Only extract  options if the program will continue after parseCommandLine()
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	// Extract option values.
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
	if (!isSet(parser, "k"))
	{
		options.stellarKmer =  (unsigned) (double) 1 / options.epsilon - 1;
		if (options.stellarKmer > 32)
		{
			options.stellarKmer = 32;
		}
	}

	if (isSet(parser, "k"))
	{
		if (options.stellarKmer >= (unsigned) (double) 1 / options.epsilon)
		{
			std::cerr << "k need to be smaller than 1 / epsilon" << std::endl;
			return seqan::ArgumentParser::PARSE_ERROR;
		}
	}
	return seqan::ArgumentParser::PARSE_OK;
}

bool call_stellar(seqan::StringSet<QueryMatches<StellarMatch<seqan::String<seqan::Dna>, seqan::CharString> > > &matches,
							LaganOptions &options,
							seqan::String<seqan::Dna> &seqOne,
							seqan::CharString &idOne,
							seqan::StringSet<seqan::String<seqan::Dna> > &setTwo)
{
	typedef Index<StringSet<String<seqan::Dna>, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
	TQGramIndex qgramIndex(setTwo);
	resize(indexShape(qgramIndex), options.stellarKmer);
	Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);

	indexRequire(qgramIndex, QGramSADir());

	typedef Finder<seqan::String<seqan::Dna>, Swift<SwiftLocal> > TFinder;
	TFinder swiftFinder(seqOne, options.repeatLength, options.repeatPeriod);

	stellar(swiftFinder,
			swiftPattern,
			options.epsilon,
			options.minLength,
			options.xDrop,
			options.disableThresh,
			options.sortThresh,
			options.numMatches,
			true, 					//verbose
			idOne,					//ID of first sequence
			true, 					//forward search only
			matches, 				//container for matches
			AllLocal());

	return true;
}

int main(int argc, char const ** argv){

	LaganOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		return res == seqan::ArgumentParser::PARSE_ERROR;
	}
	//debugging:
	std::cout << "epsilon: " << options.epsilon << std::endl;
	std::cout << "stellarKmer: " << options.stellarKmer << std::endl;
	std::cout << "dThresh: " << options.disableThresh << std::endl;
	std::cout << "minLength: " << options.minLength << std::endl;
	std::cout << "numMatches: " << options.numMatches << std::endl;
	std::cout << "q: " << options.q << std::endl;
	std::cout << "repeatLength: " << options.repeatLength << std::endl;
	std::cout << "repeatPeriod: " << options.repeatPeriod << std::endl;
	std::cout << "sortThresh: " << options.sortThresh << std::endl;
	std::cout << "seqOne: " << options.seqOne << std::endl;
	std::cout << "seqTwo: " << options.seqTwo << std::endl;
	std::cout << "floor(1/eps): " << (unsigned) (double) 1 /options.epsilon << std::endl << std::endl;

	seqan::CharString idOne;
	seqan::String<seqan::Dna> seqOne;
	seqan::CharString idTwo;
	seqan::String<seqan::Dna> seqTwo;

	typedef Seed<Simple> SSeed;
	typedef SeedSet<Simple> SSeedSet;

	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq1;
	std::vector< std::pair<unsigned, unsigned> > seedsPosSeq2;

	//read-in of FASTA sequences
	if (!read_FASTA(options.seqOne, idOne, seqOne))
	{
		std::cerr << "error: unable to read first sequence";
		return 1;
	}
	if (!read_FASTA(options.seqTwo, idTwo, seqTwo))
	{
		std::cerr << "error: unable to read second sequence";
		return 1;
	}

	//call stellar to obtain initial seeds
	std::cout << "STELLAR: searching for seeds..." << std::endl;

	seqan::StringSet<seqan::String<seqan::Dna> > setTwo;
	seqan::appendValue(setTwo, seqTwo);

	StringSet<QueryMatches<StellarMatch<seqan::String<seqan::Dna>, seqan::CharString> > > matches;
	resize(matches, 1);

	if(!call_stellar(matches, options, seqOne, idOne, setTwo))
	{
		std::cerr << "error: Couldn't call stellar." << std::endl;
	}

	/*[deprecated]
	 * if (!write_seed_positions(seedsPosSeq1, seedsPosSeq2, argv[3])){
		std::cerr << "error: STELLAR output file not found";
		return 1;
	}
	 */
	std::cout << std::endl << "STELLAR: done." << std::endl;

	//extracting positions from stellar matches (seeds) and compute best chain
	SSeedSet seedSet;

	if (!get_stellarMatches(matches, seedSet))
	{
		return 1;
	}

	std::cout << std::endl << "Creating initial seedChain...." << std::endl;

	seqan::String<SSeed> seedChain;
	seqan::chainSeedsGlobally(seedChain, seedSet, seqan::SparseChaining());


	//old: for reading stellarmatches from file
	//writing seeds and adding them to a SeedSet
	/*for (unsigned i = 0; i < seedsPosSeq1.size(); ++i){
		SSeed seed(seedsPosSeq1[i].first,
				seedsPosSeq2[i].first,
				seedsPosSeq1[i].second,
				seedsPosSeq2[i].second);

		addSeed(seedSet, seed, Single());
	}
	//std::cout << "Trennlinie" << std::endl;
	typedef Iterator<SSeedSet >::Type SetIterator;

	for (SetIterator it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it){
		std::cout << *it;
		std::cout << std::endl;
	}
	typedef Iterator<String<SSeed> > StringIterator;

	seqan::String<SSeed> seedChain;

	chainSeedsGlobally(seedChain, seedSet, seqan::SparseChaining());
	 */

	//debugging:
	for(unsigned i = 0; i < length(seedChain); ++i)
	{
		std::cout << "Initial: seedChain[" << i << "] = " << seedChain[i] << std::endl;
	}

	std::vector<std::pair<unsigned, unsigned> > startPos;
	std::vector<std::pair<unsigned, unsigned> > endPos;

	//creating a container to define the gaps between the initial stellar matches
	if (!create_positionsContainers(startPos, endPos,
			seedChain,
			options.q,
			length(seqOne), length(seqTwo)))
	{
		return 1;
	}

	std::cout << "Filling gaps between STELLAR seeds..." << std::endl;
	//finding and ordered insertion of the new seeds into the global seedChain
	if(!insert_seeds_in_seedChain(seedChain, startPos, endPos, seqOne, seqTwo, options.q))
	{
		std::cerr << "error: couldn't insert new found seeds!";
		return 1;
	}

	for(unsigned i = 0; i < length(seedChain); ++i)
	{
		std::cout << "seedChain[" << i << "] = " << seedChain[i] << std::endl;
	}
	std::cout << "Creating alignment..." << std::endl;

	for(unsigned i = 0; i < length(seedChain) -1; ++i)
	{
		if(beginPositionH(seedChain[i]) >= beginPositionH(seedChain[i+1]) ||
			beginPositionV(seedChain[i]) >= beginPositionV(seedChain[i+1]) ||
			endPositionH(seedChain[i]) >= endPositionH(seedChain[i+1]) ||
			endPositionV(seedChain[i]) >= endPositionV(seedChain[i+1]))
		{
			std::cout << "ALARM: Irregularitaet bei seedChain[" << i << "]" << std::endl;
		}
	}

	//global alignment with seedChain
	Align<String<seqan::Dna>, ArrayGaps> alignment;
	std::cout << "bla" << std::endl;
	resize(seqan::rows(alignment), 2);
	std::cout << "bla2" << std::endl;
	Score<int, Simple> scoringFunction(5, -1, -1);
	std::cout << "bla3" << std::endl;
	seqan::assignSource(seqan::row(alignment, 0), seqOne);
	std::cout << "bla4" << std::endl;
	seqan::assignSource(seqan::row(alignment, 1), seqTwo);
	std::cout << "bla5" << std::endl;
	int result = bandedChainAlignment(alignment, seedChain, scoringFunction, 2);

	//output
	std::cout << "Creating output..." << std::endl;
	if (!output_alignment(options.outputFile, idOne, idTwo, alignment, result))
	{
		std::cerr << "Unable to open output file";
		return 1;
	}
	return 0;
}
