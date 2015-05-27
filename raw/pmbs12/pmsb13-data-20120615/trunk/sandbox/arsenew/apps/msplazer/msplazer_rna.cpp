
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/misc/misc_cmdparser.h>
#include "../../../../core/apps/stellar/stellar.h"
#include "msplazer.h"
#include "msplazer_main.h"

using namespace seqan;

// Program entry point
int main(int argc, char const ** argv)
{

	double start = sysTime();
	// command line parsing
	CommandLineParser parser("msplazer");

	_setParser(parser);
	if (!parse(parser, argc, argv)) {
		if (isSetShort(parser, 'h') || isSetShort(parser, 'v')) return 0; 
		shortHelp(parser, std::cerr);
		return 1;
	}

	/////////////////////////////////////////////////////////////////////
	//Initialisation
	MSplazerOptions msplazerOptions = MSplazerOptions();

	//stellar parameters
	//StellarOptions
	StellarOptions stellarOptions = StellarOptions();

	if (!_parseOptions(parser, stellarOptions, msplazerOptions)) {
		return 1;
	}

	msplazerOptions.databaseFile = stellarOptions.databaseFile;
	msplazerOptions.queryFile = stellarOptions.queryFile;
	
	//_writeFileNames(stellarOptions);
	//_writeParams(msplazerOptions);

    //msplazer wrapper function
	//msplazer(stellarOptions, msplazerOptions);

	//Finder
	typedef String<Dna5Q> TSequence;
	//  typedef FragmentStore<void>::TReadSeqStore TSequence;

	//Database and query ID
	typedef CharString TId;

	// import query sequences using _importSequences from Stellar
	StringSet<TSequence> queries;
	StringSet<TId> queryIDs;
	std::cout << "Loading query sequences... ";
	if (!_importSequences(stellarOptions.queryFile, "query", queries, queryIDs)) return 1;
	
	std::cout << "done" << std::endl;


	// import database sequence using _importSequences from Stellar
	StringSet<TSequence > databases;
	StringSet<TId> databaseIDs;
	std::cout << "Loading reference sequences... ";
	if (!_importSequences(stellarOptions.databaseFile, "database", databases, databaseIDs)) return 1;
	
	std::cout << "done" << std::endl;

	///////////////////////////////////////////////////////////////////////////
	//Compute Stellar Matches and their score
	//Note: Matches will be sorted when calling score function _getMatchDistanceScore
	//Note: Matches of the reverse strand will be modified in the sense that the match positions refer to the forward strand (and not the reverse strand)
	//get distance scores for all matches

	//Stellar Match Space Allocator
	typedef StringSet<QueryMatches<StellarMatch<TSequence, TId> > > TMatches;
	

	TMatches stellarMatches;
	TMatches testStMatches;
	resize(stellarMatches, length(queries));

	//Score space allocator
	typedef StellarMatch<TSequence, TId> TMatch;
	typedef String<int> TScoreAlloc;

	String<TScoreAlloc> distanceScores;
	resize(distanceScores, length(stellarMatches));

	//get Stellar matches
	bool doStellar = true;
	if(msplazerOptions.stellarInputFile != ""){
		doStellar = false;
		std::cout << " Getting STELLAR matches from file, not calling STELLAR " << std::endl;
	}	
	if(doStellar){
		std::cout << "Calling STELLAR... ";

		std::cout << "Stellar options: " << std::endl;

		_writeFileNames(stellarOptions);
		_writeSpecifiedParams(stellarOptions);
		_writeCalculatedParams(stellarOptions);
		_getStellarMatches(queries, databases, databaseIDs, stellarOptions, stellarMatches);

		for(unsigned i = 0; i < length(stellarMatches); ++i){
			for(unsigned j = 0; j < length(stellarMatches[i].matches); ++j){
				std::cout << stellarMatches[i].matches[j] << std::endl;
			}
		
		} 
		std::cout << "done" << std::endl;
	}
	else{
		std::cout << "Importing STELLAR matches from file " << msplazerOptions.stellarInputFile << std::endl;
		double startST = sysTime();
		//_StellarMatchesFromFile(queries, queryIDs, databases, databaseIDs, stellarMatches, msplazerOptions.stellarInputFile);
		//_getStellarMatchesFromFile(queries, queryIDs, databases, databaseIDs, testStMatches, distTestScores, msplazerOptions.stellarInputFile);
		std::cout << "done" << std::endl;
		std::cout << "TIME importing stellar matches " <<	(sysTime() - startST) << "s" << std::endl;
	}

	_getMatchDistanceScore(stellarMatches, distanceScores, doStellar);



	//Graph structure for stellar matches
	typedef Graph<Directed<int> > TGraph;//TRowSize> > TGraph;
	//static_cast<Nothing>(TGraph());
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef Iterator<TGraph, VertexIterator>::Type TEdgeIterator;

	//Breakpoint property map
	typedef SparsePropertyMap < Breakpoint < TSequence, TId >, unsigned > TSparsePropertyMap;
	typedef String<TMatch> TMatchAlloc;

	//Container for graph/chain, scores and start and end vertices
	typedef MSplazerChain< TGraph, TVertexDescriptor, TScoreAlloc, TSparsePropertyMap, TMatchAlloc > TMSplazerChain;
	String< TMSplazerChain> queryChains;

	std::cout << "Constructing graphs... ";
	_chainQueryMatches(stellarMatches, distanceScores, queryChains, queryIDs, queries, msplazerOptions);
	std::cout << "done" << std::endl;

	/////////////////////////////////////////////////////////////////////////
	//Analyze chains

	std::cout << "Analyzing graphs... ";
	_analyzeChains(queryChains);
	std::cout << "done" << std::endl;


	typedef CharString TId;
	
	
	typedef Breakpoint<TSequence,TId> TBreakpoint;
	String<TBreakpoint> globalBreakpoints;
    String<TBreakpoint>spliceSite;
	String<TBreakpoint> globalStellarIndels;
	_findAllBestChains(queryChains, stellarMatches, queries, queryIDs, globalBreakpoints, globalStellarIndels, msplazerOptions);
	
	String<TSequence>genome = databases[0];

	for(unsigned i = 0; i < length(globalBreakpoints); ++i)
	{
	    std::cout << globalBreakpoints[i] << std::endl;

		if(getSVType(globalBreakpoints[i]) == "deletion" )
		{
			std::cout << "deletion" << std::endl;
			if(genome[(globalBreakpoints[i]).startSeqPos+1] == "G" && genome[(globalBreakpoints[i]).startSeqPos+2] == "T")
			{
					if (genome[(globalBreakpoints[i]).endSeqPos-1] == "G" && genome[(globalBreakpoints[i]).endSeqPos-2] == "A")
					{	
						std::cout << globalBreakpoints[i] << std::endl;

					    _insertBreakpoint(spliceSite, globalBreakpoints[i]);
						  
					}
			
		
			}
		}	
	}


	_writeFileNames(stellarOptions);
	_writeParams(msplazerOptions);
	

	_writeGlobalBreakpoints(globalBreakpoints, msplazerOptions, msplazerOptions.support, "breakpoints.gff");
	_writeGlobalBreakpoints(globalStellarIndels, msplazerOptions, msplazerOptions.support, "indels.gff");

	std::cout << "TIME all " <<	(sysTime() - start) << "s" << std::endl;

	
	return 0;
}
