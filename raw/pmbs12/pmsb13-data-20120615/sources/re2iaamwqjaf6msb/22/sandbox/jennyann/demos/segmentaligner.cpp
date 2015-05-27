#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

int main() {

	typedef String<Dna5> TSequence;

	// import query sequences
	StringSet<TSequence> queries;
	StringSet<CharString> queryIDs;
    if (!_importSequences(options.queryFile, "query", queries, queryIDs)) return 1;

	// import database sequence
    StringSet<TSequence > databases;
    StringSet<CharString> databaseIDs;
    if (!_importSequences(options.databaseFile, "database", databases, databaseIDs)) return 1;


	typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
    TFinder swiftFinder(database);

	stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop, 
            options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
            databaseID, databaseStrand, matches, AllLocal());


}