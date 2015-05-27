#include <iostream>
#include <seqan/sequence.h>
#include "D:/SeqAn/core/apps/stellar/stellar.h"
#include "D:/SeqAn/core/apps/stellar/stellar.cpp"
#include "../../../core/apps/stellar/stellar.h"

using namespace seqan;

int main() {

	typedef String<char> Tfiles;
	Tfiles file1 = "NC_001474.fasta";
	Tfiles file2 = "NC_001477.fasta";

	typedef String<Dna5> TSequence;

	// import query sequences
	StringSet<TSequence> queries;
	StringSet<CharString> queryIDs;
    if (!_importSequences(file1, "query", queries, queryIDs)) return 1;

	// import database sequence
    StringSet<TSequence > databases;
    StringSet<CharString> databaseIDs;
    if (!_importSequences(file2, "database", databases, databaseIDs)) return 1;

	double epsilon;
	double minLength;

	unsigned epsilon = 0.1;
	unsigned minLength = 50;
	unsigned maxValue = (unsigned)-1;
    unsigned compactThresh = 1000;
	TId id = "db";

	typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
    TQGramIndex qgramIndex(queries);
    Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);

	

	typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
    TFinder swiftFinder(databases, 1, 1000);

	//stellar(file1, file2, epsilon, minLength, AllLocal());

	stellar(swiftFinder, swiftPattern, epsilon, minLength, xDrop,
            maxValue, compactThresh, maxValue, 0, id, true,
            matches, AllLocal());


}