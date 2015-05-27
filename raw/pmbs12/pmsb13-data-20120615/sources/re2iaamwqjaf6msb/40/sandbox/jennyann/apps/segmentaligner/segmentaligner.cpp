#include <iostream>
#include <seqan/sequence.h>
#include "D:/SeqAn/core/apps/stellar/stellar.h"
#include "D:/SeqAn/core/apps/stellar/stellar.cpp"
#include "../../../core/apps/stellar/stellar.h"

using namespace seqan;

int main() {

	typedef String<Dna> TSeuquence;
	TSequence file1 = "AGTTGTTAGTCTACGTGGACCGACAAAGACAGATTCTTTGAGGGAGCTAAGCTCAACGTAGTTCTAACAG";
	TSequence file2 = "AGTTATTAGTCTACGTGGACCGACAAAGAAAGATTCTTTGAGGGAGCTAAGCTCACAATAGTTCTAACAG";

	//typedef String<Dna5> TSequence;

	// import query sequences
	//StringSet<TSequence> queries;
	//StringSet<CharString> queryIDs;
    //if (!_importSequences(file1, "query", queries, queryIDs)) return 1;

	// import database sequence
    //StringSet<TSequence > databases;
    //StringSet<CharString> databaseIDs;
    //if (!_importSequences(file2, "database", databases, databaseIDs)) return 1;

	double epsilon;
	double minLength;

	unsigned epsilon = 0.01;
	unsigned minLength = 20;

	typedef String<StellarMatch<TSequence, TId> matches;

	//typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
    //TQGramIndex qgramIndex(queries);
    //Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);

	//typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
    //TFinder swiftFinder(databases, 1, 1000);

	stellar(file1, file2, epsilon, minLength, xDrop, matches, AllLocal());


}