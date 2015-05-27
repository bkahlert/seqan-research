#include <fstream>
#include <iostream>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/graph_msa.h>
#include <seqan/index.h>

#include "../../../../core/apps/stellar/stellar.h"
#include "../../../../core/apps/stellar/stellar_types.h"
#include "../../../../core/apps/stellar/stellar_output.h"

using namespace seqan;

////////////////reading file
template <typename TSequence, typename TId>
inline bool seq_parser(const char* & fileName,
	                   TSequence & seq,
	                   TId & id) {
// FRAGMENT(read-sequences)
    // Open file and create RecordReader.
    std::fstream fasta(toCString(fileName), std::ios_base::in | std::ios_base::binary);
    if (!fasta.good()) {
       std::cerr << "Failed to open " << fileName << " file." << std::endl;
	   return false;
	}

    RecordReader<std::fstream, SinglePass< > > reader(fasta);

    while (!atEnd(reader))
    {
        if (readRecord(id, seq, reader, Fasta()) != 0)
        {
            std::cerr << "ERROR reading FASTA." << std::endl;
            return false;
        }
		return true;
    }
	return false;
}

template <typename StMatch, typename TAlign>
void mergeAlignments(String<StringSet<QueryMatches<StMatch> > > & DBmatches, String<TAlign> & MergedAlignment) {
	
	typedef typename Iterator<String<StMatch> >::Type TIterator;
	
	 for (int x = 0; x < length(DBmatches); x++) {
		 for (int y = 0; y < length(DBmatches[x]); y++) {
			QueryMatches<StMatch> &queryMatches = value(DBmatches[x], y);
			TIterator it = begin(queryMatches.matches);
			TIterator itEnd = end(queryMatches.matches);
			//std::cout << length(DBmatches[x]) << std::endl;
			int counter = 0;
			while (it < itEnd) {
			TAlign align;
			appendValue(align.data_rows, (*it).row1);
			appendValue(align.data_rows, (*it).row2);
			std::cout << "row1" << (*it).row1 << std::endl;
			std::cout << "row2" << (*it).row2 << std::endl;
			appendValue(MergedAlignment,align);
			//std::cout << "length" << length(dbAlignment) << std::endl;
			it++;
			//std::cout << align << std::endl;
			//std::cout << ++counter << std::endl;
			}
			//for(int o = 0; o < length(MergedAlignment); ++o)
			//std::cout << o << MergedAlignment[o] << std::endl;
			//std::cout <<  "-------------\n";
		 }
	 }
}

///////////////////////////calls stellar
template<typename TSequence, typename Tmatches>
void stellar_find_matches(TSequence & database, StringSet<TSequence> & queries, Tmatches & matches, String<Tmatches> &DBmatches) {
	double epsilon = 0.05;
	int minLength = 20;
	unsigned xDrop = 2;
	unsigned qGram = 3;

	clear(matches);
    resize(matches, length(queries));

	//pattern
	typedef Index<StringSet<TSequence>, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
    TQGramIndex qgramIndex(queries);
	resize(indexShape(qgramIndex), qGram);
    Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);

	typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
    TFinder swiftFinder(database);

	stellar(swiftFinder, swiftPattern, epsilon, minLength, xDrop, matches, AllLocal() );
		appendValue(DBmatches, matches);
}

//////////////////////////////merges Alingments
template<typename TSequence, typename Tmatches, typename TAlign>
void concatMatches(StringSet<TSequence> & seq, String<Tmatches> & DBmatches, String<TAlign> & complete)
{
	TSequence database;
	TSequence query;
	StringSet<TSequence> queries;
	Tmatches matches;

    for(unsigned int i = 0; i < length(seq)-1; i++)
    {
			database = seq[i];
        for(unsigned int j = i+1; j < length(seq); j++)
        {
			 query = seq[j];
			appendValue(queries, query);
		 }
		std::cout << "query" << length(queries) << std::endl;
		stellar_find_matches(database, queries,matches,DBmatches);
		String<TAlign> MergedAlignment;
		mergeAlignments(DBmatches, MergedAlignment);
		append(complete, MergedAlignment);

		clear(DBmatches);
		clear(matches);
		clear(database);
		clear(queries);
    }
}

//////////////////////////////builds outGraph
template <typename TStringSet, typename TAlignmentGraph>
void buildoutGraph (TAlignmentGraph & ali_graph, TStringSet & seq_set, Graph<Tree<double> > &guideTree, TAlignmentGraph & outGraph) {
	
	tripletLibraryExtension(ali_graph);
	
	typedef String<double> TDistanceMatrix;
	TDistanceMatrix distanceMatrix;
	getDistanceMatrix(ali_graph, distanceMatrix, KmerDistance());

	// Guide Tree
	upgmaTree(distanceMatrix, guideTree);
	clear(distanceMatrix);

	// Perform a progressive alignment
	progressiveAlignment(ali_graph, guideTree, outGraph);
	clear(ali_graph);
	clear(guideTree);
	std::cout << outGraph << std::endl;
}

//////////////////////////////////test function for matchRefignment//////////////////
template <typename TSequence, typename TAlign>
void test(TSequence & seqs,String<TAlign> &matches) {
	
	for(unsigned int i = 0; i < length(seqs)-1; ++i)
	{
		for(unsigned int j = i+1; j < length(seqs); ++j)
		{
			TAlign ali;
			resize(rows(ali), 2);
			setSource(row(ali, 0), seqs[i]);
			setSource(row(ali, 1), seqs[j]);
			LocalAlignmentFinder<> finder(ali);
			Score<int> scoring(3, -2, -1, -5);
			localAlignment(ali, finder, scoring, 0, WatermanEggert());			
			appendValue(matches, ali);
			clear(rows(ali));
		}	
	}
}

/////////////////////main function ///////////////////////
int main(int argc, char const * argv[])
{
	if (argc < 2)
    {
        std::cerr << "ERROR: Invalid argument count." << std::endl
                  << "USAGE: " << argv[0] << " FILE" << std::endl;
		return 1;
    }

	int numSequences = argc-1;

	typedef String<Dna5> TSequence;
	typedef StringSet<TSequence> TStringSet;

	TSequence sequence;
	CharString sequenceID;
	TStringSet seq_set;

	for(int i = 1; i < numSequences+1; i++)
	{
		if (!seq_parser(argv[i], sequence, sequenceID)) return 1;
		appendValue(seq_set, sequence);
	}

	typedef String<char> TId;
	typedef StellarMatch<TSequence, TId> StMatch;

	// container for eps-matches
	typedef String<char> TId;
	typedef StringSet<QueryMatches<StMatch> > Tmatches;
	String<Tmatches> DBmatches;

	typedef Align<TSequence, ArrayGaps> TAlign;
	String<TAlign> complete;

	concatMatches(seq_set,DBmatches,complete);
	//std::cout << length(DBmatches) << std::endl;

	for(int b = 2; b < length(complete); b++)
		std::cout << "here \n" << complete[b] << std::endl;

	clear(DBmatches);
	//std::cout << length(alignmentmerged) << std::endl;
		/*for(int b = 2; b < length(alignmentmerged); b++)
	{
		std::cout << "here \n" << alignmentmerged[b] << std::endl;
	}*/

	//test(seq_set,alignmentmerged);

	typedef StringSet<TSequence,Dependent<> > TDepSequence;
	typedef Graph<Alignment<TDepSequence> > TAlignmentGraph;
    TAlignmentGraph ali_graph(seq_set);

	/*matchRefinement(alignmentmerged, seq_set, ali_graph);

	clear(alignmentmerged);
	
	typedef Graph<Tree<double> > TGuideTree;
	TGuideTree guideTree;
	TAlignmentGraph outGraph(seq_set);

	buildoutGraph(ali_graph,seq_set,guideTree, outGraph);*/
	return 0;
}