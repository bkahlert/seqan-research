#include <iostream>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <vector>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

using namespace seqan;
using namespace std;

void assignmentIndex2(){
	String<char> text = "TTATTAAGCGTATAGCCCTATAAATATAA";
	StringSet<String<char> > set;

	String<char> pattern1 = "TA";
	String<char> pattern2 = "GCG";
	String<char> pattern3 = "GCC";

	appendValue(set, pattern1);
	appendValue(set, pattern2);
	appendValue(set, pattern3);

	Index<String<char>, IndexEsa<> > esaIndex(text);
	Finder<Index<String<char>, IndexEsa<> > > esaFinder(esaIndex);

	typedef Iterator<StringSet<String<char> > >::Type TStringSetIterator;
	for (TStringSetIterator it = begin(set); it != end(set); ++it){
		clear(esaFinder);
		std::cout << value(it) << ": ";
		while(find(esaFinder, getValue(it))){
			std::cout<<position(esaFinder) << " ";
		}
		std::cout << "\n";
	}
}

void assignmentIndexIterators2(){
	String<char> text = "tobeornottobe";
	typedef Index<String<char> > TIndex;
	TIndex index(text);

	Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);
	bool stop = false;
	while(!stop){
		cout << representative(it) << endl;
		if (!goDown(it))
			while(!goRight(it)){
				goUp(it);
				if (isRoot(it)){
					stop = true;
					break;
				}
			}
	}
}

void assignment1qGram(){
	typedef Index<DnaString, IndexQGram<OneGappedShape> > TIndex;
	TIndex index("CATGATTACATA");
	stringToShape(indexShape(index), "1101");
	hash(indexShape(index), "AT-A");
	for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i){
		std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

	}
}

void assignment1Graphs(){
	typedef unsigned int TCargo;
	typedef unsigned int TSpec;
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TGraph G;

	vector<int> edges;
	edges.push_back(1); edges.push_back(0); edges.push_back(0); edges.push_back(4);
	edges.push_back(2); edges.push_back(1); edges.push_back(4); edges.push_back(1);
	edges.push_back(5); edges.push_back(1); edges.push_back(6); edges.push_back(3);
	edges.push_back(2); edges.push_back(2); edges.push_back(3); edges.push_back(7);
	edges.push_back(3); edges.push_back(7); edges.push_back(5); edges.push_back(4);
	edges.push_back(6); edges.push_back(5); edges.push_back(5); edges.push_back(6);
	edges.push_back(7); edges.push_back(6); edges.push_back(7); edges.push_back(7);

	addEdges(G, edges, 14);

	FILE* strmWrite = fopen("graph.dot", "w");
	write(strmWrite, G, DotDrawing());
	fclose(strmWrite);

}

int main(){
	typedef seqan::Seed<seqan::Simple>    TSeed;
	seqan::Dna5String sequenceH = "CGAATCCATCCCACACA";
	    seqan::Dna5String sequenceV = "GGCGATNNNCATGGCACA";
	    seqan::Score<int, seqan::Simple> scoringSchemeAnchor(0, -1, -1);
	    seqan::Score<int, seqan::Simple> scoringSchemeGap(2, -1, -1, -2);
	    seqan::String<TSeed> seedChain;
	    seqan::appendValue(seedChain, TSeed(0, 2, 5, 6));
	    seqan::appendValue(seedChain, TSeed(6, 9, 9, 12));
	    seqan::appendValue(seedChain, TSeed(11, 14, 17, 16));
	    seqan::Align<seqan::Dna5String, seqan::ArrayGaps> alignment;
	    seqan::resize(seqan::rows(alignment), 2);
	    seqan::assignSource(seqan::row(alignment, 0), sequenceH);
	    seqan::assignSource(seqan::row(alignment, 1), sequenceV);
	    seqan::AlignConfig<true, false, false, true> alignConfig;
	    int result = seqan::bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2);
	    std::cout << "Score: " << result << std::endl;
	    std::cout << alignment << std::endl;
	    return 0;
}
