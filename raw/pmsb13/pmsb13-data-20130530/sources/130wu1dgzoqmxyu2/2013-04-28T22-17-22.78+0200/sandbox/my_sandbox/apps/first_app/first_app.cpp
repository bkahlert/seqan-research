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

template <typename TStringSet, typename TIndexSpec>
void qgramCounting(TStringSet &set, TIndexSpec)
{
    typedef Index<TStringSet, TIndexSpec> TIndex;
    typedef typename Fibre<TIndex, QGramCounts>::Type TCounts;
    typedef typename Fibre<TIndex, QGramCountsDir>::Type TCountsDir;
    typedef typename Value<TCountsDir>::Type TDirValue;
    typedef typename Iterator<TCounts, Standard>::Type TIterCounts;
    typedef typename Iterator<TCountsDir, Standard>::Type TIterCountsDir;

    TIndex index(set);
    indexRequire(index, QGramCounts());

    // initialize distance matrix
    int seqNum = countSequences(index);
    Matrix<int, 2> distMat;
    setLength(distMat, 0, seqNum);
    setLength(distMat, 1, seqNum);
    resize(distMat, 0);

    std::cout << std::endl << "Length of the CountsDir fibre: " << length(indexCountsDir(index)) << std::endl;
    TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
    TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
    TIterCounts itCountsBegin = begin(indexCounts(index), Standard());
    // for each bucket count common q-grams for each sequence pair
      TDirValue bucketBegin = *itCountsDir;
      for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir)
      {
          TDirValue bucketEnd = *itCountsDir;

          // q-gram must occur in at least 2 different sequences
          if (bucketBegin != bucketEnd)
          {
              TIterCounts itA = itCountsBegin + bucketBegin;
              TIterCounts itEnd = itCountsBegin + bucketEnd;
              for(; itA != itEnd; ++itA)
                  for(TIterCounts itB = itA; itB != itEnd; ++itB)
                      distMat((*itA).i1, (*itB).i1) += _min((*itA).i2, (*itB).i2);
          }
          bucketBegin = bucketEnd;
      }

      std::cout << std::endl << "Common 5-mers for Seq_i, Seq_j" << std::endl;
      std::cout << distMat;
}


int main(){
	String<char> window1 = "AAAACACGC";
	String<char> window2 = "GGTCGACCGT";

	StringSet<String<char> > stringSet;
	appendValue(stringSet, window1);
	appendValue(stringSet, window2);

	Index< StringSet<String<char> >, IndexQGram<UngappedShape<2> > > TQGramIndex;
	TQGramIndex myIndex(stringSet);
	DnaString kmer = "CG";
	unsigned hashValue = hash(indexShape(index), begin(kmer));
	for(unsigned i = dirAt(hashValue, index); i < dirAt(hashValue +1, index); ++i){
	  std::cout << getSeqNo(saAt(i,index)) << ":" << getSeqOffset(saAt(i,index)) << ::std::endl;
	}
	  /*StringSet< String<char> > mySet;
	resize(mySet, 4);
	mySet[0] = "tobeornottobe";
	mySet[1] = "thebeeonthecomb";
	mySet[2] = "hellobebe";
	mySet[3] = "beingjohnmalkovich";

	typedef Index< StringSet<String<char> >, IndexQGram<UngappedShape<2> > > TIndex;
	typedef Infix<Fibre<TIndex, QGramCounts>::Type const>::Type TCounts;

	TIndex myIndex(mySet);
	std::cout << "Number of sequences: " << countSequences(myIndex) << std::endl;
	hash(indexShape(myIndex), "be");
	TCounts cnts = countOccurrencesMultiple(myIndex, indexShape(myIndex));
	for (unsigned i = 0; i < length(cnts); ++i){
		std::cout << cnts[i].i2 << " occurrences in sequence " << cnts[i].i1  << std::endl;
	}

	String<double> distMat;
	getKmerSimilarityMatrix(myIndex,distMat);

	for( unsigned i=0; i < length(distMat); ++i){
		std::cout << distMat[i] << " ";
	}
	std::cout << std::endl;

	return 0;
	*/
	/*typedef seqan::Seed<seqan::Simple>    TSeed;
	seqan::Dna5String sequenceH = "CGAATCCATCCCACACA";
	    seqan::Dna5String sequenceV = "GGCGATNNNCATGGCACA";
	    seqan::Score<int,  StringSet< String<char> > mySet;
9	    resize(mySet, 4);
10	    mySet[0] = "tobeornottobe";
11	    mySet[1] = "thebeeonthecomb";
12	    mySet[2] = "hellobebe";
13	    mySet[3] = "beingjohnmalkovich";::Simple> scoringSchemeAnchor(0, -1, -1);
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
	    */
}
