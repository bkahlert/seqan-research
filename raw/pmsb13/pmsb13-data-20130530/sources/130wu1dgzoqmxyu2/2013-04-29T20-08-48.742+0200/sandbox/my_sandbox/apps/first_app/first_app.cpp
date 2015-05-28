#include <iostream>
#include <seqan/index.h>
#include <vector>

using namespace seqan;
using namespace std;


int main(){
	DnaString window1 = "AAAACACGC";
	DnaString window2 = "GGTCGACCGT";

	StringSet<DnaString > stringSet;
	appendValue(stringSet, window1);
	appendValue(stringSet, window2);

	typedef Index< StringSet<DnaString>, IndexQGram<UngappedShape<2> > > qGramIndex;
	typedef typename Fibre<qGramIndex, QGramCountsDir>::Type TCountsDir;
	typedef typename Iterator<TCountsDir, Standard>::Type TIterCountsDir;

	qGramIndex index(stringSet);
	Finder<qGramIndex> myFinder(index);

	while (find(myFinder, "CG")){
		std::cout << position(myFinder) << "  ";
	}
	//TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
	//TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
	std::cout << std::endl << countSequences(index) << std::endl;
	String<Dna> kmer = "CG";
	unsigned hashValue = hash(indexShape(index), begin(kmer));
	for(unsigned i = dirAt(hashValue, index); i < dirAt(hashValue +1, index); ++i)
		std::cout << getSeqNo(saAt(i,index)) << ":" << getSeqOffset(saAt(i,index)) << std::endl;
	/*for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir){
		std::cout << *itCountsDir << std::endl;
	}*/
	return 0;
}
