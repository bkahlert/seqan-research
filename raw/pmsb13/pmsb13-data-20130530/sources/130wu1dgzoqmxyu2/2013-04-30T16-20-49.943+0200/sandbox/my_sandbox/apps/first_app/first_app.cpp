#include <iostream>
#include <seqan/index.h>
#include <vector>

using namespace seqan;
using namespace std;


int main(){
	DnaString window1 = "AAAACACGC";
	DnaString window2 = "GGTCGACCGT";

	typedef Index< DnaString, IndexQGram<UngappedShape<2> > > qGramIndex;

	qGramIndex index(window1);
	Finder<qGramIndex> myFinder(index);
	std::cout << length(window2) << std::endl;
	for (unsigned i = 0; i < length(window2) - 2; ++i){
		DnaString qGram = infix(window2, i, i+2);
		std::cout << qGram << std::endl;
		while (find(myFinder, qGram)){
			std::cout << position(myFinder) << std::endl;
		}
		clear(myFinder);
	}
	//TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
	//TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
	//std::cout << std::endl << countSequences(index) << std::endl;
	String<Dna> kmer = "CG";
	unsigned hashValue = hash(indexShape(index), begin(kmer));
	for(unsigned i = dirAt(hashValue, index); i < dirAt(hashValue +1, index); ++i){
		std::cout << getSeqNo(saAt(i,index)) << ":" << getSeqOffset(saAt(i,index)) << std::endl;
	}
	/*for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir){
		std::cout << *itCountsDir << std::endl;
	}*/
	return 0;
}
