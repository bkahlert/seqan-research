#include <iostream>
#include <seqan/index.h>
#include <vector>

using namespace seqan;
using namespace std;


int main(){
	DnaString window1 = "AAAACACGCGGTC";
	DnaString window2 = "GGTCGACCGT";
	unsigned q = 4;
	typedef Index< DnaString, IndexQGram<SimpleShape > > qGramIndex;
	qGramIndex index(window1);
	resize(indexShape(index), q);
	Finder<qGramIndex> myFinder(index);
	std::cout << length(window2) << std::endl;
	for (unsigned i = 0; i < length(window2) - (q - 1); ++i){
		DnaString qGram = infix(window2, i, i+q);
		std::cout << qGram << std::endl;
		while (find(myFinder, qGram)){
			std::cout << position(myFinder) << std::endl;
		}
		clear(myFinder);
	}
	return 0;
}
