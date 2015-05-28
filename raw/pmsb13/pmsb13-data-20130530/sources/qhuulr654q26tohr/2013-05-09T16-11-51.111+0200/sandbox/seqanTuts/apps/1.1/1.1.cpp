#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main ()
{
	StringSet< String<Dna> > mySet;
	resize(mySet, 2);
	mySet[0] = "AAAAAAAAAAAAAA";
	mySet[1] = "AAAAATAAA";
	
	typedef Index< StringSet<String<Dna> >, IndexQGram<UngappedShape<2> > > TIndex;
	typedef Infix<Fibre<TIndex, QGramCounts>::Type const>::Type TCounts;
	
	TIndex myIndex(mySet);
	
	String<double> distMat;
	getKmerSimilarityMatrix(myIndex,distMat);
	
	for( unsigned i=0; i < length(distMat); ++i)
		std::cout << distMat[i] << " ";
	std::cout << std::endl;
	
	return 0;
}
