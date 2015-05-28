#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main ()
{
	StringSet< DnaString > mySet;
	resize(mySet, 2);
	mySet[0] = "AAAACCCCTTTTTT";
	mySet[1] = "GGCCCCC";
	
	typedef Index< StringSet<DnaString >, IndexQGram<UngappedShape<4> > > TIndex;
	
	TIndex myIndex(mySet);
	String<double> distMat;
	getKmerSimilarityMatrix(myIndex,distMat);
	for( unsigned i=0; i < length(distMat); ++i)
		std::cout << distMat[i] << " ";
	std::cout << std::endl;
	return 0;
}



