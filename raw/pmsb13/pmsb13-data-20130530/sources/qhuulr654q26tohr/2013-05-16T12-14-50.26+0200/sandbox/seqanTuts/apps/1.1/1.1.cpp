#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main ()
{
	StringSet< CharString > mySet;
	resize(mySet, 2);
	mySet[0] = "DDD";
	mySet[1] = "DDD";
	
	typedef Index< StringSet<CharString >, IndexQGram<UngappedShape<8> > > TIndex;
	
	TIndex myIndex(mySet);
	
	String<double> distMat;
	getKmerSimilarityMatrix(myIndex,distMat);
	
	for( unsigned i=0; i < length(distMat); ++i)
		std::cout << distMat[i] << " ";
	std::cout << std::endl;
	
	return 0;
}
