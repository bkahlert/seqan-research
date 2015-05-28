#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main ()
{
	StringSet< CharString > mySet;
	resize(mySet, 2);
	mySet[0] = "QuerySeq";
	mySet[1] = "ClusterSeq";
	
	typedef Index< StringSet<CharString >, IndexQGram<UngappedShape<8> > > TIndex;
	
	TIndex myIndex(mySet);
	String<double> distMat;
	getKmerSimilarityMatrix(myIndex,distMat);
	
	return 0;
}


//for( unsigned i=0; i < length(distMat); ++i)
//	std::cout << distMat[i] << " ";
//std::cout << std::endl;
