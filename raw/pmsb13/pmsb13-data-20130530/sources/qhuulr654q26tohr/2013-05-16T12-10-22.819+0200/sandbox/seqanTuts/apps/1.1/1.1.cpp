#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main ()
{
	StringSet< String<CharString> > mySet;
	resize(mySet, 2);
	mySet[0] = "DDDDtesttestQuerySeqDDDDtesttestQuerySeq";
	mySet[1] = "DDDDtesttestClusterSeqDDDDtesttestClusterSeq";
	
	typedef Index< StringSet<String<CharString> >, IndexQGram<UngappedShape<8> > > TIndex;
	typedef Infix<Fibre<TIndex, QGramCounts>::Type const>::Type TCounts;
	
	TIndex myIndex(mySet);
	
	String<double> distMat;
	getKmerSimilarityMatrix(myIndex,distMat);
	
	for( unsigned i=0; i < length(distMat); ++i)
		std::cout << distMat[i] << " ";
	std::cout << std::endl;
	
	return 0;
}
