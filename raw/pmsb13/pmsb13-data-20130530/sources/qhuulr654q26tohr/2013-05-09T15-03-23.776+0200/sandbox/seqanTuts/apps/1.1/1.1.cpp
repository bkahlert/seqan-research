#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

int main ()
{
	StringSet< String<Dna> > mySet;
	resize(mySet, 2);
	mySet[0] = "ACTGGCGCGTATAGCCCGGATAACGCCGATAATCGCG";
	mySet[1] = "AGGGCCGCGTATAGCTAGCTGAGATCCCGCTATGTGATCA";
	
	typedef Index< StringSet<String<Dna> >, IndexQGram<UngappedShape<3> > > TIndex;
	typedef Infix<Fibre<TIndex, QGramCounts>::Type const>::Type TCounts;
	
	TIndex myIndex(mySet);
	
	std::cout << "Number of sequences: " << countSequences(myIndex) << std::endl;  
	hash(indexShape(myIndex), "AGT");
	TCounts cnts = countOccurrencesMultiple(myIndex, indexShape(myIndex));
	for (unsigned i = 0; i < length(cnts); ++i)
		std::cout << cnts[i].i2 << " occurrences in sequence " << cnts[i].i1  << std::endl;
	
	
	String<double> distMat;
	getKmerSimilarityMatrix(myIndex,distMat);
	
	for( unsigned i=0; i < length(distMat); ++i)
		std::cout << distMat[i] << " ";
	std::cout << std::endl;
	
	return 0;
}
