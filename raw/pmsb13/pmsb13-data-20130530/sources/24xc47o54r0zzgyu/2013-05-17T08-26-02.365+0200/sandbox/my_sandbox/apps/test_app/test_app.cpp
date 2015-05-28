// Program part for parsing the input
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(){
	StringSet<DnaString> stringSet;
	DnaString str0 = "TATA";
	DnaString str1 = "CGCG";
	appendValue(stringSet, str0);
	appendValue(stringSet, str1);

	unsigned id0 = positionToId(stringSet, 0);
	unsigned id1 = positionToId(stringSet, 1);

	std::cout << id0 << " " << id1 << "\n";
	
	removeValueById(stringSet, id0);

	std::cout << positionToId(stringSet, 0) << "\n";
	removeValueById(stringSet, id1);
	return 0;
}
