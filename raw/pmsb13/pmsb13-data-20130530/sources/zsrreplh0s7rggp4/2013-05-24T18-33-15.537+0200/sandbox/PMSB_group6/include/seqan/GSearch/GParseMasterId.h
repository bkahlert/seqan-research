//Autor:Hannes

#ifndef GINGER_GSEARCH_GPARSEMASTERID_H_
#define GINGER_GSEARCH_GPARSEMASTERID_H_

#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

/**
 * \brief gets the masterId out of a clustered database Id.
 * 
 * if there is no masterId in the given reference the function calls an error.
 */
int GParseMasterId(String<char> & masterId/**<[out]object in which master id is written*/, String<char> & RefSeq/**<[in]id from which masterId is read*/) {
	masterId="";
	int dollarsFound=0;
	for (int i=0; i<length(RefSeq); ++i) {
		if (RefSeq[i]=='$') {
			++dollarsFound;
			continue;
		}
		if (dollarsFound==1)
			appendValue(masterId,RefSeq[i]);
		if (dollarsFound>2)
			cout << "Warning: wrong MasterId has been parsed -> Ids should not contain '$'" << endl ;  
	  
	}
	if(masterId==""){
	  cerr << "no masterId found in " << RefSeq << endl;
	  return 1;}
	return 0;
}


#endif  // GINGER_GSEARCH_GPARSEMASTERID_H_