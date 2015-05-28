#ifndef GINGER_GSEARCH_GPARSEMASTERID_H_
#define GINGER_GSEARCH_GPARSEMASTERID_H_

#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

int GParseMasterId(String<char> & masterId, String<char> & RefSeq) {
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
	  cerr << "no masterId found in " RefSeq << endl;
	  return 1;}
	return 0;
}


#endif  // GINGER_GSEARCH_GPARSEMASTERID_H_