#ifndef GINGER_SEARCH_GPARSEMASTERID_H_
#define GINGER_SEARCH_GPARSEMASTERID_H_

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
	}
	return 0;
}


#endif  // GINGER_SEARCH_GPARSEMASTERID_H_