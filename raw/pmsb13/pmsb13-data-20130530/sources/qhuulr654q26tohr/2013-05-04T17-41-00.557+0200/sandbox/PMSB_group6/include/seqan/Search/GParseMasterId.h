#ifndef GINGER_SEARCH_GPARSEMASTERID_H_
#define GINGER_SEARCH_GPARSEMASTERID_H_

#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

CharString GParseMasterId(String<char> & RefSeq) {
	CharString MasterId;
	int dollarsFound=0;
	for (int i=0; i<length(RefSeq); ++i) {
		if (RefSeq[i]=='$') {
			++dollarsFound;
			continue;
		}
		if (dollarsFound==1)
			appendValue(MasterId,RefSeq[i]);
	}
	return MasterId;
}


#endif  // GINGER_SEARCH_GPARSEMASTERID_H_