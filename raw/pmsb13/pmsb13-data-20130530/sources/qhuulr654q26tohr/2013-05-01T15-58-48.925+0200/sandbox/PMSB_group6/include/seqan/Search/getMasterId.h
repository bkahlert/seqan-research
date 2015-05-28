#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_GETMASTERID_H_
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_GETMASTERID_H_

#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

CharString getMasterId(String<char> & RefSeq) {
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


#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_GETMASTERID_H_