#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_GEvaluateScore_H_
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_GEvaluateScore_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>

using namespace seqan;
using namespace std;


template<typename TScore>
int GScore(GScoreStorage & storage, TScore & score, GFastaRecord & referenceRecord)
{
	if (score >= storage.maxScore) {
		storage.maxScore=score;
		storage.maxId=referenceRecord.id;
	}
	return 0;
}

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_GEvaluateScore_H_