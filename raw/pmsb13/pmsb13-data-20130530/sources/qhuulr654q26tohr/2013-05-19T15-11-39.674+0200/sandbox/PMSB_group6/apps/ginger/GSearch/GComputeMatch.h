#ifndef GINGER_GSEARCH_GCOMPUTEMATCH_H_
#define GINGER_GSEARCH_GCOMPUTEMATCH_H_

#include "GStructs/GFastaRecord.h"
#include "GStructs/GMatch.h"
#include "GSearch/GScore.h"

#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence, typename TScoringScheme>
int GComputeMatch(String<GMatch<TScore> > & matches, GFastaRecord<TSequence> & queryRecord, GFastaRecord<TSequence> & clusterRecord, TScoringScheme & scoringScheme, TScore threshold, string & scoreMode)
{
	TScore score;
	GScore(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	if (score >= threshold){
		GMatch<TScore> match(clusterRecord.id, score);
		match.queryId = queryRecord.id;
		appendValue(matches, match);
	}
	return 0;
}

#endif  // GINGER_GSEARCH_GCOMPUTEMATCH_H_