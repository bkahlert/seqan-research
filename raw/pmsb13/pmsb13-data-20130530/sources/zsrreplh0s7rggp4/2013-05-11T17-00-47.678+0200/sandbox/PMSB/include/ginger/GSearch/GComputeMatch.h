#ifndef GINGER_GSEARCH_GCOMPUTEMATCH_H_
#define GINGER_GSEARCH_GCOMPUTEMATCH_H_

#include <ginger/GStructs/GFastaRecord.h>
#include <ginger/GStructs/GMatch.h>
#include <ginger/GSearch/GScore.h>
#include <ginger/sequence.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence, typename TScoringScheme>
int GComputeMatch(String<GMatch<TScore> > & matches, GFastaRecord<TSequence> & queryRecord, GFastaRecord<TSequence> & clusterRecord, TScoringScheme & scoringScheme, string & scoreMode)
{
	TScore score;
	GScore(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	GMatch<TScore> match(clusterRecord.id, score);
	appendValue(matches, match);
	return 0;
}

#endif  // GINGER_GSEARCH_GCOMPUTEMATCH_H_