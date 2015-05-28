//Autor:Hannes

#ifndef GINGER_GSEARCH_GCOMPUTEMATCH_H_
#define GINGER_GSEARCH_GCOMPUTEMATCH_H_

#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/GStructs/GMatch.h>
#include <seqan/GSearch/GScore.h>
#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence, typename TScoringScheme>
int GComputeMatch(String<GMatch<TScore, TSequence> > & matches, GFastaRecord<TSequence> & queryRecord, GFastaRecord<TSequence> & clusterRecord, TScoringScheme & scoringScheme, TScore threshold, string & scoreMode)
{
	TScore score;
	GScore(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	
	if (score >= threshold){
		GMatch<TScore, TSequence> match(clusterRecord.id, score);
		match.queryId = queryRecord.id;
		match.querySequence=queryRecord.seq;
		match.targetSequence=clusterRecord.seq;
		appendValue(matches, match);
	}
	return 0;
}

#endif  // GINGER_GSEARCH_GCOMPUTEMATCH_H_