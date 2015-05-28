#ifndef GINGER_SEARCH_GCOMPUTEMATCH_H_
#define GINGER_SEARCH_GCOMPUTEMATCH_H_

#include <seqan/structs/GFastaRecord.h>
#include <seqan/structs/GMatch.h>
#include <seqan/Search/GScore.h>
#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

template<typename TScore>
int GComputeMatch(String<GMatch> & matches, GFastaRecord & queryRecord, GFastaRecord & clusterRecord)
{
	TScore score;
	GScore(score, queryRecord, clusterRecord, scoringScheme);
	GMatch match(clusterRecord.id, score);
	appendValue(matches, match);
	return 0;
}

#endif  // GINGER_SEARCH_GCOMPUTEMATCH_H_