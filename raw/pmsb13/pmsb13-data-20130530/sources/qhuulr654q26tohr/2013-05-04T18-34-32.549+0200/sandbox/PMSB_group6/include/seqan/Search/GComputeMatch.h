#ifndef GINGER_SEARCH_GCOMPUTEMATCH_H_
#define GINGER_SEARCH_GCOMPUTEMATCH_H_

#include <seqan/structs/GFastaRecord.h>
#include <seqan/structs/GMatch.h>
#include <seqan/Search/GScore.h>
#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence, typename TScoringScheme>
int GComputeMatch(String<GMatch<TScore> > & matches, GFastaRecord<TSequence> & queryRecord, GFastaRecord<TSequence> & clusterRecord, TScoringScheme scoringScheme)
{
	TScore score;
	GScore(score, queryRecord, clusterRecord, scoringScheme);
	GMatch<TScore> match(clusterRecord.id, score);
	appendValue(matches, match);
	cout << match.score << endl;
	return 0;
}

#endif  // GINGER_SEARCH_GCOMPUTEMATCH_H_