#ifndef GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_
#define GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_

#include <seqan/GStructs/GMatch.h>
#include <seqan/sequence.h>
#include <algorithm>

using namespace seqan;
using namespace std;

template<typename TScore>
bool myVergleich (GMatch<TScore> const & i,GMatch<TScore> const & j) {
	return (i.score >= j.score);
}

template<typename TScore>
int sortiereMatches (String<GMatch<TScore> > & unsortetResults) {
	std::sort (begin(unsortetResults), end(unsortetResults), myVergleich);
	return 0;
}

template<typename TScore>
int GPostProcessMatches(String<GMatch<TScore> > & matches)
{
	//sort (begin(matches), end(matches), myVergleich);
	sortiereMatches(matches);
	return 0;
}

#endif  // GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_