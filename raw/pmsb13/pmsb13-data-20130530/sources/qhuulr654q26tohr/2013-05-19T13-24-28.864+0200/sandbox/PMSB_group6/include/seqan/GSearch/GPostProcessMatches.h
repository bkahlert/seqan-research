#ifndef GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_
#define GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_

#include <seqan/GStructs/GMatch.h>
#include <seqan/sequence.h>
#include <algorithm>

using namespace seqan;
using namespace std;

bool myVergleich (GMatch<int> const & i,GMatch<int> const & j) {
	return (i.score >= j.score);
}

int sortiereMatches (String<GMatch<int> > & unsortetResults) {
	std::sort (begin(unsortetResults), end(unsortetResults), myVergleich);
	return 0;
}

template<typename TScore>
int GPostProcessMatches(String<GMatch<TScore> > & matches)
{
	sort (begin(unsortetResults), end(unsortetResults), myVergleich);
	//sortiereMatches(matches);
	return 0;
}

#endif  // GINGER_GSEARCH_GPOSTPROCESSMATCHES_H_