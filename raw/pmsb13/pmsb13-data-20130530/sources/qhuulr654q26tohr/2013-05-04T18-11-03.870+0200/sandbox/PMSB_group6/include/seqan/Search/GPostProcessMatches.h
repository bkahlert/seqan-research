#ifndef GINGER_SEARCH_GPOSTPROCESSMATCHES_H_
#define GINGER_SEARCH_GPOSTPROCESSMATCHES_H_

#include <seqan/structs/GMatch.h>
#include <seqan/sequence.h>
#include <algorithm>

using namespace seqan;
using namespace std;

template<typename TScore>
bool myVergleich (GMatch<TScore> i,GMatch<TScore> j) {
	return (i.score<TScore> < j.score<TScore>);
}

template<typename TScore>
int sortiereMatches (String<GMatch<TScore> > & unsortetResults) {
	std::sort (begin(unsortetResults), end(unsortetResults), myVergleich);
	return 0;
}

template<typename TScore>
int GPostProcessMatches(String<GMatch<TScore> > & matches)
{
	sortiereMatches(matches);
	return 0;
}

#endif  // GINGER_SEARCH_GPOSTPROCESSMATCHES_H_