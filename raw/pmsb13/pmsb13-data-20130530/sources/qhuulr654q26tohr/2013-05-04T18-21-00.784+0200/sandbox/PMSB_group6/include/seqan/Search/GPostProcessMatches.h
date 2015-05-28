#ifndef GINGER_SEARCH_GPOSTPROCESSMATCHES_H_
#define GINGER_SEARCH_GPOSTPROCESSMATCHES_H_

#include <seqan/structs/GMatch.h>
#include <seqan/sequence.h>
#include <algorithm>

using namespace seqan;
using namespace std;

bool myVergleich (GMatch<int> i,GMatch<int> j) {
	return (i.score < j.score);
}

int sortiereMatches (String<GMatch<int> > & unsortetResults) {
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