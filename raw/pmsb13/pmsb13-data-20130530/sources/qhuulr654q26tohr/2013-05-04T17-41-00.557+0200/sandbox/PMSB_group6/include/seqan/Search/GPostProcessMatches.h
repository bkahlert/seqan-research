#ifndef GINGER_SEARCH_GPOSTPROCESSMATCHES_H_
#define GINGER_SEARCH_GPOSTPROCESSMATCHES_H_

#include <seqan/structs/GMatch.h>
#include <seqan/sequence.h>
#include <algorithm>

using namespace seqan;
using namespace std;

bool myVergleich (GMatch i,GMatch j) {
	return (i.score<j.score);
}

int sortiereMatches (String<GMatch> & unsortetResults) {
	std::sort (begin(unsortetResults), end(unsortetResults), myVergleich);
	return 0;
}

int GComputeMatch(String<GMatch> & matches)
{
	sortiereMatches(matches);
	return 0;
}

#endif  // GINGER_SEARCH_GPOSTPROCESSMATCHES_H_