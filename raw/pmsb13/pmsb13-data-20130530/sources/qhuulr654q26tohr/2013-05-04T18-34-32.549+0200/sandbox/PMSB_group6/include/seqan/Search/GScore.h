#ifndef GINGER_SEARCH_GSCORE_H_
#define GINGER_SEARCH_GSCORE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/structs/GFastaRecord.h>

using namespace seqan;
using namespace std;


template<typename TScore, typename TSequence, typename TScoringScheme>
int GScore(TScore &score, GFastaRecord<TSequence> &newRecord, GFastaRecord<TSequence> &clusterRecord,
				TScoringScheme &scoringScheme)
{
	score=0;
	LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, 1);
	Align < String<Dna5> > align;
	resize(rows(align), 2);
	assignSource(row(align, 0), newRecord.seq);
	assignSource(row(align, 1), clusterRecord.seq);
	
	while (nextLocalAlignment(align, enumerator))
	{
		score +=getScore(enumerator);
	}
	return 0;
}

#endif  // GINGER_SEARCH_GSCORE_H_