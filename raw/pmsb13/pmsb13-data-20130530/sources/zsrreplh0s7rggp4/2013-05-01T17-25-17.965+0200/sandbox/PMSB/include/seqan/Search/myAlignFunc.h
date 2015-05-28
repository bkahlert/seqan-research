#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_MYALIGNFUNC_H_
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_MYALIGNFUNC_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>

using namespace seqan;
using namespace std;


template<typename TScore, typename TRecord, typename TScoringScheme>
int myAlignFunc(TScore &score, TRecord &newRecord, TRecord &clusterRecord,
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

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_MYALIGNFUNC_H_