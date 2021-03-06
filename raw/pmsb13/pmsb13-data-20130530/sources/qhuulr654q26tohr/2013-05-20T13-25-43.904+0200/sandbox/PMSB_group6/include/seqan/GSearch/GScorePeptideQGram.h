//Autor:Jakob

#ifndef GINGER_GSEARCH_GSCOREPEPTIDEQGRAM_H_
#define GINGER_GSEARCH_GSCOREPEPTIDEQGRAM_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/index.h>

#include <seqan/GStructs/GFastaRecord.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence, typename TScoringScheme>
int GScorePeptideQGram(TScore &score, GFastaRecord<TSequence> &newRecord, GFastaRecord<TSequence> &clusterRecord,
           TScoringScheme &scoringScheme, string &scoreMode)
{
	// LOCAL_ALIGNMENT_MODE
	if (scoreMode=="LOCAL_ALIGNMENT_MODE") {
		score=0;
		LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, 20);
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

// GLOBAL_ALIGNMENT_MODE
if (scoreMode=="GLOBAL_ALIGNMENT_MODE") {
	Align < String<Dna5> > align;
	resize(rows(align), 2);
	assignSource(row(align, 0), newRecord.seq);
	assignSource(row(align, 1), clusterRecord.seq);
	score = globalAlignment(align,Score<int,Simple>(1,-1,0),Hirschberg());
	return 0;
	}
	
	// COMMON_QGRAM_MODE
    if (scoreMode=="COMMON_QGRAM_MODE") {
        typedef StringSet<TSequence > TSet;
        typedef Index< TSet, IndexQGram<UngappedShape<3> > > TIndex;
        TSet mySet;
        appendValue(mySet, clusterRecord.seq);
        appendValue(mySet, newRecord.seq);
        TIndex myIndex(mySet);
        String<double> distMat;
        getKmerSimilarityMatrix(myIndex,distMat);
        score=(int) (distMat[1]*1000);
        return 0;
    }
    
    // NO_MODE
    cerr << "ERROR: Called GScorePeptideQGram with scoreMode " << scoreMode << endl;
    return 1;
}

#endif  // GINGER_GSEARCH_GSCOREPEPTIDEQGRAM_H_