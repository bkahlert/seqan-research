//Autor:Hannes(local und global aligment)/Jakob(q-gram)

#ifndef GINGER_GSEARCH_GSCORE_H_
#define GINGER_GSEARCH_GSCORE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/index.h>

#include <seqan/GStructs/GFastaRecord.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence, typename TScoringScheme>
/**
 * \brief scoring function for DNA and RNA mode
 * This function can compute a similarity score of the two sequences of the two fasta records in three different ways.
 * It is strongly recommended to use COMMON_QGRAM_MODE because the other two modi are significantly slower and were mainly used as a first approach to scoring.
 * LOCAL_ALIGNMENT_MODE:
 * Computes all local alignments whose score is 20 or higher and sums up those scores.
 * GLOBAL_ALIGNMENT_MODE:
 * Computes a global alignment with match benefit 1, mismatch penalty -1 and no gap costs.
 * COMMON_QGRAM_MODE:
 * Counts all 8-grams in the two sequences. The fraction of 8-grams from the smaller sequence which also occur in the bigger sequence is muliplied with 1000 and stored as an integer score value.
 * 
 * Please note that this function is only called on DnaString and RnaString. For Peptide see function GScorePeptideQGram.
 */
int GScore(TScore &score/**<[out] computed score*/, GFastaRecord<TSequence> &newRecord/**<[in] one fasta record*/, GFastaRecord<TSequence> &clusterRecord/**<[in] another fasta record*/,
           TScoringScheme &scoringScheme/**<[in] scoring scheme used for LOCAL_ALIGNMENT_MODE*/, string &scoreMode/**<[in] scoring mode to use*/)
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
        typedef Index< TSet, IndexQGram<UngappedShape<8> > > TIndex;
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
    cerr << "Unknown score mode: " << scoreMode << endl;
    return 1;
}

#endif  // GINGER_GSEARCH_GSCORE_H_