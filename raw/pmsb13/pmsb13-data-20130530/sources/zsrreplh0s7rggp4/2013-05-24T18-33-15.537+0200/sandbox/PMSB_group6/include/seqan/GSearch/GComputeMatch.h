//Autor:Hannes

#ifndef GINGER_GSEARCH_GCOMPUTEMATCH_H_
#define GINGER_GSEARCH_GCOMPUTEMATCH_H_

#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/GStructs/GMatch.h>
#include <seqan/GSearch/GScore.h>
#include <seqan/GSearch/GScorePeptideQGram.h>
#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence, typename TScoringScheme>
/**
 * \brief Computes the Match of a query and a clustered record
 * 
 * This function Computes the Match of a queryRecord and a clusterRecord, using either Peptide or RNA/DNA and therefore using either the GScorePeptideQGram or the GScore function.
 * 
 */
int GComputeMatch(String<GMatch<TScore, TSequence> > & matches/**<[out]object in which the match is written*/, GFastaRecord<TSequence> & queryRecord/**<[in]the query*/, GFastaRecord<TSequence> & clusterRecord/**<[in]the entry of the clustered database*/, TScoringScheme & scoringScheme/**<[in]the used scoring scheme*/, TScore threshold/**<[in]the threshold which determines if a match is valid*/, string & scoreMode/**<[in]the Score mode used by GScore*/)
{
	TScore score;
	if (IsSameType<TSequence, Peptide>::VALUE){
		GScorePeptideQGram(score, queryRecord, clusterRecord, scoringScheme, scoreMode);}
	else
		GScore(score, queryRecord, clusterRecord, scoringScheme, scoreMode);
	
	if (score >= threshold){
		GMatch<TScore, TSequence> match(clusterRecord.id, score);
		match.queryId = queryRecord.id;
		match.querySequence=queryRecord.seq;
		match.targetSequence=clusterRecord.seq;
		appendValue(matches, match);
	}
	return 0;
}

#endif  // GINGER_GSEARCH_GCOMPUTEMATCH_H_