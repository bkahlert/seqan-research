//Autor:Jakob

#ifndef GINGER_GSEARCH_GEVALUATESCORE_H_
#define GINGER_GSEARCH_GEVALUATESCORE_H_

#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/GStructs/GScoreStorage.h>

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence>
/** \brief stores maximum score and corresponding fasta record id.
 * The function checks the given score whether it is higher than the current maximum score.
 * If so it updates the maximum score and the corresponding (given) sequence id.
 * These information are stored in a GScoreStorage object.
 */
int GEvaluateScore(GScoreStorage<TScore> & storage/**<[out]stores maxScore and corresponding id*/, TScore & score/**<[score which is compared to current maximum]*/, GFastaRecord<TSequence> & referenceRecord/**<[fasta record from the masterFile, which caused the score]*/)
{
    if (storage.scoreInitialized) {  
      if (score > storage.maxScore) {
            storage.maxScore=score;
            storage.maxId=referenceRecord.id;
        }
    }
    else {
        storage.maxScore=score;
        storage.maxId=referenceRecord.id;
		storage.scoreInitialized=true;
    }
    if (storage.maxId=="")
      {cerr << "Warning: empty storage.maxId... may lead to empty results!" << endl;
	return 1;}
      return 0;
}

#endif  // GINGER_GSEARCH_GEVALUATESCORE_H_