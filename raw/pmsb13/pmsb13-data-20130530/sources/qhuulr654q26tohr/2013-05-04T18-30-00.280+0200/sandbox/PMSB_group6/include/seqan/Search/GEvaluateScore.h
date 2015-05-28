#ifndef GINGER_SEARCH_GEVALUATESCORE_H_
#define GINGER_SEARCH_GEVALUATESCORE_H_

#include <seqan/structs/GFastaRecord.h>
#include <seqan/structs/GScoreStorage.h>

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence>
int GEvaluateScore(GScoreStorage<TScore> & storage, TScore & score, GFastaRecord<TSequence> & referenceRecord)
{
    if (storage.scoreInitialized) {
        if (score >= storage.maxScore) {
            storage.maxScore=score;
            storage.maxId=referenceRecord.id;
        }
    }
    else {
        storage.maxScore=score;
        storage.maxId=referenceRecord.id;
    }
    return 0;
}

#endif  // GINGER_SEARCH_GEVALUATESCORE_H_