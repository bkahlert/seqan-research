#ifndef GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_
#define GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_

#include <seqan/structs/GFastaRecord.h>
#include <seqan/structs/GScoreStorage.h>
#include <seqan/Search/GParseMasterId.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence>
bool GCheckClusterMaster(GScoreStorage<TScore> & storage, GFastaRecord<TSequence> & clusterRecord)
{
	CharString masterId;
	GParseMasterId(masterId, clusterRecord.id);
    return (storage.maxId == masterId);
}

#endif  // GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_