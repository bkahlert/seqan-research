#ifndef GINGER_SEARCH_GCHECKCLUSTERMASTER_H_
#define GINGER_SEARCH_GCHECKCLUSTERMASTER_H_

#include <seqan/structs/GFastaRecord.h>
#include <seqan/structs/GScoreStorage.h>
#include <seqan/Search/GParseMasterId.h>

using namespace seqan;
using namespace std;

bool GCheckClusterMaster(GScoreStorage & storage, GFastaRecord & clusterRecord)
{
	CharString masterId;
	GParseMasterId(masterId, clusterRecord.id);
    return (storage.maxId == masterId);
}

#endif  // GINGER_SEARCH_GCHECKCLUSTERMASTER_H_