#ifndef GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_
#define GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_

#include "../GStructs/GFastaRecord.h"
#include "../GStructs/GScoreStorage.h"
#include "GParseMasterId.h"

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence>
bool GCheckClusterMaster(GScoreStorage<TScore> & storage, GFastaRecord<TSequence> & clusterRecord, bool MASTER_FILE_GIVEN)
{
	if (!MASTER_FILE_GIVEN)
		return true;
	CharString masterId;
	GParseMasterId(masterId, clusterRecord.id);
	if (masterId==""){
	  cerr << "masterId is empty"<< endl;
	  return false;}
	if (storage.maxId==""){
	  cerr << "storage.maxId is empty"<< endl;
	  return false;}  
	return (storage.maxId == masterId);
}

#endif  // GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_