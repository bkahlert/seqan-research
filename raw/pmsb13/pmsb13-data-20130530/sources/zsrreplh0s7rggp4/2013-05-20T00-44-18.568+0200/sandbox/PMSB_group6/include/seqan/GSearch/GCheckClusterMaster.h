//Autor:Hannes
#ifndef GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_
#define GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_

#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/GStructs/GScoreStorage.h>
#include <seqan/GSearch/GParseMasterId.h>

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
	  //cout << storage.maxId << " == " << masterId << endl;
	return (storage.maxId == masterId);
}

#endif  // GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_