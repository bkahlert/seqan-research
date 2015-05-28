//Autor:Hannes
#ifndef GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_
#define GINGER_GSEARCH_GCHECKCLUSTERMASTER_H_

#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/GStructs/GScoreStorage.h>
#include <seqan/GSearch/GParseMasterId.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence>
/**
 * \brief checks if the two given masterIds are the same
 * 
 * The function compares the master sequences out of the two goven objects and returns weather the master sequences are the same. 
 * The function creates an error if there are empty master Ids.
 */

bool GCheckClusterMaster(GScoreStorage<TScore> & storage/**<[in] object in which the first master Id is stored */, GFastaRecord<TSequence> & clusterRecord/**<[in] the record which contains the second MasterId */, bool MASTER_FILE_GIVEN/**<[in] contains weather the cluster record contains a master or not */)
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