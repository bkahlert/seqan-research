#ifndef GINGER_TEST_GSEARCH_GCHECKCLUSTERMASTER_H_
#define GINGER_TEST_GSEARCH_GCHECKCLUSTERMASTER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GCheckClusterMaster_dummy)
{
  GScoreStorage<int> storage;
  GFastaRecord<CharString> clusterRecord;
  bool MASTER_FILE_GIVEN=true;
  
  storage.maxId="CCC";
  clusterRecord.id="CDAB$CCC$";
  
  bool res = GCheckClusterMaster(storage,clusterRecord,MASTER_FILE_GIVEN);
	
  SEQAN_ASSERT(res);
}


SEQAN_DEFINE_TEST(test_GCheckClusterMaster_halfMaster)
{
  GScoreStorage<int> storage;
  GFastaRecord<CharString> clusterRecord;
  bool MASTER_FILE_GIVEN=true;
  
  storage.maxId="ABCD$CCC$";
  clusterRecord.id="CDDA$ABCD$";
  
  bool res = GCheckClusterMaster(storage,clusterRecord,MASTER_FILE_GIVEN);
	
  SEQAN_ASSERT(res);
}

SEQAN_DEFINE_TEST(test_GCheckClusterMaster_withoutMaster)
{
  GScoreStorage<int> storage;
  GFastaRecord<CharString> clusterRecord;
  bool MASTER_FILE_GIVEN=true;
  
  storage.maxId="ABCD";
  clusterRecord.id="ABCD";
  
  bool res = GCheckClusterMaster(storage,clusterRecord,MASTER_FILE_GIVEN);
	
   SEQAN_ASSERT(!res);
}

SEQAN_DEFINE_TEST(test_GCheckClusterMaster_empty)
{
  GScoreStorage<int> storage;
  GFastaRecord<CharString> clusterRecord;
  bool MASTER_FILE_GIVEN=true;
  
  storage.maxId="";
  clusterRecord.id="";
  
  bool res = GCheckClusterMaster(storage,clusterRecord,MASTER_FILE_GIVEN);
	
   SEQAN_ASSERT(!res);
}
#endif  // GINGER_TEST_GCHECKCLUSTERMASTER_H_