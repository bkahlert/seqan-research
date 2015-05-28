#ifndef GINGER_TEST_GSEARCH_GCHECKCLUSTERMASTER_H_
#define GINGER_TEST_GSEARCH_GCHECKCLUSTERMASTER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GClusterMaster.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GCheckClusterMaster_empty)
{
  GScoreStorage<TScore> storage;
  GFastaRecord<TSequence> clusterRecord;
  bool MASTER_FILE_GIVEN=TRUE;
  
  storage.maxId="ABCD$CCC$";
  clusterRecord.id="CDAB$CCC$";
  
  int res = GCheckClusterMaster(storage,clusterRecord,MASTER_FILE_GIVEN);
	
  SEQAN_ASSERT_EQ(res,0);
}
  