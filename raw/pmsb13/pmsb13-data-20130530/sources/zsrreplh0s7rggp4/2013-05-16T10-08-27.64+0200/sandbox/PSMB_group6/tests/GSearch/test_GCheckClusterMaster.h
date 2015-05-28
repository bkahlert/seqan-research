#ifndef GINGER_TEST_GSEARCH_GCHECKCLUSTERMASTER_H_
#define GINGER_TEST_GSEARCH_GCHECKCLUSTERMASTER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GCheckClusterMaster_empty)
{
  GScoreStorage<int> storage;
  GFastaRecord<CharString> clusterRecord;
  bool MASTER_FILE_GIVEN=true;
  
  storage.maxId="ABCD$CCC$";
  clusterRecord.id="CDAB$CCC$";
  
  int res = GCheckClusterMaster(storage,clusterRecord,MASTER_FILE_GIVEN);
	
  SEQAN_ASSERT_EQ(res,0);
}

#endif  // GINGER_TEST_GCHECKCLUSTERMASTER_H_