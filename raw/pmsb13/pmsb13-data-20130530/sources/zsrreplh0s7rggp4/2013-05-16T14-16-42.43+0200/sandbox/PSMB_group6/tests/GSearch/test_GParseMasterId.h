#ifndef GINGER_TEST_GSEARCH_GPARSEMASTERID_H_
#define GINGER_TEST_GSEARCH_GPARSEMASTERID_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GParseMasterID_dummy)
{
  String<char> masterId= "$ABCD$";
  String<char> refSeq= "ABCDE$CDE$"
 
  
  storage.maxId="CCC";
  clusterRecord.id="CDAB$CCC$";
  
  String<char> res = GParseMasterID(masterId,refSeq);
	
  SEQAN_ASSERT_EQ(res=="$CDE$");
}

#endif  // GINGER_TEST_GPARSEMASTERID_H_ 
