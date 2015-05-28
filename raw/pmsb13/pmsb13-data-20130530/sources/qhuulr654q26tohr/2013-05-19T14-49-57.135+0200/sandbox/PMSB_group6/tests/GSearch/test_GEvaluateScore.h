#ifndef GINGER_TEST_GSEARCH_GEVALUATESCORE_H_
#define GINGER_TEST_GSEARCH_GEVALUATESCORE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GEvaluateScore_dummy)
{
  GScoreStorage<int> storage;
  int score=5;
  GFastaRecord<Dna5> referenceRecord;
  referenceRecord.seq="ATGGCTAGAAGCT";
  referenceRecord.id="NUMBER1";
  storage.maxId="NUMBER2";
  storage.maxScore=3;

  int res=GEvaluateScore(storage,score,referenceRecord);
  
  SEQAN_ASSERT_EQ(res,0);
  SEQAN_ASSERT_EQ(storage.maxId,"NUMBER1");
  SEQAN_ASSERT_EQ(storage.maxScore,5);
}

#endif  // GINGER_TEST_GEVALUATESCORE_H_ 
