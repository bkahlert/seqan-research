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
  int score;
  GFastaRecord<Dna5> referenceRecord;

  res=GEvaluateScore(storage,score,referenceRecord);
  
  SEQAN_ASSERT_EQ(res,)
}

#endif  // GINGER_TEST_GEVALUATESCORE_H_ 
