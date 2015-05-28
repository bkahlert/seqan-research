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
  storage.scoreInitialized=true;
  GFastaRecord<Dna5> referenceRecord;
  referenceRecord.id="NUMBER1";
  storage.maxId="NUMBER2";
  storage.maxScore=3;

  int res=GEvaluateScore(storage,score,referenceRecord);
  
  SEQAN_ASSERT_EQ(res,0);
  SEQAN_ASSERT_EQ(storage.maxId,"NUMBER1");
  SEQAN_ASSERT_EQ(storage.maxScore,5);
}

SEQAN_DEFINE_TEST(test_GEvaluateScore_sameScore)
{
  GScoreStorage<int> storage;
  int score=5;
  storage.scoreInitialized=true;
  GFastaRecord<Dna5> referenceRecord;
  referenceRecord.id="NUMBER1";
  storage.maxId="NUMBER2";
  storage.maxScore=5;

  int res=GEvaluateScore(storage,score,referenceRecord);
  
  SEQAN_ASSERT_EQ(res,0);
  SEQAN_ASSERT_EQ(storage.maxId,"NUMBER1");
  SEQAN_ASSERT_EQ(storage.maxScore,5);
}

SEQAN_DEFINE_TEST(test_GEvaluateScore_sameId)
{
  GScoreStorage<int> storage;
  int score=5;
  storage.scoreInitialized=true;
  GFastaRecord<Dna5> referenceRecord;
  referenceRecord.id="NUMBER1";
  storage.maxId="NUMBER1";
  storage.maxScore=3;

  int res=GEvaluateScore(storage,score,referenceRecord);
  
  SEQAN_ASSERT_EQ(res,0);
  SEQAN_ASSERT_EQ(storage.maxId,"NUMBER1");
  SEQAN_ASSERT_EQ(storage.maxScore,5);
}

SEQAN_DEFINE_TEST(test_GEvaluateScore_noId)
{
  GScoreStorage<int> storage;
  int score=5;
  storage.scoreInitialized=true;
  GFastaRecord<Dna5> referenceRecord;
  referenceRecord.id="";
  storage.maxId="CAATTGG";
  storage.maxScore=3;

  int res=GEvaluateScore(storage,score,referenceRecord);
  
  SEQAN_ASSERT_EQ(res,1);
  SEQAN_ASSERT_EQ(storage.maxId,"");
  SEQAN_ASSERT_EQ(storage.maxScore,5);
}

SEQAN_DEFINE_TEST(test_GEvaluateScore_dummy2)
{
  
  GScoreStorage<int> storage;
  storage.scoreInitialized=true;
  int score=3;
  GFastaRecord<Dna5> referenceRecord;
  referenceRecord.id="NUMBER1";
  storage.maxId="NUMBER2";
  storage.maxScore=5;

  int res=GEvaluateScore(storage,score,referenceRecord);
  
  SEQAN_ASSERT_EQ(res,0);
  SEQAN_ASSERT_EQ(storage.maxId,"NUMBER2");
  SEQAN_ASSERT_EQ(storage.maxScore,5);
}

SEQAN_DEFINE_TEST(test_GEvaluateScore_noScore)
{
  
  GScoreStorage<int> storage;
  int score=3;
  GFastaRecord<Dna5> referenceRecord;
  referenceRecord.id="NUMBER1";
  storage.maxId="NUMBER2";
  storage.maxScore=5;

  int res=GEvaluateScore(storage,score,referenceRecord);
  
  SEQAN_ASSERT_EQ(res,0);
  SEQAN_ASSERT_EQ(storage.maxId,"NUMBER1");
  SEQAN_ASSERT_EQ(storage.maxScore,3);
}
#endif  // GINGER_TEST_GEVALUATESCORE_H_ 
