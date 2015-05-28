#ifndef GINGER_TEST_GSEARCH_GPOSTPROCESSMATCHES_H_
#define GINGER_TEST_GSEARCH_GPOSTPROCESSMATCHES_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GPostProcessMatches_dummy)
{
 String<GMatch<int> > unsortetResults;
 unsortetResults[0].score=2;
 unsortetResults[1].score=3;
 unsortetResults[2].score=7;
 unsortetResults[3].score=13;
 unsortetResults[4].score=5;
 
 int res=GPostProcessMatches(unsortetResults);
 SEQAN_ASSERT_EQ(unsortetResults[0].score,2);
 SEQAN_ASSERT_EQ(unsortetResults[1].score,3);
 SEQAN_ASSERT_EQ(unsortetResults[2].score,5);
 SEQAN_ASSERT_EQ(unsortetResults[3].score,7);
 SEQAN_ASSERT_EQ(unsortetResults[4].score,13);
 SEQAN_ASSERT_EQ(res,0);
  
}

#endif  // GINGER_TEST_GPOSTPROCESSMATCHES_H_ 
 
