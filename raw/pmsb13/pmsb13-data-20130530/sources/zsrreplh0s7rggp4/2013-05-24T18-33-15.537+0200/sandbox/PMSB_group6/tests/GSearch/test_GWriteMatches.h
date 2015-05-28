//autor:Hannes

#ifndef GINGER_TEST_GSEARCH_GWRITEMATCHES_H_
#define GINGER_TEST_GSEARCH_GWRITEMATCHES_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/GSearch.h>
#include <seqan/GSearch/GWriteMatches.h>


using namespace seqan;
using namespace std;

SEQAN_DEFINE_TEST(test_GWriteMatches_dummy)
{
  String<GMatch<int, DnaString> > matches;
  GFastaRecord<DnaString> queryRecord;
  queryRecord.id="testQueryId";
  queryRecord.seq="TTTATTGCTAAGCCAAA";
  resize(matches,3);
  matches[0].id="ABC";
  matches[1].id="DEF";
  matches[2].id="GHI";
  matches[0].querySequence="GGGGGGG";
  matches[1].querySequence="AAAAAAA";
  matches[2].querySequence="TTTTTTT";
  matches[0].targetSequence="GGAGGAG";
  matches[1].targetSequence="AACAACA";
  matches[2].targetSequence="TTCTTCT";
  matches[0].queryId="JKL";
  matches[1].queryId="MNO";
  matches[2].queryId="PQR";
  matches[0].score=1;
  matches[1].score=2;
  matches[2].score=3;
  fstream outputStream(toCString("../../seqan-trunk/sandbox/PMSB_group6/tests/autoTestData/mutable/outFileGWriteMatches.tsv"));
  
  int res=GWriteMatches(outputStream, matches, queryRecord);
  SEQAN_ASSERT_EQ(res,0);
  SEQAN_ASSERT(outputStream.good());
}

#endif  // GINGER_TEST_GWRITEMATCHES_H_ 

