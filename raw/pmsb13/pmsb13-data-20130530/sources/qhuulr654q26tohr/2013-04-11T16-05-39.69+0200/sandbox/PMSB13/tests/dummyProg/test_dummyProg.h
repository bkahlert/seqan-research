#ifndef SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_
#define SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include </home/development/seqan-trunk/sandbox/PMSB13/include/seqan/dummyProg.h>

using namespace seqan;
using namespace std;

typedef Iterator<Index<String<char>, IndexEsa<> >, TopDown<ParentLinks<> > >::Type TIteratorIndex;


TIteratorIndex getIt(){
    String<char> haystack = "mississippi";
    TIndex index(haystack);
    TIteratorIndex it(index);
    return it;
}

// A test for strings.
SEQAN_DEFINE_TEST(test_dummyProg_strings_example1)
{
    TIteratorIndex it = getIt();
    recGoDown(it);
    SEQAN_ASSERT(true);
}

#endif  // SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_
