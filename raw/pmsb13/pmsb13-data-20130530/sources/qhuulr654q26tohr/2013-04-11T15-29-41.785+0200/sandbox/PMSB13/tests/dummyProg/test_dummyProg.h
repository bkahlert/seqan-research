#ifndef SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_
#define SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <../../sandbox/PMSB13/apps/dummyprog/dummyprog.cpp>

// A test for strings.
SEQAN_DEFINE_TEST(test_dummyProg_strings_example1)
{
    using namespace seqan;
    using namespace std;

    SEQAN_ASSERT(true);
}

#endif  // SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_
