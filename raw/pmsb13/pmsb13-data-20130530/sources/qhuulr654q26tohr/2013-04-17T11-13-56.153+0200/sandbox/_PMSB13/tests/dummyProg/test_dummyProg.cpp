#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_dummyProg.h"


SEQAN_BEGIN_TESTSUITE(test_dummyProg)
{
    // Call tests.
	SEQAN_CALL_TEST(test_dummyProg_strings_example1);
}
SEQAN_END_TESTSUITE
