#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1
#include "own_functions.h"

int hash(int x,int y,int z){
	int erg = (x*16)+(y*4)+z;
	return erg;
}


SEQAN_DEFINE_TEST(test_my_app_hash)
{
	int result = hash(1,0,3);
	SEQAN_ASSERT_IN_DELTA(result,31.5,31.5);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_hash);
}
SEQAN_END_TESTSUITE