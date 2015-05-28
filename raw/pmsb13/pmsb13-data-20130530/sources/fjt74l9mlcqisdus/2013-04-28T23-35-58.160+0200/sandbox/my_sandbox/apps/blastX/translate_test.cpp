#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1
#include "own_functions.h"

int hash(int x,int y,int z){
	int erg = (x*16)+(y*4)+z;
	if (erg>=0 && erg<=63) return erg;
}


SEQAN_DEFINE_TEST(test_my_app_hash_normal)
{
	int result = hash(1,0,3);
	SEQAN_ASSERT_IN_DELTA(result,31.5,31.5);
}

SEQAN_DEFINE_TEST(test_my_app_hash_zu_hoch)
{
	int result = hash(4,7,900);
	SEQAN_ASSERT_IN_DELTA(result,31.5,31.5);
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_hash_normal);
	SEQAN_CALL_TEST(test_my_app_hash_zu_hoch);
}
SEQAN_END_TESTSUITE