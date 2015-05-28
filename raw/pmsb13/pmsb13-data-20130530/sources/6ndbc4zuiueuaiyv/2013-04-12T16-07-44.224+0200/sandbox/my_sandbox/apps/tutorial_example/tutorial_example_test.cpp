#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>
//#include <impl.h>


SEQAN_DEFINE_TEST(test_my_app_funcs_hello)
{
    double x = quadrat(3.0);
    SEQAN_ASSERT_EQ(x, 9.0);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_funcs_hello);
}
SEQAN_END_TESTSUITE

