#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>

double quadrat(double x) {
    return x*x;
}

SEQAN_DEFINE_TEST(test_1)
{
    double x = quadrat(3.0);
    SEQAN_ASSERT_EQ(x,9.0);
}
SEQAN_DEFINE_TEST(test_2)
{
    double x = quadrat(1.2);
    SEQAN_ASSERT_EQ(x,1.44);
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_1);
    SEQAN_CALL_TEST(test_2);
}
SEQAN_END_TESTSUITE