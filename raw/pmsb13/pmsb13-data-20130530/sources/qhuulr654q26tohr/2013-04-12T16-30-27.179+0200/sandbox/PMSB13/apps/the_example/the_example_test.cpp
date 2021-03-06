#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/file.h>

double quadrat(double x) {
    return x*x;
}

void iota(seqan::String<int> & result, int begin , int end){
    resize(result, end - begin, 0);
    for (int i = begin, k = 0; i<end; ++k, ++i)
	result[k] = i;
}

SEQAN_DEFINE_TEST(test_iota)
{
    seqan::String<int> result;
    iota(result, 0, 3);
    
    SEQAN_ASSERT_EQ(length(result),3);
    SEQAN_ASSERT_EQ(result[0],0);
    SEQAN_ASSERT_EQ(result[1],1);
    SEQAN_ASSERT_EQ(result[2],2);
}

SEQAN_DEFINE_TEST(test_1)
{
    double x = quadrat(3.0);
    SEQAN_ASSERT_IN_DELTA(x,9.0,0.00000001);
}
SEQAN_DEFINE_TEST(test_2)
{
    double x = quadrat(1.2);
    SEQAN_ASSERT_IN_DELTA(x,1.44,0.00000001);
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_1);
    SEQAN_CALL_TEST(test_2);
    SEQAN_CALL_TEST(test_iota);
}
SEQAN_END_TESTSUITE