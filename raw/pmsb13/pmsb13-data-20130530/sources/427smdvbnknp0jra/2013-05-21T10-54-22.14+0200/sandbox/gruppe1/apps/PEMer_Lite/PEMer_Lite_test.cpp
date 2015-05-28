#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include "PEMer_Lite.h"

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input)
{
	candidates test;
	modifier opt;
	statisticals stats
	char *file ="/Informatik/Development/test.sam";
	SEQAN_ASSERT(!saminput(test,opt,stats));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input2)
{
	candidate test;
	char *file =" ";
	SEQAN_ASSERT(saminput(test,file));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input3)
{
	candidate test;
	char *file ="2579816ß15362466183b513ü6ß148z5!§%%!&°%(!%§(%§c465cn4ü457346511425129562e85z148z5!§%%!&°%(!%§(%§&&&&&(";
	SEQAN_ASSERT(saminput(test,file));
}
/*
SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_devide)
{
	candidates test,result;
	modifier opt;
	for(unsigned i=1;i<11;i++)
	{
		test.len.push_back(i);
		test.ref.push_back(i);
		test.pos.push_back(i);
		test.type.push_back(candidate::u);
	}
	opt.e=5;
	opt.sd=2;
	devide(test,result,opt);

	SEQAN_ASSERT_EQ(result.type[0],1);
	SEQAN_ASSERT_EQ(result.ref[0],1);
	SEQAN_ASSERT_EQ(result.type[1],1);
	SEQAN_ASSERT_EQ(result.ref[1],2);
	SEQAN_ASSERT_EQ(result.type[2],2);
	SEQAN_ASSERT_EQ(result.ref[2],8);
	SEQAN_ASSERT_EQ(result.type[3],2);
	SEQAN_ASSERT_EQ(result.ref[3],9);
	SEQAN_ASSERT_EQ(result.type[4],2);
	SEQAN_ASSERT_EQ(result.ref[4],10);

}
*/
SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_median)
{
	medianTree median;
	for(unsigned i=0;i<100;i++)
	{
		median.add(i);
	}
	SEQAN_ASSERT_EQ(median.getMedian(),50);
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_clustering1)
{
    candidates c(5);

    c[0].addvalues(0,159,417,candidate::insertion);
    c[1].addvalues(1,206,446,candidate::insertion);
    c[2].addvalues(2,836,324,candidate::deletion);
    c[3].addvalues(3,846,552,candidate::deletion);
    c[4].addvalues(4,800,477,candidate::deletion);

    std::vector<std::vector<unsigned> > cluster;

    findCluster(c, cluster, candidate::insertion, 0.5);
    findCluster(c, cluster, candidate::deletion, 0.5);

    SEQAN_ASSERT_EQ(cluster.size(),2);
    SEQAN_ASSERT_EQ(cluster[0].size(),2);
    SEQAN_ASSERT_EQ(cluster[1].size(),3);
    SEQAN_ASSERT_EQ(cluster[0][0],0);
    SEQAN_ASSERT_EQ(cluster[0][1],1);
    SEQAN_ASSERT_EQ(cluster[1][0],2);
    SEQAN_ASSERT_EQ(cluster[1][1],3);
    SEQAN_ASSERT_EQ(cluster[1][2],4);
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_clustering2)
{
    candidates c(5);

    c[0].addvalues(0,194,580,candidate::deletion);
    c[1].addvalues(1,196,185,candidate::deletion);
    c[2].addvalues(2,590,183,candidate::deletion);
    c[3].addvalues(3,797,595,candidate::deletion);
    c[4].addvalues(4,981,60,candidate::deletion);

    std::vector<std::vector<unsigned> > cluster;

    findCluster(c, cluster, candidate::insertion, 0.5);
    findCluster(c, cluster, candidate::deletion, 0.5);

    SEQAN_ASSERT_EQ(cluster.size(),2);
    SEQAN_ASSERT_EQ(cluster[0].size(),3);
    SEQAN_ASSERT_EQ(cluster[1].size(),2);
    SEQAN_ASSERT_EQ(cluster[0][0],0);
    SEQAN_ASSERT_EQ(cluster[0][1],1);
    SEQAN_ASSERT_EQ(cluster[0][2],2);
    SEQAN_ASSERT_EQ(cluster[1][0],3);
    SEQAN_ASSERT_EQ(cluster[1][1],4);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_median);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_input);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_input2);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_input3);
	//SEQAN_CALL_TEST(test_my_app_PEMer_Lite_devide);
    SEQAN_CALL_TEST(test_my_app_PEMer_Lite_clustering1);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_clustering2);
}
SEQAN_END_TESTSUITE
