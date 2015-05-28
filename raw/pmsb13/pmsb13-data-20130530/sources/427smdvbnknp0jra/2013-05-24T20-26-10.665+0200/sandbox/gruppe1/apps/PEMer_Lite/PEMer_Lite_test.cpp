/*
Name: PEMer_Lite_test
Author: Lars Zerbe <larszerbe@live.de>
Author: Stephan Peter <s.peter@fu-berlin.de>
License: GPL v3
Copyright: 2008-2013, FU Berlin
Status: under development
*/

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include "PEMer_Lite.h"

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input)
{
	//int saminput(candidates &save, seqan::String<seqan::CharString> &header, modifier &options, statisticals &stats);

	candidates test;
	modifier opt;
	statisticals stats;
	seqan::String<seqan::CharString> header;
	char *file ="/Informatik/Development/test.sam";

	//SEQAN_ASSERT(!output(test,header,opt,stats));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input2)
{
	candidates test;
	modifier opt;
	statisticals stats;
	seqan::String<seqan::CharString> header;
	//SEQAN_ASSERT(!output(test,header,opt,stats));
	char *file =" ";

	//SEQAN_ASSERT(!output(test,header,opt,stats));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input3)
{
	candidates test;
	modifier opt;
	statisticals stats;
	seqan::String<seqan::CharString> header;
	char *file ="2579816ß15362466183b513ü6ß148z5!§%%!&°%(!%§(%§c465cn4ü457346511425129562e85z148z5!§%%!&°%(!%§(%§&&&&&(";

	//SEQAN_ASSERT(!output(test,header,opt,stats));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_sort)
{
	//int samSort(candidates &save);
	candidates test(100);
	for(unsigned i=0;i<100;i++)
	{
		test[i].addvalues(0,100-i,0,candidate::undefine);
	}
	samSort(test);
	for(unsigned i=0;i<100;i++)
	{
		SEQAN_ASSERT_EQ(test[i].pos,i+1);
	}

}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_coverage)
{
    SEQAN_FAIL("number of expected arguments: 5 != 3");
	//SEQAN_ASSERT_EQ(coverage(100,100,10000),-1.718282);
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_SD)
{
	//unsigned generateStandardDegression(candidates &input,int e_value);

	candidates test(100);
	for(unsigned i=0;i<100;i++)
	{
		test[i].addvalues(0,100,0,candidate::undefine);
	}
	generateStandardDegression(test,100);

	SEQAN_ASSERT_EQ(generateStandardDegression(test,100),0);
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_devide)
{
	//int devide(candidates &input,candidates &result,modifier &options,statisticals &stats);

	candidates test(11),result(4);
	modifier opt;
	statisticals stats;
	stats.e_value=5;
	stats.standardDegression=2;

	for(unsigned i=1;i<11;i++)
	{
		test[i].addvalues(i,i,i,candidate::undefine);
	}

	devide(test,result,opt,stats);

	SEQAN_ASSERT_EQ(result[0].type,1);
	SEQAN_ASSERT_EQ(result[0].ref,1);
	SEQAN_ASSERT_EQ(result[1].type,1);
	SEQAN_ASSERT_EQ(result[1].ref,2);
	SEQAN_ASSERT_EQ(result[2].type,2);
	SEQAN_ASSERT_EQ(result[2].ref,8);
	SEQAN_ASSERT_EQ(result[3].type,2);
	SEQAN_ASSERT_EQ(result[3].ref,9);
	SEQAN_ASSERT_EQ(result[4].type,2);
	SEQAN_ASSERT_EQ(result[4].ref,10);
}

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
	statisticals stats;
    candidates c(5);

    c[0].addvalues(0,159,417,candidate::insertion);
    c[1].addvalues(1,206,446,candidate::insertion);
    c[2].addvalues(2,836,324,candidate::deletion);
    c[3].addvalues(3,846,552,candidate::deletion);
    c[4].addvalues(4,800,477,candidate::deletion);

    std::vector<std::vector<unsigned> > cluster;

    findCluster(c, cluster, candidate::insertion, stats);
    findCluster(c, cluster, candidate::deletion, stats);

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
    statisticals stats;
    candidates c(5);

    c[0].addvalues(0,194,580,candidate::deletion);
    c[1].addvalues(1,196,185,candidate::deletion);
    c[2].addvalues(2,590,183,candidate::deletion);
    c[3].addvalues(3,797,595,candidate::deletion);
    c[4].addvalues(4,981,60,candidate::deletion);

    std::vector<std::vector<unsigned> > cluster;

    findCluster(c, cluster, candidate::insertion, stats);
    findCluster(c, cluster, candidate::deletion, stats);

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
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_sort);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_SD);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_coverage);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_devide);
    SEQAN_CALL_TEST(test_my_app_PEMer_Lite_clustering1);
	SEQAN_CALL_TEST(test_my_app_PEMer_Lite_clustering2);
}
SEQAN_END_TESTSUITE
