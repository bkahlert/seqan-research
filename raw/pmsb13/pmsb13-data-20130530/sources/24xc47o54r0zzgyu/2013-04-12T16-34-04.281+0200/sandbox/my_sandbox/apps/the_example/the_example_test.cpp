#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1 //unbedingt am Anfang vor allen Inlcudes

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h> 
//Hier wir eigentlich noch das zu testende Modul included

double quadrat(double x)	//Zu Testende Funkion(nen); Eigentlich in eigener datei
{
		return x*x;
}

void iota(seqan::String<int> & result, int begin, int end)
{
	resize(result, end - begin, 0);
	for (int i = begin, k=0; i<end; ++k,++i)
		result[k]=i;
}


SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_3)	//Test definieren
{
    double x = quadrat(3.0);
	SEQAN_ASSERT_IN_DELTA(x, 9.0,0.01); //IN_DELTA wegen eventueller Ungenauigkeiten (0.01 Abweichung)
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_12)	//anderen Test definieren
{
    double x = quadrat(1.2);
	SEQAN_ASSERT_IN_DELTA(x, 1.44, 0.01);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_quadrat_iota)	//anderen Test definieren
{
	seqan::String<int> result;
	iota(result, 0,3);

	SEQAN_ASSERT_EQ(length(result), 3);
	SEQAN_ASSERT_EQ(result[0],0);
	SEQAN_ASSERT_EQ(resulst[2],2);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_funcs_quadrat_3);	//Test aufrufen
	SEQAN_CALL_TEST(test_my_app_funcs_quadrat_12);
	SEQAN_CALL_TEST(test_my_app_funcs_quadrat_iota);
}
SEQAN_END_TESTSUITE