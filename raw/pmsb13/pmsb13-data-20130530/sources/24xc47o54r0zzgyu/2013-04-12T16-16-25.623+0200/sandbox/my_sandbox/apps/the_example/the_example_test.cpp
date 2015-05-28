#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1 //unbedingt am Anfang vor allen Inlcudes

#include <seqan/basic.h>
#include <seqan/file.h>
//Hier wir eigentlich noch das zu testende Modul included

double quadrat(double x)	//Zu Testende Funkion(nen); Eigentlich in eigener datei
{
		return 9.0;
}

SEQAN_DEFINE_TEST(test_my_app_funcs_hello)	//Test definieren
{
    double x = quadrat(3.0);
	SEQAN_ASSERT_EQ(x, 9.0);
}

SEQAN_DEFINE_TEST(test_my_app_funcs_hello2)	//anderen Test definierne
{
    double x = quadrat(1.2);
	SEQAN_ASSERT_EQ(x, 1.44);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_funcs_hello);	//Test aufrufen
	SEQAN_CALL_TEST(test_my_app_funcs_hello2);
}
SEQAN_END_TESTSUITE