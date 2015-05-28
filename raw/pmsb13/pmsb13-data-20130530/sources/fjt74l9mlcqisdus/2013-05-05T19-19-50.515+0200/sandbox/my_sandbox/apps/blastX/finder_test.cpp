// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1
#include "own_functions.h"


// TEST APPEND_TO_MATCH_FOUND ------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_append_to_match_found)
{
	
}
// TEST APPEND_TO_MATCH_FOUND ------------------------------------------------------------------------------------------

// TEST APPEND_TO_MATCH_FOUND ------------------------------------------------------------------------------------------
SEQAN_DEFINE_TEST(test_my_app_find_matches)
{
	
}
// TEST APPEND_TO_MATCH_FOUND ------------------------------------------------------------------------------------------

// TEST DURCHLAUF ------------------------------------------------------------------------------------------------------
// Normale eingabe
SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(test_my_app_append_to_match_found);
	SEQAN_CALL_TEST(test_my_app_find_matches);
	
}
SEQAN_END_TESTSUITE
// TEST DURCHLAUF ------------------------------------------------------------------------------------------------------
