#ifndef SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_
#define SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include </home/development/seqan-trunk/sandbox/PMSB13/include/seqan/dummyProg.h>

using namespace seqan;
using namespace std;

typedef Iterator<Index<String<char>, IndexEsa<> >, TopDown<ParentLinks<> > >::Type TIteratorIndex;


TIteratorIndex getIt(){
    String<char> haystack = "mississippi";
    TIndex index(haystack);
    TIteratorIndex it(index);
    return it;
}

// A test for strings.
SEQAN_DEFINE_TEST(test_dummyProg_strings_example1)
{
    TIteratorIndex it = getIt();
    
    std::stringstream redirectStream;
    std::cout.rdbuf( redirectStream.rdbuf() );
    
    std::cout << "Hello1\n";
    std::cout << "Hello2\n";
    
    while(std::getline(redirectStream, str))
    {
	// This does not work - as the contents of redirectStream 
	// do not include the '\n' - I only see "Hello1Hello2"
    }
    
    //recGoDown(it);
    SEQAN_ASSERT(true);
}

#endif  // SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_
