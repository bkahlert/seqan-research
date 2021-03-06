#ifndef SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_
#define SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include </home/development/seqan-trunk/sandbox/PMSB13/include/seqan/dummyProg.h>

using namespace seqan;
using namespace std;

typedef Iterator<Index<String<char>, IndexEsa<> >, TopDown<ParentLinks<> > >::Type TIteratorIndex;


// A test for strings.
SEQAN_DEFINE_TEST(test_dummyProg_strings_example1)
{
    String<char> haystack = "mississippi";
    Index<String<char>, IndexEsa<> > index(haystack);
    TIteratorIndex it(index);
    
    std::stringstream redirectStream;
    std::cout.rdbuf( redirectStream.rdbuf() );
    
    std::cout << "Hello1\n";
    std::cout << "Hello2\n";
    
    string str;
    CharString res;
    while(std::getline(redirectStream, str))
    {
	cout << str;
	append(res,str);
	append(res,'\n');
    }
    
    recGoDown(it);
    SEQAN_ASSERT(res==" i p pi ppi s si ssi \n");
}

#endif  // SANDBOX_PMSB13_TESTS_DUMMYPROG_TEST_DUMMYPROG_H_
