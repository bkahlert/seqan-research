#include <seqan/sequence.h>
#include <seqan/index.h>

#include <../../sandbox/PMSB13/apps/dummyProg/dummyProg.h>
using namespace seqan;
using namespace std;



int main()
{
    myIndexExample<Index<String<char>, IndexEsa<> > >();
    myIndexExample<Index<String<char>, IndexWotd<> > >();
    return 0;
}