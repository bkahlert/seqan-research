#include </home/development/seqan-trunk/sandbox/PMSB13/include/seqan/dummyProg.h>
using namespace seqan;
using namespace std;



int main()
{
    myIndexExample<Index<String<char>, IndexEsa<> > >();
    myIndexExample<Index<String<char>, IndexWotd<> > >();
    return 0;
}