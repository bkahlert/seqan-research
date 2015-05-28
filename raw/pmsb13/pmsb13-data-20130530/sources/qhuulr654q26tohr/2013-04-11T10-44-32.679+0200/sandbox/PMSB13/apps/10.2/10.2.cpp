#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

//template<typename TIterator>
int recGoDown(Iterator<Index<String<char> >, TopDown<> >::Type it){
    while(goDown(it));
    cout << representative(it) << " ";
    while(!goRight(it))
	if (!goUp(it))
	    return 0;
    return recGoDown(it);
}

int main()
{
    String<char> haystack = "tobeornottobe";
    
    typedef  TIndex;
    TIndex index(haystack);
    
    typedef Iterator<TIndex, TopDown<> >::Type TIteratorIndex;
    TIteratorIndex it(index);
    
    recGoDown(it);
    cout << endl;
	
    return 0;
}