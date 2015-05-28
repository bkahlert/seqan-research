#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

template<typename TIterator>
int recGoDown(TIterator it){
    cout << representative(it) << " ";
    while(goDown(it))
	cout << representative(it) << " ";
    while(!goRight(it)){
	if (!goUp(it))
	    return 0;
	cout << representative(it) << " ";
    }
    cout << representative(it) << " ";
    return recGoDown(it);
}

int main()
{
    String<char> haystack = "tobeornottobe";
    
    typedef Index<String<char> > TIndex;
    TIndex index(haystack);
    
    typedef Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIteratorIndex;
    TIteratorIndex it(index);
    
    recGoDown(it);
    cout << endl;
	
    return 0;
}