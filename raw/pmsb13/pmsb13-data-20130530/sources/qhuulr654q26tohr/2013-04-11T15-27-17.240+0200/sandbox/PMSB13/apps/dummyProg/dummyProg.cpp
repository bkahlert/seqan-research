#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

template<typename TIterator>
int recGoDown(TIterator it){
    if (repLength(it)<=3){
	cout << representative(it) << " ";
	while(goDown(it) && repLength(it)<=3)
	    cout << representative(it) << " ";
    }
    while(!goRight(it)){
	if (!goUp(it))
	    return 0;
    }
    return recGoDown(it);
}

template < typename TIndex >
int myIndexExample(){
    String<char> haystack = "tobeornottobe";
    
    TIndex index(haystack);
    
    typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIteratorIndex;
    TIteratorIndex it(index);
    
    recGoDown(it);
    cout << endl;
    return 0;
}

int main()
{
    myIndexExample<Index<String<char>, IndexEsa<> > >();
    myIndexExample<Index<String<char>, IndexWotd<> > >();
    return 0;
}