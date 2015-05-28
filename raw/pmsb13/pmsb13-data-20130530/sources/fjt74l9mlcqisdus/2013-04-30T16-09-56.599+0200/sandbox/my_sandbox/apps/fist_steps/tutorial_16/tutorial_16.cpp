#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;


template <typename Tit>
void once(Tit it){
	do{
		if (length(representative(it))<=3)cout << representative(it) << endl;
        if (goDown(it)==0 && goRight(it)==0)
			while (goUp(it) && goRight(it)==0) ;
	}while(isRoot(it)==0);
}
/*
template <typename >
? twice (? ?,? ?){

}
*/


int main ()
{
	typedef Index<CharString> TIndex;
    TIndex index ("tobeornottobe");
	Iterator< TIndex, TopDown<ParentLinks<>> >::Type it(index);
	once(it);
	
	IndexWotd<>
  	return 0;
}