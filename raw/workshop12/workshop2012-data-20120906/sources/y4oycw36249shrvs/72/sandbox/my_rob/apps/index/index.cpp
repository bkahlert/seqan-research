#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main (){
    typedef Index<CharString> TIndex;
    
    TIndex index("dsfjsdaoijfosaj");
    Index<String<char> > esaIndex(index);
    Iterator< TIndex, TopDown<PreOrder> >::Type it(index);
    std::cout << "JJJJJ"<< std::endl;
    /*
    do{
		std::cout << representative(it) << "gggg" << std::endl;
        if (!goDown(it) && !goRight(it))	
		  while(goUp(it) && !goRight(it));
	}while (!isRoot(it));
	*/
    return 0;
}