#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main (){
    typedef Index<CharString> TIndex;
    TIndex index("dsfjsdaoijfosaj");
    Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);
    std::cout << "JJJJJ"<< std::endl;
    while (!isRoot(it)){
		std::cout << representative(it) << "gggg" << std::endl;
        if (!goDown(it) && !goRight(it))	
		  while(goUp(it) && !goRight(it));
	};
    return 0;
}