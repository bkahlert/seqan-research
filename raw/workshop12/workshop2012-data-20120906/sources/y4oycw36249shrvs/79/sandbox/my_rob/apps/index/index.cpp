#include <iostream>
#include <seqan/index.h>

using namespace seqan;
int main (){
    typedef Index<CharString> TIndex;
    TIndex index("dsfjsdaoijfosaj");
    Iterator< Index<String<char> >, TopDown<Preorder> >::Type it(esaIndex);
    std::cout << "JJJJJ"<< std::endl;
    do{
		std::cout << representative(it) << "gggg" << std::endl;
      //  if (!goDown(it) && !goRight(it))	
	//	  while(goUp(it) && !goRight(it));
	}while (!atEnd(it));
    return 0;
}