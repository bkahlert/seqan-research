#ifndef SANDBOX_PMSB13_INCLUDE_SEQAN_DUMMYPROG_DUMMYPROG_BASE_H_
#define SANDBOX_PMSB13_INCLUDE_SEQAN_DUMMYPROG_DUMMYPROG_BASE_H_

namespace seqan {
    using namespace std;
    
    template<typename TIterator>
    int recGoDown(CharString output, TIterator it){
	if (repLength(it)<=3){
	    output += representative(it) + " ";
	    while(goDown(it) && repLength(it)<=3)
		output += representative(it) + " ";
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
}  // namespace seqan

#endif  // SANDBOX_PMSB13_INCLUDE_SEQAN_DUMMYPROG_DUMMYPROG_BASE_H_
