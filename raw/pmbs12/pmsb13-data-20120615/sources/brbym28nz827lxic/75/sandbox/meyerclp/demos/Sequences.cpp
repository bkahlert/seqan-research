#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/find_motif.h>
#include <seqan/file.h>
#include <iostream>

using namespace seqan;

template<typename TString>//String<char>
void countOneMers(TString const& str){
	Iterator<TString>::Type StringIterator = begin(str);
	
	Iterator<TString>::Type EndIterator = end(str);
	String<unsigned int> counter;
	resize(counter, 26,0);//26 = AlphSize
	while(StringIterator != EndIterator){


		++value(counter,(ordValue(*StringIterator)-ordValue('a')));
		++StringIterator;
	}
}

int main(){

	String<char> str = "hello world";
	countOneMers(str);


	return 0;
}