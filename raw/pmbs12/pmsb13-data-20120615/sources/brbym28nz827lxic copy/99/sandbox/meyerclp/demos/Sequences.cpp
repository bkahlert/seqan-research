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
	unsigned int normalize =ordValue('a');
	unsigned int a=0;
	while(StringIterator != EndIterator){
		a=  ordValue(*StringIterator);
		std::cout<<a-normalize<<std::endl;
		++value(counter,(a-normalize));
		++StringIterator;
	}
	StringIterator = begin(str);
	//Iterator<String<unsigned int> >::Type countIterator = begin(counter);
	int i=0;
	while(i<26){
		if(counter[i]>0)
			std::cout<<char(i+'a')<<" "<<counter[i]<<std::endl;
		++i;
		
	}
}
void replaceAs(String<char>& str){
	Iterator<String<char>,Rooted >::Type it = begin(str);
	while(!atEnd(it)){
		if( *it=='a')
			*it='X';
		++it;
	}


}

int main(){

	String<char> str = "hello Seqan!";
	char arr[12] = "hallo welt!";
	//countOneMers(str);
	//replaceAs(str);
	//std::cout<<str;

	Iterator<String<char> >::Type it = begin(arr);
	Iterator<String<char> >::Type it2 = end(arr);
	int a=int(it);
	int b=int(it2);
	std::cout<<a<<" "<<b;
	//Length --> b-a
	//Access --> it+i

	return 0;
}