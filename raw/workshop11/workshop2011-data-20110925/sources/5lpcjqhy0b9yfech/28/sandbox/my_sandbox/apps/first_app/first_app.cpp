#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <iostream>
using namespace seqan;


void replaceAs(String<char> str){

	std::cout << str << std::endl;

	typedef Iterator<String<char> >::Type TIterator;
	for (TIterator it = begin(str); it != end(str); ++it)
	
		if(*it == 'a'){
			assignValue(it,'X');
		}
		
	}

	std::cout << str << std::endl;

}


int main(int, char **) {

	String<char> str1 = "abcdefghijklmnopqrstuvxyz";
	String<char> str2 = "Hello SeqAn!";
	String<char> str3 = "Hello Seqan!";

	replaceAs(str1);
	replaceAs(str2);
	replaceAs(str3);

    return 1;
}
