#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <iostream>
//using namespace seqan;
int main(int, char **) {

	String<Char> str1 = "hello world";
	String<Char> str2 = "banana";
	String<Char> str3 = "mississippi";

	countOneMers(str1);
	countOneMers(str2);
	countOneMers(str3);

    return 1;
}

void countOneMers(str) {

	std::cout << "word: " << str << std::endl;
	String<int> letterNum;
	
	for(int i=0; i<length(str);i++) {
		char thisLetter = str[i];
		letterNum[thisLetter]++;
	}
	for(int i=0; i<length(letterNum);i++) {
		int thisNum = letterNum[i];
		if(thisNum>0){
			std::cout << "letter " << (char)i << ": " << thisNum << std::endl;
		}
	}
	

}