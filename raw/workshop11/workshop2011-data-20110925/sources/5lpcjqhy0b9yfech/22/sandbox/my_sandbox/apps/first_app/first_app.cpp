#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <iostream>
using namespace seqan;


void printPermutations(int len) {

	SimpleType<char, Finite<len> > alphPart;
	std::cout << "chosen alphabet: " << alphPart << std::endl;

}


template <typename T> 
void countOneMers(T & str) {

	std::cout << "word: " << str << std::endl;
	String<int> letterNum;
	resize(letterNum, 256, 0);
	
	for(int i=0; i<length(str);i++) {
		typename Value<T>::Type thisLetter = str[i];
		letterNum[thisLetter]++;
	}

	for(int i=0; i<length(letterNum);i++) {
		int thisNum = letterNum[i];
		if(thisNum>0){
			std::cout << "letter " << (typename Value<T>::Type)i << ": " << thisNum << std::endl;
		}
	}
	

}

int main(int, char **) {

	printPermutations(3);

    return 1;
}
