#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <iostream>
using namespace seqan;

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

	String<char> str1 = "hello world";
	String<Dna> str2 = "TATACGCTA";
	String<AminoAcid> str3 = "MQDRVKRPMNAFIVWSRDQRRKMALEN";

	countOneMers(str1);
	countOneMers(str2);
	countOneMers(str3);

    return 1;
}
