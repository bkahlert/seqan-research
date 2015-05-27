#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <iostream>
using namespace seqan;


void permutations2(int len, String<char> str){
	if(len==2){
		for(char i = 'a'; i<= 'z',i++){
			append(str, i);
			permutations2(1,str);
		}
	}else{
		// now want third letter
		for(char i = 'a'; i<= 'z',i++){
			std::cout << str << i << ","; 
		}
		std::cout << std::endl;
	}
}

void printPermutations(int len) {

	for(char i = 'a'; i<= 'z',i++){
		
		String<char> myStr;
		append(myStr,i);
		permutations2(2,myStr);

	}

}


int main(int, char **) {

	printPermutations(3);

    return 1;
}
