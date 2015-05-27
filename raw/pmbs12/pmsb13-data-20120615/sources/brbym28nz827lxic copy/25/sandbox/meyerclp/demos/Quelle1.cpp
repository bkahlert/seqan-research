#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/find_motif.h>
#include <iostream>

using namespace seqan;

int main(){

String<char> str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
Iterator<String<char> >::Type StringIterator = begin(str);
Iterator<String<char> >::Type it2 = end(str);
while(StringIterator != it2){
	if(*StringIterator == 'R')
		*StringIterator = 'A';
	++StringIterator;
}
std::cout<<str;
return 0;
}