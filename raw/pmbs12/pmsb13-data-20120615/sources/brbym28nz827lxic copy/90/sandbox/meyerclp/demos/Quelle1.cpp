#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/find_motif.h>
#include <seqan/file.h>
#include <iostream>

using namespace seqan;



template<typename TAlphabet>
void countAllLetters(TAlphabet const&, String<AminoAcid> str){
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	unsigned int a=0;
	String<unsigned int> counter;// besser tSize
	resize(counter,alphSize,0);
	Iterator<String<AminoAcid> >::Type StringIterator = begin(str);
	Iterator<String<AminoAcid> >::Type it2 = end(str);
	while(StringIterator != it2){
		a=ordValue(*StringIterator);
		++value(counter,a);
		++StringIterator;
	}
	StringIterator=begin(str);
	Iterator<String<unsigned int> >::Type counterIterator = begin(counter);
	int i =0;
	while(i<alphSize){

		std::cout<<str[i]<<" "<<counter[i]<<std::endl;
		++i;
		
	}
}



int main(){

String<AminoAcid> str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
Iterator<String<AminoAcid> >::Type StringIterator = begin(str);
Iterator<String<AminoAcid> >::Type it2 = end(str);
while(StringIterator != it2){
	if(*StringIterator == 'R')
		*StringIterator = 'A';
	++StringIterator;
}
std::cout<<str;

//String<AminoAcid> test = "MAL";

std::cout<<std::endl;


countAllLetters(AminoAcid(),str);


/*
FrequencyDistribution<AminoAcid> F;
Iterator<String<char> >::Type it1 = begin(str);

absFreqOfLettersInSeq(F,it1,it2);

Iterator<FrequencyDistribution<AminoAcid> >::Type DistributionIterator = begin(F);
Iterator<FrequencyDistribution<AminoAcid> >::Type it4 = end(F);
while(DistributionIterator != it4){
	std::cout<< *DistributionIterator;
	++DistributionIterator;

}
*/


return 0;
}