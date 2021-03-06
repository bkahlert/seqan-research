#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/find_motif.h>
#include <seqan/file.h>
#include <iostream>

using namespace seqan;



template<typename TAlphabet>
void countAllLetters(TAlphabet const&, String<AmicoAcid> str){
	//typedef typename Size<TAlphabet>::Type TSize;
	//TSize alphSize = ValueSize<TAlphabet>::VALUE;
	
	Iterator<String<AmicoAcid> >::Type StringIterator = begin(str);
	Iterator<String<AmicoAcid> >::Type it2 = end(str);
	while(StringIterator != it2){
		std::cout<<ordValue(value(StringIterator))<< " ";
		++StringIterator;
	}


}



int main(){

String<char> str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
//Iterator<String<AminoAcid> >::Type StringIterator = begin(str);
//Iterator<String<AminoAcid> >::Type it2 = end(str);
//while(StringIterator != it2){
//	if(*StringIterator == 'R')
//		*StringIterator = 'A';
//	++StringIterator;
//}
//std::cout<<str;

String<AminoAcid> test = "MAL";




countAllLetters(AminoAcid(),test);


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